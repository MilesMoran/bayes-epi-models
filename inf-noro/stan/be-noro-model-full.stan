functions 
{
    vector inv_ilr_dirichlet_prior_lp(vector y, vector shapes) {
        int N = rows(y) + 1;
        vector[N - 1] ns = linspaced_vector(N - 1, 1, N - 1);
        vector[N - 1] w = y ./ sqrt(ns .* (ns + 1));
        vector[N] z = append_row(reverse(cumulative_sum(reverse(w))), 0) - append_row(0, ns .* w);
        real r = log_sum_exp(z);
        vector[N] x = exp(z - r);
        target += 0.5 * log(N);
        target += sum(z) - N * r;
        target += dirichlet_lpdf(x | shapes[1:N]);
        return x;
    }
}
data
{
    int<lower=1> T,G,I;             // data dimensions
    matrix<upper=0>[G,I] log_e;     // PUMA-by-age population fractions
    matrix[G,G] D;                  // adjacency-order "distance" between PUMAs
    array[T,G,I] int<lower=0> y;    // incident counts (assumed eq. to prevalence)
    vector<lower=0,upper=1>[T] xmas; // indicator for week of christmas shock
    vector<lower=-1,upper=1>[T] sin_omega, cos_omega; // endemic seasonality components
}
transformed data
{
    int N = G*I;
    int TN = T*G*I;

    matrix[G,G] logD = log(D);
    array[TN] int y_arr = to_array_1d(y);

    // because r[t,g,i] is degenerately equal to 0 whenever y[t-1,g,i]=0,
    // we have to pre-determine which y[t,g,i] will end-up zero or nonzero.
    // Here, `n_nonzero_y` is the number of nonzero y's. The indices
    // for those nonzero elements flattened and stored as `idx_t_raw[idx]`, 
    // `idx_n_raw[idx]`, etc. and then used to sample ONLY those specific elements
    // of the r array that are non-degenerate. All other elements of `r` are 
    // set equal to 0 manually.
    
    array[TN] int idx_t_raw, idx_g_raw, idx_i_raw;
    int n_nonzero_y = 0;
    for(t in 2:T) {
    for(g in 1:G) {
    for(i in 1:I) {
        if(y[(t-1),g,i] > 0) {
            n_nonzero_y += 1;
            idx_t_raw[n_nonzero_y] = t;
            idx_g_raw[n_nonzero_y] = g;
            idx_i_raw[n_nonzero_y] = i;
        }
    }}}
    array[n_nonzero_y] int idx_t = idx_t_raw[1:n_nonzero_y];
    array[n_nonzero_y] int idx_g = idx_g_raw[1:n_nonzero_y];
    array[n_nonzero_y] int idx_i = idx_i_raw[1:n_nonzero_y];

    vector[n_nonzero_y] log_y_prev_vec;
    for (idx_nonzero in 1:n_nonzero_y) {
        int t0 = idx_t[idx_nonzero];
        int g0 = idx_g[idx_nonzero];
        int i0 = idx_i[idx_nonzero];
        log_y_prev_vec[idx_nonzero] = log(y[(t0-1),g0,i0]);
    }

    //************************************************************************//    
    //*** quantities for the prior on W_strata ***// 

    array[I] int strata_dof;
    for (i in 1:I) {
        strata_dof[i] = I - i;
    }
    array[I] int strata_start;
    strata_start[1] = 1;
    for (i in 2:I) {
        strata_start[i] = strata_start[i-1] + strata_dof[i-1];
    }

    real w_diag_prop = 0.4; // expected value of W_strata[i,i]'s 
    real concentration = 1.8;
    real shape_diag = I * w_diag_prop * concentration;
    real shape_offd = (1.0*I/(I-1)) * (1-w_diag_prop) * concentration;
    vector[I] shapes = append_row(shape_diag, rep_vector(shape_offd, I-1));
    vector[I] gamma_shapes;
    for (i in 1:I) {
        gamma_shapes[i] = shape_diag + shape_offd*(I-i);
    }
}
parameters
{
    real end_1;                     // endemic reference group baseline  (\beta_0         in manuscript)
    real end_christmas;             // endemic christmas shock effect    (\beta^{xmas}    in manuscript)
    row_vector[I-1] end_strata_raw; // endemic per-strata additive shift (\beta^{(I)}_{i} in manuscript)
    vector[G-1] end_geogs_raw;      // endemic per-geog additive shift   (\beta^{(G)}_{g} in manuscript)
    row_vector[I] end_sin;          // endemic per-strata sine effect    (\beta^{(S)}_{i} in manuscript)
    row_vector[I] end_cos;          // endemic per-strata cosine effect  (\beta^{(C)}_{i} in manuscript)
    
    real ne_log_pop;                // epidemic population-size effect    (\eta^{(pop)}   in manuscript)
    real ne_1;                      // epidemic reference group baseline  (\eta_0         in manuscript)
    row_vector[I-1] ne_strata_raw;  // epidemic per-strata additive shift (\eta^{(I)}_{i} in manuscript)
    vector[G-1] ne_geogs_raw;       // epidemic per-geog. additive shift  (\eta^{(G)}_{g} in manuscript)

    real neweights_logd;            // geog. mixing weights distance decay rate (\log() of \rho in manuscript)
    real<lower=0> overdisp;         // obsv.-level dispersion param (negbin. inv-size, \psi in manuscript)
    real<lower=0> theta;            // latent-level dispersion param (gamma scale, \theta in manuscript)

    vector[n_nonzero_y] log_r_raw;  // noncentered parameterization for the ~LogNormal approx.
                                    // of r[t,g,i], which is assumed ~Gamma in reality
    vector[(I*(I-1))%/%2] C_comp;   // params that form column-compositions of C
    vector<lower=0>[I] u;           // params that form column-scales of C
}
transformed parameters 
{
    real inv_overdisp = inv(overdisp);
    matrix[I,I] W_strata;

    {
        //*** Construct W_strata ***//
        matrix[I,I] C;
        for (i in 1:(I-1)) {
            int lo = strata_start[i];
            int hi = (strata_start[i+1] - 1);
            C[i:I,i] = u[i] * inv_ilr_dirichlet_prior_lp(C_comp[lo:hi], shapes);
            for (j in 1:(i-1))
                C[j,i] = C[i,j];
        }
        for (j in 1:(I-1)) {
            C[j,I] = C[I,j];
        }
        C[I,I] = u[I];        
        
        W_strata[,1] = C[,1] / u[1]; // u[1] is the 1st column scale; u[i] scale subsequent sub-columns
        for (i in 2:I) {
            W_strata[,i] = C[,i] / (sum(C[1:(i-1),i]) + u[i]);
        }
    }
}
model
{
    u ~ gamma(gamma_shapes, 2*I);
    
    end_1 ~ normal(3,2);
    end_christmas  ~ normal(0,3);
    end_strata_raw ~ normal(0,3);
    end_geogs_raw  ~ normal(0,3);

    end_sin ~ normal(0,3);
    end_cos ~ normal(0,3);
    ne_1 ~ normal(2,5);
    ne_strata_raw ~ normal(0,3);
    ne_geogs_raw  ~ normal(0,3);

    ne_log_pop ~ normal(0,2);  
    neweights_logd ~ std_normal(); 
    overdisp ~ cauchy(0,1);
    theta ~ cauchy(0,1);
    
    log_r_raw ~ std_normal();
    array[T] matrix[G,I] r;
    {
        //*** Construct r[t,g,i] ***//
        vector[n_nonzero_y] vr_log_r = log1p(theta*exp(-log_y_prev_vec));
        vector[n_nonzero_y] r_vec    = exp(log_y_prev_vec + -0.5*vr_log_r + sqrt(vr_log_r).*log_r_raw);
        for (t in 1:T) {
            r[t,,] = rep_matrix(0,G,I);
        }
        for (idx_non_zero in 1:n_nonzero_y) {
            int t0 = idx_t[idx_non_zero];
            int g0 = idx_g[idx_non_zero];
            int i0 = idx_i[idx_non_zero];
            r[t0,g0,i0] = r_vec[idx_non_zero];
        }
    }
    {
        row_vector[I] end_strata = append_col(0, end_strata_raw);
        vector[G] end_geogs  = append_row(0, end_geogs_raw);
        row_vector[I] ne_strata = append_col(0, ne_strata_raw);
        vector[G] ne_geogs  = append_row(0, ne_geogs_raw);

        matrix[I,I] W_strata_t = W_strata';
        matrix[G,G] W_geog;
        matrix[G,G] eta = -exp(neweights_logd)*logD;
        for (g_ in 1:G) { // g_ = col = origin
            W_geog[,g_] = softmax(eta[,g_]);
        }

        row_vector[N] mu_ne_mult;
        mu_ne_mult = to_row_vector((ne_log_pop*log_e + rep_matrix(ne_geogs,I) + rep_matrix(ne_strata,G))');
        mu_ne_mult = exp(ne_1 + mu_ne_mult);
        
        row_vector[N] mu_end_base = to_row_vector((log_e + rep_matrix(end_geogs,I) + rep_matrix(end_strata,G))');
        vector[T]   mu_end_xmas = (end_christmas * xmas);
        matrix[T,I] mu_end_ssnl = ((sin_omega * end_sin) + 
                                   (cos_omega * end_cos));
        matrix[T,N] mu_end_ssnl_extended;
        for(g in 1:G) {
            mu_end_ssnl_extended[,((g-1)*I+1):(g*I)] = mu_end_ssnl;
        }
        matrix[T,N] mu_end = exp(end_1 + rep_matrix(mu_end_base,T) + 
                                         rep_matrix(mu_end_xmas,N) + 
                                         mu_end_ssnl_extended);

        matrix[T,N] wr;
        wr[1,] = rep_row_vector(0,N);
        for (t in 2:T) {
            wr[t,] = to_row_vector((W_geog * (r[t,,] * W_strata_t))');
        }

        matrix[T,N] mu_y_t = mu_end + (wr .* rep_matrix(mu_ne_mult,T));

        y_arr ~ neg_binomial_2(to_array_1d(mu_y_t'), inv_overdisp);
    }
}