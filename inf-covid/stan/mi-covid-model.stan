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
    array[] int nonzero(array[] int x) {
        int n_el = size(x);
        array[n_el + 1] int idx;
        int dummy_i = 0;
        for (i in 1:n_el) {
            if (x[i] != 0) {
                dummy_i += 1;
                idx[1 + dummy_i] = i;
            }
        }
        idx[1] = dummy_i;
        return(idx);
    }
    real partial_sum(array[] matrix r_slice,
                     int start, int end,
                     data array[,] int dy,
                     row_vector delta,
                     matrix W_geog,
                     matrix W_strata_t,
                     real k,
                     data vector E_inv,
                     data array[,] int x_hat,
                     data int G, data int I) {
        int n_t = end - start + 1;
        int N = size(delta);

        matrix[n_t,N] rW;
        for (t in 1:n_t) {
            // Since I<<G in most instances, compute r[t,,]*W_strata first
            matrix[G,I] rW_t = W_geog * (r_slice[t,,] * W_strata_t);
            for (g in 1:G) {
            for (i in 1:I) {
                rW[t, (g-1)*I + i] = rW_t[g,i];
            }}
        }
        matrix[N, n_t] k_p = k * exp(-(rep_matrix(delta, n_t) + diag_post_multiply(rW, E_inv))');
        // likelihood
        return(beta_binomial_lpmf(to_array_1d(dy[start:end]) | to_array_1d(x_hat[start:end]),
                                                          to_array_1d(k - k_p),
                                                          to_array_1d(k_p)));
    }
}
data 
{
    int<lower=1> T,G,I;           // data dimensions
    array[G,I] int<lower=1> E;    // PUMA-by-age population counts
    matrix[G,G] D;                // "human distance" between PUMAs
    array[T,G,I] int<lower=0> dy; // incident counts
    int<lower=1> grainsize;       // (only use if parallelizing this script)
}
transformed data 
{
    int N = G*I;
    int TN = T*G*I;

    array[T,N] int x_hat; // estimate of susceptible pool size in SIR-like trajectory
    array[T,N] int dy_arr;
    
    // because r[t,g,i] is degenerately equal to 0 whenever Y_hat[t,g,i]=0,
    // we have to pre-determine which Y_hat[t,g,i] will end-up zero or nonzero.
    // Here, `n_nonzero_y` is the number of nonzero Y_hat estimates. The indices
    // for those nonzero elements flattened and stored as `idx_t_raw[idx]`, 
    // `idx_n_raw[idx]`, etc. and then used to sample ONLY those specific elements
    // of the r array that are non-degenerate. All other elements of `r` are 
    // set equal to 0 manually.
    
    array[TN] int idx_t_raw, idx_n_raw, idx_g_raw, idx_i_raw;
    int n_nonzero_y = 0;
    x_hat[1,] = to_array_1d(E); // array[,] -> array[] conversion is always Row-major
    for (g in 1:G)
    for (i in 1:I)
       dy_arr[1,(g-1)*I + i] = dy[1,g,i];
    for(t in 2:T) {
    for(g in 1:G) {
    for(i in 1:I) {
        int n = (g-1)*I + i; // flatten in ROW-major order
        int cumulative_incidence = sum( dy[1:(t-1),g,i] );
        x_hat[t,n] = (E[g,i] - cumulative_incidence);
        dy_arr[t,n] = dy[t,g,i];

        if(cumulative_incidence > 0) {
            n_nonzero_y += 1;
            idx_t_raw[n_nonzero_y] = t;
            idx_n_raw[n_nonzero_y] = n;
            idx_g_raw[n_nonzero_y] = g;
            idx_i_raw[n_nonzero_y] = i;
        }
    }}}
    array[n_nonzero_y] int idx_t = idx_t_raw[1:n_nonzero_y];
    array[n_nonzero_y] int idx_n = idx_n_raw[1:n_nonzero_y];
    array[n_nonzero_y] int idx_g = idx_g_raw[1:n_nonzero_y];
    array[n_nonzero_y] int idx_i = idx_i_raw[1:n_nonzero_y];
    array[n_nonzero_y] int idx_tmo;
    for(idx_nonzero in 1:n_nonzero_y) {
        idx_tmo[idx_nonzero] = (idx_t[idx_nonzero] - 1);
    }

    vector[T-1] TMinusTwoToZero;
    for(t in 1:(T-1)) { TMinusTwoToZero[t] = (T - t - 1); }

    // Most of this stuff (e.g. log_dy_mat) is just precomputation of static 
    // elements so that FLOPs involving only the data are not performed for
    // each parameter draw

    array[N,T] int  dy_mat;
    matrix[T,N] log_dy_mat;
    array[N,T] int n_nonzero_by_n;
    array[N, T, T] int idx_nonzero_by_n;
    vector[N] E_inv;
    for (g in 1:G) {
    for (i in 1:I) {
        int n = (g-1)*I + i; // flatten in ROW-major order
        E_inv[n] =  inv(E[g,i]);
        dy_mat[n,] = dy[,g,i];
        log_dy_mat[,n] = log(to_vector(dy_mat[n,]));
        for (t in 2:T) {
          array[t] int nonzero_idx = nonzero(dy_mat[n,1:(t-1)]);
          int n_nonzero_i = nonzero_idx[1];
          n_nonzero_by_n[n,t] = n_nonzero_i;
          if (n_nonzero_i > 0) {
               idx_nonzero_by_n[n,t,1:n_nonzero_i] = nonzero_idx[2:(n_nonzero_i + 1)];
          }
        }
    }}

    array[T*N] int dy_1d_arr = to_array_1d(dy);
    array[T*N] int x_hat_1d_arr = to_array_1d(x_hat);

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
    row_vector<lower=0>[N] delta;  // endemic addend in linear predictor
    real<lower=0> theta;           // latent-level dispersion param (gamma scale)
    real<lower=0> k;               // obsv.-level dispersion param (beta-binom precision)
    real<lower=0> gamma;           // recovery rate
    vector[T-1] log_R;             // latent mean scalar, akin to effective repro. number
    vector[n_nonzero_y] log_r_raw; // noncentered parameterization for the ~LogNormal approx.
                                   // of r[t,g,i], which is assumed ~Gamma in reality
    vector[(I*(I-1))%/%2] C_comp;  // params that form column-compositions of C
    vector<lower=0>[I] u;          // params that form column-scales of C
    vector[G-1] eta_log_tau;       // spatial weights destination-attraction parameter
    vector[G] eta_log_rho;         // ncp "raw" special weights distance-decay parameters
    real mean_log_rho;             // mean of rho random effects
    real<lower=0> sd_log_rho;      // sd of rho random effects
}
transformed parameters 
{
    matrix[G,G] W_geog;
    matrix[I,I] W_strata;
    vector[I] log_alpha;
    vector[G] rho = exp(mean_log_rho + sd_log_rho * eta_log_rho);
    vector[G] log_tau = append_row(eta_log_tau,[0]');

    {
        //*** Construct W_geog ***//
        matrix[G,G] D_rho = diag_post_multiply(D, rho);
        for (g_ in 1:G) { // g_ = col = origin
            W_geog[,g_] = softmax(log_tau - D_rho[,g_]);
        }
        
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
        log_alpha[1] = log(u[1]);
        
        for (i in 2:I) {
            real alpha_i = sum(C[1:(i-1),i]) + u[i];
            W_strata[,i] = C[,i] / alpha_i;
            log_alpha[i] = log(alpha_i);
        }
    }
}
model 
{
    matrix[I,I] W_strata_t = W_strata';
    log_R ~ normal(0,2);
    theta ~ cauchy(0,5);
    k ~ cauchy(0,1e4);
    gamma ~ cauchy(0,1); 
    eta_log_tau ~ normal(0,3);
    eta_log_rho ~ std_normal();
    delta ~ std_normal();
    u ~ gamma(gamma_shapes, 2*I);
    mean_log_rho ~ normal(0,2);
    sd_log_rho ~ normal(0,1);

    log_r_raw ~ std_normal();
    array[T] matrix[G,I] r;
    {
        //*** Construct y_hat_vec[idx] ***//
        vector[T-1] log_decay_weights = -gamma * TMinusTwoToZero;
        matrix[T,N] log_y_hat_mat;
        for (t in 2:T) {
          vector[t-1] tail_weights = tail(log_decay_weights,t-1);
          for (n in 1:N) {
            int n_nonzero_i = n_nonzero_by_n[n,t];
            if (n_nonzero_i > 0) {
               array[n_nonzero_i] int nonzero_idx_i = idx_nonzero_by_n[n,t,1:n_nonzero_i];
               log_y_hat_mat[t,n] = log_sum_exp(tail_weights[nonzero_idx_i]
                                                  + (log_dy_mat[1:(t-1),n])[nonzero_idx_i]);
              }
          }
        }

        vector[n_nonzero_y] y_hat_vec;
        for (idx_nonzero in 1:n_nonzero_y) {
            int t0 = idx_t[idx_nonzero];
            int n0 = idx_n[idx_nonzero];
            y_hat_vec[idx_nonzero] = log_y_hat_mat[t0,n0];
        }

        //*** Construct r[t,g,i] ***//
        vector[n_nonzero_y] log_mu_r = log_R[idx_tmo] + log_alpha[idx_i] + y_hat_vec;
        vector[n_nonzero_y] vr_log_r = log1p(theta*exp(-log_mu_r));
        vector[n_nonzero_y] r_vec    = exp(log_mu_r + -0.5*vr_log_r + sqrt(vr_log_r).*log_r_raw);
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
    target += reduce_sum(partial_sum, r, grainsize, dy_arr, delta, W_geog, W_strata_t, k, E_inv, x_hat, G, I);
}
