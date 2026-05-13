data
{
    int<lower=1> T,G,I;             // data dimensions
    matrix<upper=0>[G,I] log_e;     // PUMA-by-age population fractions
    matrix[G,G] D;                  // adjacency-order "distance" between PUMAs
    array[T,G,I] int<lower=0> y;    // incident counts (assumed eq. to prevalence)
    vector<lower=0,upper=1>[T] xmas; // indicator for week of christmas shock
    vector<lower=-1,upper=1>[T] sin_omega, cos_omega; // endemic seasonality components
    
    // note: for the W_strata used in the Berlin Norovirus analysis, all eigen-parts were 
    //       real-valued. If your eigen-parts are complex, this simply won't work!
    vector[I] EValC; // eigenvalues of (pre-row-normalized) contact matrix C
    matrix[I,I] EVecC; // eigenvectors of (pre-row-normalized) contact matrix C
}
transformed data
{
    int N = G*I;
    int TN = T*G*I;

    matrix[G,G] logD = log(D);
    matrix[I,I] Inv_EVecC = inverse(EVecC);

    array[TN] int y_arr = to_array_1d(y);
    array[T] matrix[G,I] y_mats;
    for(t in 1:T) {
        y_mats[t,,] = to_matrix(y[t,,]);
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

    real logpower;                  // age-group mixing / eigen-deformation parameter (\log() of \kappa in manuscript)
}
transformed parameters 
{
    real inv_overdisp = inv(overdisp);
    real kappa = exp(logpower);
    matrix[I,I] W_strata;

    {
        //*** Construct W_strata via eigenvalue matrix power deformation ***//
        matrix[I,I] W_strata_raw = (EVecC * diag_matrix(pow(EValC,kappa)) * Inv_EVecC)';
        for(c in 1:I) {
            for(r in 1:I) {
                if(W_strata_raw[r,c] < 0) W_strata_raw[r,c] = 0;
            }
            real colSum = sum(W_strata_raw[,c]);
            if(colSum != 1) {
                W_strata_raw[,c] = W_strata_raw[,c]/colSum;
            }
            W_strata[,c] = W_strata_raw[,c];
        }        
    }
}
model
{
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
    logpower ~ normal(0,0.75);

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

        matrix[T,N] wY_prev;
        wY_prev[1,] = rep_row_vector(0,N);
        for (t in 2:T) {
            wY_prev[t,] = to_row_vector((W_geog * (y_mats[t-1,,] * W_strata_t))');
        }

        matrix[T,N] mu_y_t = mu_end + (wY_prev .* rep_matrix(mu_ne_mult,T));

        // likelihood
        y_arr ~ neg_binomial_2(to_array_1d(mu_y_t'), inv_overdisp);
    }
}