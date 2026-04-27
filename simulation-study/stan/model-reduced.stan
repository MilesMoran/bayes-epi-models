data
{
    int<lower=1> T,H,G,I;
    matrix<upper=0>[G,I] log_e;
    matrix[G,G] D; 
    array[T,G,I] int<lower=0> y;
    array[H,G,I] int<lower=0> y_test;
    vector<lower=0,upper=1>[T] xmas_train;
    vector<lower=-1,upper=1>[T] sin_omega_train, cos_omega_train;
    vector<lower=0,upper=1>[H] xmas_test;
    vector<lower=-1,upper=1>[H] sin_omega_test, cos_omega_test;
    
    // note: for the W_strata used in the research study, all eigen-parts were 
    //       real-valued. If your eigen-parts are complex, this won't work!
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
    real end_1;
    real end_christmas;
    row_vector[I-1] end_strata_raw; // alpha_strata[i]
    vector[G-1] end_geogs_raw;  // alpha_geogs[g]

    row_vector[I] end_sin;          // gamma[i]
    row_vector[I] end_cos;          // delta[i]
    
    real ne_1;
    row_vector[I-1] ne_strata_raw;
    vector[G-1] ne_geogs_raw;   // log(phi_geogs[g])

    real ne_log_pop;     // tau
    real neweights_logd; // log(rho)
    real nlog_overdisp;  // -log(psi)
    
    real logpower; // log(kappa)
}
transformed parameters 
{
    real overdisp = exp(-nlog_overdisp);
    real inv_overdisp = exp(nlog_overdisp);
    real kappa = exp(logpower);

    column_stochastic_matrix[I,I] W_strata;
    {
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
    end_christmas  ~ normal(0,3); // beta
    end_strata_raw ~ normal(0,3); // alpha_strata[i]
    end_geogs_raw  ~ normal(0,3); // alpha_geogs[g]

    end_sin ~ normal(0,3);        // gamma[i]
    end_cos ~ normal(0,3);        // delta[i]
    ne_1 ~ normal(2,5);
    ne_strata_raw ~ normal(0,3);  // log(phi_strata[i])
    ne_geogs_raw  ~ normal(0,3);  // log(phi_geogs[g])

    ne_log_pop ~ normal(0,2);     // tau
    neweights_logd ~ std_normal(); // log(rho)
    nlog_overdisp ~ normal(-0.5,1); // -log(psi)
    logpower ~ normal(0,0.75); // log(kappa)

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
        vector[T]   mu_end_xmas = (end_christmas * xmas_train);
        matrix[T,I] mu_end_ssnl = ((sin_omega_train * end_sin) + 
                                   (cos_omega_train * end_cos));
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

        y_arr ~ neg_binomial_2(to_array_1d(mu_y_t'), inv_overdisp);
    }
}
generated quantities
{
    array[H,G,I] int y_pred;
    array[H] matrix[G,I] logscore_y_test;
    array[H] matrix[G,I] logscore_y_test_O;
    array[H] matrix[G,I] logscore_y_test_C;
    array[H] matrix[G,I] logscore_y_pred;
    array[H] matrix[G,I] logscore_y_pred_O;
    array[H] matrix[G,I] logscore_y_pred_C;
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
        vector[H]   mu_end_xmas = (end_christmas * xmas_test);
        matrix[H,I] mu_end_ssnl = ((sin_omega_test * end_sin) +
                                   (cos_omega_test * end_cos));
        matrix[H,N] mu_end_ssnl_extended;
        for(g in 1:G) {
            mu_end_ssnl_extended[,((g-1)*I+1):(g*I)] = mu_end_ssnl;
        }
        matrix[H,N] mu_end = exp(end_1 + rep_matrix(mu_end_base,H) +
                                         rep_matrix(mu_end_xmas,N) +
                                         mu_end_ssnl_extended);

        for(h in 1:H) {

            matrix[G,I] y_hmo = ((h==1) ? y_mats[T,,] : to_matrix(y_pred[h-1,,]));
            row_vector[N] wY_hmo = to_row_vector((W_geog * (y_hmo * W_strata_t))');
            row_vector[N] mu_y_h = mu_end[h,] + (wY_hmo .* mu_ne_mult);
            
            for(g in 1:G) {
            for(i in 1:I) {
                int n = (g-1)*I + i;
                real mu = mu_y_h[n];
                real log_p_z = neg_binomial_2_lpmf(0 | mu, inv_overdisp); // log(Pr(y[h,g,i] = 0 | draw))
                real log_p_nz = log1m_exp(log_p_z); // log(Pr(y[h,g,i] > 0 | draw))
    
                logscore_y_test[h,g,i] = neg_binomial_2_lpmf(y_test[h,g,i] | mu, inv_overdisp); 
                if (y_test[h,g,i] > 0) {
                    logscore_y_test_O[h,g,i] = log_p_nz;
                    logscore_y_test_C[h,g,i] = logscore_y_test[h,g,i] - log_p_nz;
                } else {
                    logscore_y_test_O[h,g,i] = log_p_z;
                    logscore_y_test_C[h,g,i] = 0;
                }
    
                int y_rng = neg_binomial_2_rng(mu, inv_overdisp);
                y_pred[h,g,i] = y_rng;
    
                logscore_y_pred[h,g,i] = neg_binomial_2_lpmf(y_pred[h,g,i] | mu, inv_overdisp); 
                if (y_pred[h,g,i] > 0) {
                    logscore_y_pred_O[h,g,i] = log_p_nz;
                    logscore_y_pred_C[h,g,i] = logscore_y_pred[h,g,i] - log_p_nz;
                } else {
                    logscore_y_pred_O[h,g,i] = log_p_z;
                    logscore_y_pred_C[h,g,i] = 0;
                }
    
            }}
        } // end h-loop
    } // end scope
}
