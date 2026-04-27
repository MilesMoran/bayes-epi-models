data
{
    int<lower=1> T,H,G,I;
    matrix<upper=0>[G,I] log_e;
    array[T,G,I] int<lower=0> y;
    array[H,G,I] int<lower=0> y_test;
    vector<lower=0,upper=1>[T] xmas_train;
    vector<lower=-1,upper=1>[T] sin_omega_train, cos_omega_train;
    vector<lower=0,upper=1>[H] xmas_test;
    vector<lower=-1,upper=1>[H] sin_omega_test, cos_omega_test;
}
transformed data
{
    int N = G*I;
    int TN = T*G*I;
    array[TN] int y_arr = to_array_1d(y);
}
parameters
{
    real end_1;
    real end_christmas;
    row_vector[I-1] end_strata_raw; // alpha_strata[i]
    vector[G-1] end_geogs_raw;  // alpha_geogs[g]

    row_vector[I] end_sin;          // gamma[i]
    row_vector[I] end_cos;          // delta[i]

    real nlog_overdisp;  // -log(psi)
}
transformed parameters 
{
    real overdisp = exp(-nlog_overdisp);
    real inv_overdisp = exp(nlog_overdisp);
}
model
{
    end_1 ~ normal(3,2);
    end_christmas  ~ normal(0,3); // beta
    end_strata_raw ~ normal(0,3); // alpha_strata[i]
    end_geogs_raw  ~ normal(0,3); // alpha_geogs[g]

    end_sin ~ normal(0,3);        // gamma[i]
    end_cos ~ normal(0,3);        // delta[i]

    nlog_overdisp ~ normal(-0.5,1); // -log(psi)

    {
        row_vector[I] end_strata = append_col(0, end_strata_raw);
        vector[G] end_geogs  = append_row(0, end_geogs_raw);

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

        y_arr ~ neg_binomial_2(to_array_1d(mu_end'), inv_overdisp);
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
        for(g in 1:G) {
        for(i in 1:I) {
            int n = (g-1)*I + i;
            real mu = mu_end[h,n];
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

        }}}
    } // end scope
}
