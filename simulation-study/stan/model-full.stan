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
    int<lower=1> T,H,G,I;
    matrix<upper=0>[G,I] log_e;
    matrix[G,G] D; 
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

    matrix[G,G] logD = log(D);

    array[TN] int y_arr = to_array_1d(y);
    
    int n_nonzero_y = 0;
    array[TN] int idx_t_raw, idx_g_raw, idx_i_raw;
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
    for (i in 2:I)
      strata_start[i] = strata_start[i-1] + strata_dof[i-1];
      
    real w_diag_prop = 0.4;
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
    real nlog_theta;     // -log(theta)
    
    vector[n_nonzero_y] log_r_raw;
    
    // real<lower=0,upper=1> w_diag_prop;
    // real<lower=0> concentration;
    vector[(I*(I-1))%/%2] C_comp;
    vector<lower=0>[I] u;
}
transformed parameters 
{
    real overdisp = exp(-nlog_overdisp);
    real inv_overdisp = exp(nlog_overdisp);
    real theta = exp(-nlog_theta);
    
    matrix[I,I] C;
    matrix[I,I] W_strata;

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

    W_strata[,1] = C[,1] / u[1];
    for (i in 2:I)
        W_strata[,i] = C[,i] / (sum(C[1:(i-1),i]) + u[i]);
}
model
{
    // w_diag_prop ~ beta(4,6);
    // concentration ~ normal(0,1.25);
    u ~ gamma(gamma_shapes, 2*I);
    
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
    nlog_theta ~ normal(-2,1); // -log(theta)
    
    log_r_raw ~ std_normal();
    array[T] matrix[G,I] r;
    {
        //*** Construct r[t,g,i] ***//
        vector[n_nonzero_y] vr_log_r = log1p_exp(-(nlog_theta+log_y_prev_vec));
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

        matrix[T,N] wr;
        wr[1,] = rep_row_vector(0,N);
        for (t in 2:T) {
            wr[t,] = to_row_vector((W_geog * (r[t,,] * W_strata_t))');
        }

        matrix[T,N] mu_y_t = mu_end + (wr .* rep_matrix(mu_ne_mult,T));

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
            
            int nnz_y_gq = 0;
            array[N] int idx_g_raw_gq, idx_i_raw_gq;
            for(g in 1:G) {
            for(i in 1:I) {
                if(((h==1) && (y[T,g,i] > 0)) ||
                   ((h >1) && (y_pred[h-1,g,i] > 0))) {
                    nnz_y_gq += 1;
                    idx_g_raw_gq[nnz_y_gq] = g;
                    idx_i_raw_gq[nnz_y_gq] = i;
                }
            }}
            array[nnz_y_gq] int idx_g_gq = idx_g_raw_gq[1:nnz_y_gq];
            array[nnz_y_gq] int idx_i_gq = idx_i_raw_gq[1:nnz_y_gq];

            vector[nnz_y_gq] log_yp_vec_gq;
            for (idx_nonzero in 1:nnz_y_gq) {
                int g0 = idx_g_gq[idx_nonzero];
                int i0 = idx_i_gq[idx_nonzero];
                if(h==1) {
                    log_yp_vec_gq[idx_nonzero] = log(y[T,g0,i0]);
                } else {
                    log_yp_vec_gq[idx_nonzero] = log(y_pred[h-1,g0,i0]);
                }
            }

            vector[nnz_y_gq] log_r_raw_gq = to_vector(
                normal_rng(rep_vector(0,nnz_y_gq), rep_vector(1,nnz_y_gq))
            );
            vector[nnz_y_gq] vr_log_r = log1p_exp(-(nlog_theta+log_yp_vec_gq));
            vector[nnz_y_gq] r_vec    = exp(log_yp_vec_gq + -0.5*vr_log_r + sqrt(vr_log_r).*log_r_raw_gq);
            matrix[G,I] r_pred_h = rep_matrix(0,G,I);
            for (idx_non_zero in 1:nnz_y_gq) {
                int g0 = idx_g_gq[idx_non_zero];
                int i0 = idx_i_gq[idx_non_zero];
                r_pred_h[g0,i0] = r_vec[idx_non_zero];
            }

            row_vector[N] wr = to_row_vector((W_geog * (r_pred_h * W_strata_t))');
            row_vector[N] mu_y_h = mu_end[h,] + (wr .* mu_ne_mult);
            
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
