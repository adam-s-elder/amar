cv_test(obs_data = MASS::mvrnorm(n = 100, mu = rep(0,6), Sigma = diag(6)),
        param = "pearson", pos_lp_norms = c(1, 2, 3, 4, "max"),
        ts_ld_bs_samp = 250, num_folds = 1,
        n_bs_smp = 250, nrm_type = "lp",
        f_cv_summary = mean, test_stat_func = "mag",
        incl_chsn_norm = TRUE, test_type = "perm", perf_meas = "mag")
