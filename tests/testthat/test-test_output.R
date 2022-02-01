test_that("Test returns expected value", {
  dat <- matrix(rnorm(80), ncol = 4)
  tst_res <- mv_pn_test(dat, ic.pearson, test.control(
    pos_lp_norms = c(1, 2, 4, 6, "max"),
    ts_ld_bs_samp = 20,
    more_info = TRUE
  ))
  expect_type(tst_res, "list")
  expect_equal(length(tst_res$chosen_norm), 5)
  expect_equal(length(tst_res$test_st_eld), 20)
  expect_equal(tst_res$var_mat, NULL)
  tst_res_2 <- mv_pn_test(dat, ic.pearson, test.control(
    pos_lp_norms = c(1, 2, 3, "max"),
    ts_ld_bs_samp = 20,
    more_info = FALSE
  ))
  expect_type(tst_res_2, "double")

  tst_res_3 <- mv_pn_test(dat, ic.pearson, test.control(
    pos_lp_norms = c(1, 2, 3, "max"),
    ts_ld_bs_samp = 20,
    more_info = TRUE,
    ret_cov_mat = TRUE
  ))
  expect_type(tst_res_3, "list")
  expect_true("matrix" %in% class(tst_res_3$var_mat))
})
