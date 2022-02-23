test_that("Test returns expected value", {
  dat <- matrix(rnorm(80), ncol = 4)
  tst_res <- mv_pn_test(dat, ic.pearson, test.control(
    pos_lp_norms = c(1, 2, 4, 6, "max"),
    ts_ld_bs_samp = 20
  ))
  expect_type(tst_res, "list")
  expect_equal(length(tst_res$chosen_norm), 5)
  expect_equal(length(tst_res$test_st_eld), 20)
  expect_equal(tst_res$var_mat, NULL)

  tst_res_2 <- mv_pn_test(dat, ic.pearson, test.control(
    pos_lp_norms = c(1, 2, 3, "max"),
    ts_ld_bs_samp = 20,
    other_output = "var_est"
  ))
  expect_type(tst_res_2, "list")
  expect_true("matrix" %in% class(tst_res_2$var_mat))
})
