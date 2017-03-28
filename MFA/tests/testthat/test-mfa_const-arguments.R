context("mfa_const arguments")

test_that("check_dataset with ok datasets", {

  expect_true(check_dataset(matrix(rnorm(200), ncol = 10, byrow = TRUE)))
  expect_true(check_dataset(as.data.frame(matrix(rnorm(200), ncol = 10, byrow = TRUE))))

})

test_that("check_dataset fails with invalid classes", {

  expect_error(check_dataset(c('one', 'two', 'three')))
  expect_error(check_dataset(TRUE))

})

test_that("check_sets with ok sets", {

  expect_true(check_sets(list(1:3, 4:5, 6:10)))
  expect_true(check_sets(list(c('one', 'two', 'three'))))

})

test_that("check_sets fails with invalid classes", {

  expect_error(check_sets(list(1:3, 4:7, TRUE)))

})

test_that("check_ncomps with okay values", {

  expect_true(check_ncomps(5))
  expect_true(check_ncomps(5.0))

})

test_that("check_ncomps fails with certain values", {

  expect_error(check_ncomps(c("one")))
  expect_error(check_ncomps(5.4))

})

test_that("check_center with okay values", {

  expect_true(check_center(TRUE))
  expect_true(check_center(c(5, 4, 6)))

})

test_that("check_center fails with certain values", {

  expect_error(check_center(c("one")))

})

test_that("check_scale with okay values", {

  expect_true(check_scale(TRUE))
  expect_true(check_scale(c(5, 4, 6)))

})

test_that("check_scale fails with certain values", {

  expect_error(check_scale(c("one")))

})
