context("contribution arguments")

test_that("check_table with ok tables", {

  expect_true(check_table(matrix(rnorm(200), ncol = 10, byrow = TRUE)))

})

test_that("check_table fails with invalid classes", {

  expect_error(check_table(c('one', 'two', 'three')))
  expect_error(check_table(TRUE))

})

test_that("check_rows with ok rows", {

  expect_true(check_rows(matrix(rnorm(200), nrow = 10),
                         matrix(rnorm(100), nrow = 10)))

})

test_that("check_rows fails with invalid row numbers", {

  expect_error(check_rows(matrix(rnorm(200), nrow = 20),
                         matrix(rnorm(100), nrow = 10)))

})
