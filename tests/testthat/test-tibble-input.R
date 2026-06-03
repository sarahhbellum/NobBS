test_that("NobBS accepts tibble input", {
  skip_if_not_installed("tibble")

  data(denguedat)

  expect_no_error(
    NobBS(
      data = tibble::as_tibble(denguedat),
      now = as.Date("1990-04-09"),
      units = "1 week",
      onset_date = "onset_week",
      report_date = "report_week",
      specs = list(
        nSamp = 10,
        nBurnin = 10,
        nAdapt = 10,
        nChains = 1
      ),
      quiet = TRUE
    )
  )
})


test_that("NobBS.strat accepts tibble input", {
  skip_if_not_installed("tibble")

  data(denguedat)

  expect_no_error(
    NobBS.strat(
      data = tibble::as_tibble(denguedat),
      now = as.Date("1990-02-05"),
      units = "1 week",
      onset_date = "onset_week",
      report_date = "report_week",
      strata = "gender",
      specs = list(
        nSamp = 10,
        nBurnin = 10,
        nAdapt = 10,
        nChains = 1
      ),
      quiet = TRUE
    )
  )
})
