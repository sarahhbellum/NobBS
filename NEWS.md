
# NobBS v. 1.0.0

Updates to v. 0.1.0 (Initial Release) are major fixes to remove
deprecated `tidyverse` functions, replacing `group_by_` with
`group_by(!!sym())` and `vars_select` with standard `select`.

# NobBS v. 0.1.0 (Initial Release)

See `README.md` for a full tutorial on using the NobBS package.

Updates to the development version include:

- A function for stratified Bayesian nowcasting, `NobBS.strat`
- An argument to specify the `proportion_reported` and inflate the
  estimated cases by that factor, important in outbreaks with severe
  under-reporting or asymptomatic cases.
- Additional checks and warnings issued for impossible temporal
  combinations of moving window sizes, maximum delay D, and the length
  of the time series.
- Warning issued when a case onset week is skipped in the time series or
  is missing reports altogether.
- Plot examples

# NobBS v. 0.0.0.9000 (Development)

NobBS (Nowcasting by Bayesian Smoothing) is a new R package that
produces smooth Bayesian nowcasts of incomplete, time-stamped reporting
data implemented from [McGough et
al. 2019](https://www.biorxiv.org/content/10.1101/663823v1.full).
