# NEWS

## TRMF 0.1.5

### Enhancements

  - Added retrain() function to allow for warm starts and restarts.
  - Setting numit = 0 in train() now initializes the computation, but doesn't do any iterations


### Bug fixes

  - Fixed a bug where weights for regression with factorization were incorrectly included in the calculation.
  


## TRMF 0.1.2

### Enhancements

  - Updated some of the internal linear algebra for 2-5x speedup
  
  - TRMF_summary now prints a weighted R2 for weighted data.


### Bug fixes

  - Fixed a bug in the component function which only returned the Xm component.
  
  - Fixed a bug in NormalizeMatrix which returns NA's in some cases.
  
### Changes

  - Renamed TRMF_coefficients() to TRMF_columns() to make less confusing. TRMF_coefficients() will be deprecated future releases.
