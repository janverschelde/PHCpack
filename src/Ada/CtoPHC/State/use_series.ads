with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_series ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway from C to the power series Newton method.
--   Before calling this procedure, the systems and solutions containers
--   must be initialized properly.  The resulting series are stored in
--   the systems pool, corresponding to each precision level.
--   The verbose level is given by the parameter vrblvl.

-- ON ENTRY :
--   job =  0 : power series Newton method in standard double precision,
--              the input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the maximal degree of the series,
--              a[2] : the number of steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  1 : power series Newton method in double double precision,
--              the input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the maximal degree of the series,
--              a[2] : the number of steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  2 : power series Newton method in quad double precision,
--              the input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the maximal degree of the series,
--              a[2] : the number of steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  3 : same as 0, but instead of solutions, the start terms
--              on the series are given in the systems pool.
--   job =  4 : same as 1, but instead of solutions, the start terms
--              on the series are given in the systems pool.
--   job =  5 : same as 2, but instead of solutions, the start terms
--              on the series are given in the systems pool.
--   job =  6 : power series Newton method in standard double precision,
--              for a Pade approximant.  The input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the degree of the numerator,
--              a[2] : the degree of the denominator,
--              a[3] : the number of the steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  7 : power series Newton method in double double precision,
--              for a Pade approximant.  The input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the degree of the numerator,
--              a[2] : the degree of the denominator,
--              a[3] : the number of the steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  8 : power series Newton method in quad double precision,
--              for a Pade approximant.  The input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the degree of the numerator,
--              a[2] : the degree of the denominator,
--              a[3] : the number of the steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the right range.
