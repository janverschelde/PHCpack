with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_series ( job : integer32;
                      a : C_intarrs.Pointer;
		      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer ) return integer32;

-- DESCRIPTION :
--   Provides a gateway from C to the power series Newton method.
--   Before calling this procedure, the systems and solutions containers
--   must be initialized properly.  The resulting series are stored in
--   the systems pool, corresponding to each precision level.

-- ON ENTRY :
--   job =  0 : power series Newton method in standard double precision,
--              the input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the number of steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  1 : power series Newton method in double double precision,
--              the input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the number of steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  2 : power series Newton method in quad double precision,
--              the input parameters are
--              a[0] : the index of the series parameter,
--              a[1] : the number of steps in Newton's method,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the right range.
