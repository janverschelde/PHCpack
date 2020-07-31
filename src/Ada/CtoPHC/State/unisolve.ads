with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function unisolve ( job : integer32;
                    a : C_intarrs.Pointer;
                    b : C_intarrs.Pointer;
                    c : C_dblarrs.Pointer;
                    vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Gateway to the univariate solvers in PHCpack.

-- ON ENTRY :
--   job   =   1 : standard precision univariate solver of polynomial
--                 which must be in the systems container, on input are
--                 in a: the maximum number of iterations;
--                 in c: the accuracy requirement, and on return are
--                 in b: the number of iterations spent and the solution
--                 container contains the solutions.
--   job   =   2 : double double precision univariate solver of polynomial
--                 which must be in the systems container, on input are
--                 in a: the maximum number of iterations;
--                 in c: the accuracy requirement, and on return are
--                 in b: the number of iterations spent and the solution
--                 container contains the solutions.
--   job   =   3 : quad double precision univariate solver of polynomial
--                 which must be in the systems container, on input are
--                 in a: the maximum number of iterations;
--                 in c: the accuracy requirement, and on return are
--                 in b: the number of iterations spent and the solution
--                 container contains the solutions.
--   job   =   4 : multiprecision univariate solver of polynomial
--                 which must be in the systems container, on input are
--                 in a[0]: number of decimal places in the working precision,
--                    a[1]: the maximum number of iterations,
--                 in c: the accuracy requirement, and on return are
--                 in b: the number of iterations spent and the solution
--                 container contains the solutions.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
