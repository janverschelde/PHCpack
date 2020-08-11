with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_solpool ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to the pool of solution lists.

-- ON ENTRY :
--   job      = 0 : initialize solutions pool with n = a[0];
--            = 1 : get the size of the solutions pool returned in a[0];
--            = 2 : returns the length of the k-th solution list,
--                  where k = a[0] on entry and length = b[0] on return;
--            = 3 : returns the dimension of the k-th solution list,
--                  where k = a[0] on entry and dimension = b[0] on return;
--            = 4 : appends a solution to the k-th solution list,
--                  where k = a[0], n = b[0], m = b[1], and c is an array
--                  of 2*n+5 doubles in the following order:
--                  two doubles for the complex continuation parameter t,
--                  2*n doubles for the coefficients of the solution vector,
--                  one double for the norm of last Newton update,
--                  one double for the inverse of condition# estimate,
--                  one double for the norm of the residual;
--            = 5 : retrieves the i-th solution of the k-th list,
--                  where k = a[0] and i = a[1] on entry,
--                  on return in b[0] is m and c is the array of doubles
--                  of dimension 2*n+5 to store the solution data.
-- 
--   a        memory allocated a natural number, either the size
--            of the systems pool or the index of a system in the pool;
--   b        discrete information about a solution;
--   c        continuous data about a solution.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
