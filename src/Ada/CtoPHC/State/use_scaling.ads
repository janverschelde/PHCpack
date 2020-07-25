with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_scaling ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Gateway to the coefficient and variable scaling methods in PHCpack.
--   We scale systems and solutions, in double, double double, quad double,
--   and arbitrary multiprecision arithmetic.

-- ON ENTRY :
--   job =     1 : scales the system in the standard container,
--                 on input in a is either 0, 1, or 2:
--                 0 : only equation scaling,
--                 1 : variable scaling without variability reduction
--                 2 : variable scaling with variability reduction;
--                 if a is 1 or 2, then on return c contains the
--                 scaling coefficients and the estimate for the
--                 condition number, stored as a standard complex vector,
--                 the size of this vector is 4*n + 2,
--                 where n equals the number of equations and variables;
--             2 : scales the system in the dobldobl container,
--                 on input in a is either 0, 1, or 2:
--                 0 : only equation scaling,
--                 1 : variable scaling without variability reduction
--                 2 : variable scaling with variability reduction;
--                 if a is 1 or 2, then on return c contains the
--                 scaling coefficients and the estimate for the
--                 condition number, stored as a dobldobl complex vector,
--                 the size of this vector is 8*n + 4,
--                 where n equals the number of equations and variables;
--             3 : scales the system in the quaddobl container,
--                 on input in a is either 0, 1, or 2:
--                 0 : only equation scaling,
--                 1 : variable scaling without variability reduction
--                 2 : variable scaling with variability reduction;
--                 if a is 1 or 2, then on return c contains the
--                 scaling coefficients and the estimate for the
--                 condition number, stored as a quaddobl complex vector,
--                 the size of this vector is 16*n + 8,
--                 where n equals the number of equations and variables;
--             4 : scales the system in the multprec container,
--                 on input in a is either 0, 1, or 2:
--                 0 : only equation scaling,
--                 1 : variable scaling without variability reduction
--                 2 : variable scaling with variability reduction;
--                 if a is 1 or 2, then on return c contains the
--                 scaling coefficients and the estimate for the
--                 condition number, stored as a quadobl complex vector,
--                 the size of this vector is 16*n + 8,
--                 where n equals the number of equations and variables;
--   job =     5 : scales the solutions in the standard container,
--                 in a is the number n, in b the basis, and in c
--                 there are n doubles for use as coefficients,
--                 note that n/2 is the dimension of the complex vector,
--                 and n/2 is then also the dimension of the solutions;
--             6 : scales the solutions in the dobldobl container,
--                 in a is the number n, in b the basis, and in c
--                 there are n doubles for use as coefficients,
--                 note that n/4 is the dimension of the complex vector,
--                 and n/4 is then also the dimension of the solutions;
--             7 : scales the solutions in the quaddobl container,
--                 in a is the number n, in b the basis, and in c
--                 there are n doubles for use as coefficients,
--                 note that n/8 is the dimension of the complex vector,
--                 and n/8 is then also the dimension of the solutions;

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
