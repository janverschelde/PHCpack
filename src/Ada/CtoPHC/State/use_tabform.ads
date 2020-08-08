with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_tabform ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to the systems container
--   via the tableau format of polynomial systems.

-- ON ENTRY :
--   job    =   0 : initializes the systems container with the system 
--                  defined by tableau form stored in the triplet (a, b, c):
--                  in a[0] is the number of equations, 
--                  in a[1] is the number of variables,
--                  in a[2] is the total number of monomials,
--                  in a[k+2] is the number of terms of the k-th polynomial,
--                  exponents of monomials are stored in b,
--                  the size of b is a[1]*a[2],
--                  complex coefficients of monomials are stored in c,
--                  the size of c is 2*a[2];
--          =   1 : returns in a the dimensions of the system in the
--                  container for the standard double precision systems,
--                  a[0] : the number of equations,
--                  a[1] : the number of variables, and
--                  a[2] : the total number of terms.
--                  This job is good to allocate the data to retrieve
--                  the tableau form of a system.

--
-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: indices to monomial out of range, or job not in the proper range.
