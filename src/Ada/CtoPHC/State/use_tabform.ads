with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_tabform ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to the systems container
--   via the tableau format of polynomial systems.

-- ON ENTRY :
--   job    =   0 : the triplet (a, b, c) defines a tableau form,
--                  in a[0] is the number of equations, 
--                  in a[1] is the number of variables,
--                  in a[2] is the total number of monomials,
--                  in a[k+2] is the number of terms of the k-th polynomial,
--                  exponents of monomials are stored in b,
--                  the size of b is a[1]*a[2],
--                  complex coefficients of monomials are stored in c,
--                  the size of c is 2*a[2].

--
-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: indices to monomial out of range, or job not in the proper range.
