with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Tableau_Form_Interface is

-- DESCRIPTION :
--   The functions below export the tableau form of a polynomial system.

  function Tableau_Form_Store
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Stores the polynomial system in double precision
  --   for the given tableau form.

  -- ON ENTRY :
  --   a       in a[0] is the number of equations,
  --           in a[1] is the number of variables,
  --           in a[2] is the total number of monomials,
  --           in a[k+2] is the number of terms of the k-th polynomial;
  --   b       contains the exponents of the monomials;
  --   c       contains the coefficients of the monomials;
  --   vrblvl  is the verbose level.

  function Tableau_Form_Dimensions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimensions of the tableau form of the polynomial
  --   system in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of equations,
  --           in a[1] is the number of variables,
  --           in a[2] is the total number of monomials.

end Tableau_Form_Interface;
