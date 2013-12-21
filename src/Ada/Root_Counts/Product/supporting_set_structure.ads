with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Supporting_Set_Structure is

-- DESCRIPTION :
--   A set structure models the product structure of a polynomial system.
--   The functions offered by this package allow to verify if the set
--   structure supports a given polynomial system.

-- REQUIRED : Set_Structure.Empty is false.

  function Is_Supporting
             ( p : Poly; i : natural32; verbose : boolean ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the i-th set structure supports the polynomial p.
  --   For this, every variable must occur in at least as many distinct
  --   sets as the degree of p in that variable.
  --   If verbose, then the degree of the polynomial in each variable
  --   and the count of each variable in the set structure is written.

  function Is_Supporting ( p : Poly_Sys; verbose : boolean ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the set structure supports the polynomials in p.
  --   If verbose, then degrees and counts are written to screen
  --   which allows to track the cause for a false return.

end Supporting_Set_Structure;
