with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Bracket_Polynomials;       use QuadDobl_Bracket_Polynomials;

package QuadDobl_Bracket_Systems is

-- DESCRIPTION :
--   This package contains routines to represent systems of bracket
--   polynomials, with complex coefficients in quad double precision.

  type Bracket_System is array ( integer32 range <> ) of Bracket_Polynomial;

  procedure Clear ( s : in out Bracket_System );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the bracket system.

end QuadDobl_Bracket_Systems;
