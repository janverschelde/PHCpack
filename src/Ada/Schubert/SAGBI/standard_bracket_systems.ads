with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;

package Standard_Bracket_Systems is

-- DESCRIPTION :
--   This package contains routines to manipulate systems of
--   bracket polynomials.

  type Bracket_System is array ( integer32 range <> ) of Bracket_Polynomial;

  function Straightening_Syzygies ( n,d : natural32 ) return Bracket_System;

  -- DESCRIPTION :
  --   Returns the system of straightening syzygies that forms a Groebner
  --   basis for the ideal of all d-planes in n-space.

  procedure Clear ( s : in out Bracket_System );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the bracket system.

end Standard_Bracket_Systems;
