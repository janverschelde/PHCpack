with String_Splitters;
with Standard_Complex_Poly_Systems;

package Random_Conditioned_Homotopies is

-- DESCRIPTION :
--   To test variable precision path trackers, we buid homotopies
--   that have at their middle a polynomial system with prescribed
--   condition number at some of its roots.
--   At the start and the end of a random conditioned homotopy are
--   systems with the same support but with random complex coefficients.
--   A good variable precision path tracker should be able to reach all
--   isolated solutions of the target system.

  function Strip_Semicolon ( f : string ) return string;

  -- DESCRIPTION :
  --   Removes the semicolon from the last position in the string f.

  function Random_Coefficient_System
             ( f : String_Splitters.Array_of_Strings )
             return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns a polynomial system with the same supports as f,
  --   but with random complex coefficients on the unit circle.

  function Conditioned_Homotopy
             ( f : String_Splitters.Array_of_Strings )
             return String_Splitters.Array_of_Strings;

  -- DESCRPTION :
  --   On input is a polynomial system f with a prescribed conditioning.
  --   On return is a string representation of a homotopy with 
  --   at t = 0.5, the system f given on input, and 
  --   at t = 0 and 1, random coefficient systems with same support as f.

end Random_Conditioned_Homotopies;
