with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Random_Product_Start_Systems is

-- DESCRIPTION :
--   This package constructs a random product start
--   product system for a given polynomial system.

  procedure Build_Set_Structure ( i,d : in natural32; p : in Poly );

  -- DESCRIPTION :
  --   Builds a set structure for the polynomial p.

  -- ON ENTRY :
  --   i        index of the polynomial in a system of equations;
  --   d        degree of the polynomial;
  --   p        the polynomial for the i-th equation.

  procedure Build_Set_Structure ( p : in Poly_Sys );

  -- DESCRIPTION :
  --   This is a heuristic procedure for constructing a supporting
  --   set structure of the system p.

  procedure Build_Random_Product_System ( n : in natural32 );

  -- DESCRIPTION :
  --   Based on the set structure, a random linear-product system
  --   will be constructed.  The result is stored in the internal
  --   data manage by the package Random_Product_System.

  -- REQUIRED :
  --   The set structure may not be empty.

  procedure Construct ( p : in Poly_Sys; q : in out Poly_Sys;
		        sols : in out Solution_List );

  -- DESCRIPTION :
  --   Constructs a start system q, with more or less the same 
  --   structure as p.  A heuristic procedure will be used for
  --   constructing a supporting set structure.

end Random_Product_Start_Systems;
