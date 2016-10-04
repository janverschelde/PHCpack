with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;

package Random_Polynomial_Systems is

-- DESCRIPTION :
--   This package offers interactive generators of systems of polynomials,
--   generated at random, for coefficients in double, double double,
--   quad double, and arbitrary multiprecision precision.

  procedure Save_System_to_File
              ( s : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Save_System_to_File
              ( s : in DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Save_System_to_File
              ( s : in QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Save_System_to_File
              ( s : in Multprec_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Prompts the user for a file name
  --   and writes the system to file.

  procedure Standard_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure DoblDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure QuadDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Multprec_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with complex coefficients, in double, double double,
  --   quad double, or quad double precision.
  --   The procedure is interactive, the result is returned in lp.

end Random_Polynomial_Systems;
