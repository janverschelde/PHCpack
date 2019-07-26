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

  function Standard_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function DoblDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function QuadDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Multprec_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns a randomly generated polynomial system, with double,
  --   double double, quad double, or arbitrary precision coefficients.

  -- ON ENTRY :
  --   nvr     number of variables;
  --   deg     largest degree of the monomials;
  --   mct     monomial count, if zero, then the polynomials are dense;
  --   ctp     coefficient type is one of the following:
  --             0 : on complex unit circle (recommended),
  --             1 : the constant one (useful for templates),
  --             2 : a random float between -1.0 and +1.0;
  --   verbose is the verbose flag, if true, then the generated polynomials
  --           are written to screen, if false, no output is written.

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
