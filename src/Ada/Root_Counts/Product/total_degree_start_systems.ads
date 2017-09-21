with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Total_Degree_Start_Systems is

-- DESCRIPTION :
--   This package constructs the start system based on the total degree.
--   The i-th equation in a start system based on the total degree is
--     x^(d(i)) - c(i) = 0, for i in p'range and d(i) = degree(p(i)),
--   for some polynomial system p.

  procedure Total_Degree_Info;

  -- DESCRIPTION :
  --   Displays information about the total degree on screen.

  function Degrees ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Natural_Vectors.Vector;
  function Degrees ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Natural_Vectors.Vector;
  function Degrees ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector with degrees of the polynomials in p.

  function Product ( d : Standard_Natural_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the product of the elements in d.
  --   Product(Degrees(p)) computes the total degree of the system p.

  function Total_Degree
             ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return natural32;
  function Total_Degree
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return natural32;
  function Total_Degree
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return natural32;

  -- DESCRIPTION :
  --   Returns the total degree of the polynomial system p.

-- CREATE THE SYSTEM :

  function Start_System
             ( d : Standard_Natural_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Start_System
             ( d : Standard_Natural_Vectors.Vector;
               c : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns a start system for the degrees in d and righthand side c.
  --   If c is omitted, then random coefficients are used instead.

  -- REQUIRED : d'range = c'range.

  function Start_System
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Start_System
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               c : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns a start system based on the total degree of p.

  -- ON ENTRY :
  --   p       a polynomial system;
  --   c       a vector of constants, to be used as righthand side,
  --           if omitted, random complex numbers will be generated.

  function Coefficients
             ( q : in Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Vectors.Vector;
  function Coefficients
             ( q : in DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Vectors.Vector;
  function Coefficients
             ( q : in QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of complex coefficients used to create
  --   the total degree start system q,
  --   in particular we have: Coefficients(Start_System(p,c)) = c.

-- PARTICULAR SOLVERS :

  function Eval ( d : in Standard_Natural_Vectors.Vector;
                  c,x : in Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of x^d - c, for use as residual
  --   when x is a root.

  function Root ( d : in Standard_Natural_Vectors.Vector;
                  s : in Standard_Natural_Vectors.Vector;
                  c : in Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector;
  function Root ( d : in Standard_Natural_Vectors.Vector;
                  s : in Standard_Natural_Vectors.Vector;
                  c : in DoblDobl_Complex_Vectors.Vector )
                return DoblDobl_Complex_Vectors.Vector;
  function Root ( d : in Standard_Natural_Vectors.Vector;
                  s : in Standard_Natural_Vectors.Vector;
                  c : in QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a selected root of c, the i-th entry of the vector on 
  --   return is the s(i)-th root of c(i)^(1/d(i)), for i in c'range.

  function Create ( r : Standard_Complex_Vectors.Vector )
                  return Standard_Complex_Solutions.Solution;
  function Create ( r : Standard_Complex_Vectors.Vector;
                    rcond : double_float )
                  return Standard_Complex_Solutions.Solution;
  function Create ( r : DoblDobl_Complex_Vectors.Vector )
                  return DoblDobl_Complex_Solutions.Solution;
  function Create ( r : DoblDobl_Complex_Vectors.Vector;
                    rcond : double_double )
                  return DoblDobl_Complex_Solutions.Solution;
  function Create ( r : QuadDobl_Complex_Vectors.Vector )
                  return QuadDobl_Complex_Solutions.Solution;
  function Create ( r : QuadDobl_Complex_Vectors.Vector;
                    rcond : quad_double )
                  return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns the default encapsulation of a complex vector as a solution,
  --   with optional argument the inverse for the condition number.

  function Solve ( d : in Standard_Natural_Vectors.Vector;
                   c : in Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns all solutions of the system x^d - c = 0.

-- COMPREHENSIVE SOLVERS :

  procedure Start_System 
             ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
               q : out Standard_Complex_Poly_Systems.Poly_Sys;
               c : in Standard_Complex_Vectors.Vector;
               qsols : out Standard_Complex_Solutions.Solution_List );

  -- ON ENTRY :
  --   p       the polynomial system that has to be solved;
  --   c       a vector with constants.

  -- ON RETURN :
  --   q       the start system, with for i in q'range :
  --           q(i) = x_i^Degree(p(i)) - c(i);
  --   qsols   the solutions of the start system.

  procedure Start_System 
             ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
               q : out Standard_Complex_Poly_Systems.Poly_Sys;
               qsols : out Standard_Complex_Solutions.Solution_List ); 

  -- DESCRIPTION :
  --   The meaning of the parameters is here the same, except that instead
  --   of the vector c, a random number generator will be used for choosing
  --   the c(i)'s on the unit circle.

-- RESIDUAL CALCULATION :

  function Sum_Residuals
             ( q : Standard_Complex_Poly_Systems.Poly_Sys;
               qsols : Standard_Complex_Solutions.Solution_List )
             return double_float;

  -- DESCRIPTION :
  --   Returns the sum of all residuals of the start solutions qsols,
  --   of the start system q.  The residual is the max norm of the solution
  --   vector evaluated at the system q.

  procedure Write_Residuals
              ( file : in file_type;
                q : in Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : in Standard_Complex_Solutions.Solution_List;
                sum : out double_float );

  -- DESCRIPTION :
  --   Writes the residual for each solution vector to screen.
  --   The residual is the max norm of the vector evaluation at q.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   q        a start system;
  --   qsols    solutions to th start system.

  -- ON RETURN :
  --   sum      sum of all residuals.

end Total_Degree_Start_Systems;
