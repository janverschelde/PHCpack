with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Regular_Solution_Curves_Series is

-- DESCRIPTION :
--   A solution curve is regular if the following three conditions are met:
--   (1) the defining polynomial equations form a complete intersection,
--   in particular: the n polynomials in n+1 variables define curves as
--   the only solutions; there are no higher dimensional solutios;
--   (2) the system is in Noether position: there are no solution curves
--   contained in hyperplanes perpendicular to coordinate axes;
--   (3) the initial form solutions are regular, no multiple solutions.
--   The operations in this package compute series developments for such
--   regular solution curves of Laurent systems.

  procedure Mixed_Cell_Tropisms
              ( report : in boolean;
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision;
                mv : out natural32 );

  -- DESCRIPTION :
  --   Given a system of n Laurent polynomials in n+1 variables,
  --   computes the tropisms, where the last variable is the parameter.
  --   The tropisms are computed via a regular mixed cell configuration,
  --   induced by the lifting defined by the last variable.
  --   As the lifting may differ, even if several supports would be the same,
  --   the type of mixed is assumed to be fully mixed.

  -- ON ENTRY :
  --   report   if intermediate output has to be written to screen;
  --   sup      a list of n supports in n+1 variables.

  -- ON RETURN :
  --   mcc      a mixed cell configuration induced by the lifting
  --            defined by the last exponent in each monomial of p;
  --   mv       the mixed volume of the cells in mcc.

  procedure Mixed_Cell_Tropisms
              ( file : in file_type;
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision;
                mv : out natural32 );

  -- DESCRIPTION :
  --   Given a system of n Laurent polynomials in n+1 variables,
  --   computes the tropisms, where the last variable is the parameter.
  --   The tropisms are computed via a regular mixed cell configuration,
  --   induced by the lifting defined by the last variable.
  --   As the lifting may differ, even if several supports would be the same,
  --   the type of mixed is assumed to be fully mixed.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   sup      a list of n supports in n+1 variables.

  -- ON RETURN :
  --   mcc      a mixed cell configuration induced by the lifting
  --            defined by the last exponent in each monomial of p;
  --   mv       the mixed volume of the cells in mcc.

  procedure Initial_Coefficients
              ( file : in file_type;
                p : in Laur_Sys; mic : in Mixed_Cell;
                psub : out Laur_Sys; sols : out Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver on the subsystem of p
  --   supported by the points that span the mixed cell mic.
  --
  -- ON ENTRY :
  --   file     to write output to file;
  --   p        system of n polynomials in n+1 variables;
  --   mic      a mixed cell with inner normal computed with a lifting
  --            equal to the exponent of the last variable in p.
  --
  -- ON RETURN :
  --   psub     a square subsystem of p with the last variable removed;
  --   sols     solutions to psub, as computed by the blackbox solver.

  procedure Initial_Coefficients
              ( p : in Laur_Sys; mic : in Mixed_Cell;
                psub : out Laur_Sys; sols : out Solution_List;
                report : in boolean );

  -- DESCRIPTION :
  --   Calls the blackbox solver on the subsystem of p
  --   supported by the points that span the mixed cell mic.
  --
  -- ON ENTRY :
  --   p        system of n polynomials in n+1 variables;
  --   mic      a mixed cell with inner normal computed with a lifting
  --            equal to the exponent of the last variable in p;
  --   report   if true, then output is written to screen.
  --
  -- ON RETURN :
  --   psub     a square subsystem of p with the last variable removed;
  --   sols     solutions to psub, as computed by the blackbox solver.

  procedure Shift ( p : in out Poly; verbose : in boolean );
  procedure Shift ( p : in out Laur_Sys; verbose : in boolean );

  -- DESCRIPTION :
  --   Multiplies the monomials in p so that all monomials have 
  --   nonnegative exponents.

  procedure Transform_Coordinates
              ( file : in file_type;
                p : in Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Laur_Sys );

  -- DESCRIPTION :
  --   Applies the unimodular coordinate transformation defined by v
  --   to the system p, the result is in the system q.
  --   Output is written to file.

  procedure Transform_Coordinates
              ( p : in Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Laur_Sys; report : in boolean );

  -- DESCRIPTION :
  --   Applies the unimodular coordinate transformation defined by v
  --   to the system p, the result is in the system q.

  function Initial_Residual
              ( p : in Laur_Sys;
                sol : in Standard_Complex_Vectors.Vector )
              return double_float;

  -- DESCRIPTION :
  --   Returns the residual of the solution extended with zero,
  --   evaluated at p.

  function Initial_Residuals
              ( file : in file_type; p : in Laur_Sys;
                sols : in Solution_List ) return double_float;

  -- DESCRIPTION :
  --   Returns the accumulated sum of all residuals of the solutions,
  --   extended with zero for the last coordinate.
  --   The residuals are written to file for each solution.

  function Initial_Residuals
              ( p : in Laur_Sys;
                sols : in Solution_List;
                report : in boolean ) return double_float;

  -- DESCRIPTION :
  --   Returns the accumulated sum of all residuals of the solutions,
  --   extended with zero for the last coordinate.
  --   If report, then the residuals are written to screen for each solution,
  --   otherwise, this function remains silent.

  procedure Initial
              ( file : in file_type;
                p : in Laur_Sys; mic : in Mixed_Cell;
                tsq : out Poly_Sys; sols : out Solution_List );

  -- DESCRIPTION :
  --   Solves the initial form system defined by the mixed cell mic.

  -- ON ENTRY :
  --   file     to write output;
  --   p        a Laurent system of n equations in n+1 variables,
  --            where the last variable is the lifting;
  --   mic      a mixed cell in the regular mixed cell configuration 
  --            defined by the lower hull of the supports of p.

  -- ON RETURN :
  --   tsq      transformed and shifted polynomial system;
  --   sols     solutions of tsq for t = 0.

  procedure Initial
              ( p : in Laur_Sys; mic : in Mixed_Cell;
                tsq : out Poly_Sys; sols : out Solution_List;
                report : in boolean );

  -- DESCRIPTION :
  --   Solves the initial form system defined by the mixed cell mic.

  -- ON ENTRY :
  --   p        a Laurent system of n equations in n+1 variables,
  --            where the last variable is the lifting;
  --   mic      cell in the regular mixed cell configuration defined
  --            by the lower hull of the supports of p;
  --   report   if true, then output is written to screen.

  -- ON RETURN :
  --   tsq      transformed and shifted polynomial system;
  --   sols     solutions of tsq for t = 0.

  procedure Initials
              ( file : in file_type;
                p : in Laur_Sys; mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the cells in mcc.

  -- ON ENTRY :
  --   file     to write output;
  --   p        a Laurent system of n equations in n+1 variables,
  --            where the last variable is the lifting;
  --   mcc      regular mixed cell configuration defined by the
  --            lower hull of the supports of p.

  procedure Initials
              ( p : in Laur_Sys; mcc : in Mixed_Subdivision;
                report : in boolean );

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the cells in mcc.

  -- ON ENTRY :
  --   p        a Laurent system of n equations in n+1 variables,
  --            where the last variable is the lifting;
  --   mcc      regular mixed cell configuration defined by the
  --            lower hull of the supports of p;
  --   report   if true, then output is written to screen.

end Regular_Solution_Curves_Series;
