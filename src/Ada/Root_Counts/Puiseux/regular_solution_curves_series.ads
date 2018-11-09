with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with QuadDobl_CSeries_Poly_Systems;

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

  function Tropisms ( mcc : Mixed_Subdivision; mixvol : natural32 )
                    return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the normals to the cells in the cell configuration mcc.
  --   The sum of the mixed volumes of all cells in mcc equals mixvol
  --   and 1..mixvol is the range of the vector on return.
  --   If the mixed cell has mixed volume mv, then mv copies of its
  --   inner normal are in the vector on return.

  procedure Initial_Coefficients
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Initial_Coefficients
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Initial_Coefficients
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the blackbox solver on the subsystem of p
  --   supported by the points that span the mixed cell mic,
  --   in double, double double, or quad double precision.
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
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                report : in boolean );
  procedure Initial_Coefficients
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                report : in boolean );
  procedure Initial_Coefficients
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                report : in boolean );

  -- DESCRIPTION :
  --   Calls the blackbox solver on the subsystem of p
  --   supported by the points that span the mixed cell mic,
  --   in double, double double, and quad double precision.
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

  procedure Shift ( file : in file_type;
                    p : in out Standard_Complex_Laurentials.Poly );
  procedure Shift ( file : in file_type;
                    p : in out DoblDobl_Complex_Laurentials.Poly );
  procedure Shift ( file : in file_type;
                    p : in out QuadDobl_Complex_Laurentials.Poly );
  procedure Shift ( file : in file_type;
                    p : in out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Shift ( file : in file_type;
                    p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Shift ( file : in file_type;
                    p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Multiplies the monomials in p so that all monomials have 
  --   nonnegative exponents.  Output is written to file.

  procedure Shift ( p : in out Standard_Complex_Laurentials.Poly;
                    verbose : in boolean );
  procedure Shift ( p : in out DoblDobl_Complex_Laurentials.Poly;
                    verbose : in boolean );
  procedure Shift ( p : in out QuadDobl_Complex_Laurentials.Poly;
                    verbose : in boolean );
  procedure Shift ( p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                    verbose : in boolean );
  procedure Shift ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    verbose : in boolean );
  procedure Shift ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    verbose : in boolean );

  -- DESCRIPTION :
  --   Multiplies the monomials in p so that all monomials have 
  --   nonnegative exponents.  If verbose, the output is written to screen.

  function Transform ( d,v : Standard_Integer_Vectors.Vector )
                     return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Transforms the degrees in d, using the inner normal in v.

  function Transform ( t : Standard_Complex_Laurentials.Term;
                       v : Standard_Integer_Vectors.Vector )
                     return Standard_Complex_Laurentials.Term;
  function Transform ( t : DoblDobl_Complex_Laurentials.Term;
                       v : Standard_Integer_Vectors.Vector )
                     return DoblDobl_Complex_Laurentials.Term;
  function Transform ( t : QuadDobl_Complex_Laurentials.Term;
                       v : Standard_Integer_Vectors.Vector )
                     return QuadDobl_Complex_Laurentials.Term;

  -- DESCRIPTION :
  --   Transforms the monomial defined by t, using the inner normal in v,
  --   for coefficients in double, double double, or quad double precision.

  function Transform ( p : Standard_Complex_Laurentials.Poly;
                       v : Standard_Integer_Vectors.Vector )
                     return Standard_Complex_Laurentials.Poly;
  function Transform ( p : DoblDobl_Complex_Laurentials.Poly;
                       v : Standard_Integer_Vectors.Vector )
                     return DoblDobl_Complex_Laurentials.Poly;
  function Transform ( p : QuadDobl_Complex_Laurentials.Poly;
                       v : Standard_Integer_Vectors.Vector )
                     return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Transforms the polynomial p, using the inner normal in v,
  --   for coefficients in double, double double, or quad double precision.

  function Transform ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                       v : Standard_Integer_Vectors.Vector )
                     return Standard_Complex_Laur_Systems.Laur_Sys;
  function Transform ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                       v : Standard_Integer_Vectors.Vector )
                     return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Transform ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                       v : Standard_Integer_Vectors.Vector )
                     return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Transforms the system p, using the inner normal in v,
  --   for coefficients in double, double double, or quad double precision.

  procedure Transform_Coordinates
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Transform_Coordinates
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Transform_Coordinates
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Applies the unimodular coordinate transformation defined by v
  --   to the system p, the result is in the system q.
  --   Output is written to file.

  procedure Transform_Coordinates
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                report : in boolean );
  procedure Transform_Coordinates
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean );
  procedure Transform_Coordinates
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean );

  -- DESCRIPTION :
  --   Applies the unimodular coordinate transformation defined by v
  --   to the system p, the result is in the system q.
  --   If report, then output is written to screen.

  function Initial_Residual
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sol : in Standard_Complex_Vectors.Vector )
              return double_float;
  function Initial_Residual
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector )
              return double_float;
  function Initial_Residual
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sol : in QuadDobl_Complex_Vectors.Vector )
              return double_float;

  -- DESCRIPTION :
  --   Returns the residual of the solution extended with zero,
  --   evaluated at p, in double, double double, or quad double precision.

  function Initial_Residuals
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List )
              return double_float;
  function Initial_Residuals
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List )
              return double_float;
  function Initial_Residuals
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List )
              return double_float;

  -- DESCRIPTION :
  --   Returns the accumulated sum of all residuals of the solutions,
  --   extended with zero for the last coordinate.
  --   The residuals are written to file for each solution.

  function Initial_Residuals
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                report : in boolean ) return double_float;
  function Initial_Residuals
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                report : in boolean ) return double_float;
  function Initial_Residuals
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                report : in boolean ) return double_float;

  -- DESCRIPTION :
  --   Returns the accumulated sum of all residuals of the solutions,
  --   extended with zero for the last coordinate.
  --   If report, then the residuals are written to screen for each solution,
  --   otherwise, this function remains silent.

  procedure Initial
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Initial
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Initial
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

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
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                report : in boolean );
  procedure Initial
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                report : in boolean );
  procedure Initial
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
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

  function Series ( file : file_type;
                    p : Standard_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : Standard_Complex_Vectors.Vector;
                    maxdeg,nit : integer32 )
                  return Standard_Complex_Series_Vectors.Vector;
  function Series ( file : file_type;
                    p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : DoblDobl_Complex_Vectors.Vector;
                    maxdeg,nit : integer32 )
                  return DoblDobl_Complex_Series_Vectors.Vector;
  function Series ( file : file_type;
                    p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : QuadDobl_Complex_Vectors.Vector;
                    maxdeg,nit : integer32 )
                  return QuadDobl_Complex_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies Newton's method to p, starting at xt0,
  --   using as many iterations as the value of nit,
  --   writing output to file.
  --   The vector on returns contains the truncated series
  --   expansion for the solution with leading coefficients xt0.

  function Series ( p : Standard_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : Standard_Complex_Vectors.Vector;
                    maxdeg,nit : integer32; report : boolean )
                  return Standard_Complex_Series_Vectors.Vector;
  function Series ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : DoblDobl_Complex_Vectors.Vector;
                    maxdeg,nit : integer32; report : boolean )
                  return DoblDobl_Complex_Series_Vectors.Vector;
  function Series ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : QuadDobl_Complex_Vectors.Vector;
                    maxdeg,nit : integer32; report : boolean )
                  return QuadDobl_Complex_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies Newton's method to p, starting at xt0,
  --   using as many iterations as the value of nit,
  --   writing output to screen if report is true.
  --   The vector on returns contains the truncated series
  --   expansion for the solution with leading coefficients xt0.

  function Series ( file : file_type;
                    p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32 )
                  return Standard_Complex_Series_VecVecs.VecVec;
  function Series ( file : file_type;
                    p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32 )
                  return DoblDobl_Complex_Series_VecVecs.VecVec;
  function Series ( file : file_type;
                    p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32 )
                  return QuadDobl_Complex_Series_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Applies Newton's method to compute series solutions,
  --   starting at the initial solutions in the list sols.

  -- ON ENTRY :
  --   file     to write output to;
  --   p        system with solutions for t = 0;
  --   sols     solutions of p for t = 0;
  --   maxdeg   the maximal degree of the series;
  --   nit      number of iterations of Newton's method.

  -- ON RETURN :
  --   The sequence of series with leading coefficients in sols.

  function Series ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32; report : boolean )
                  return Standard_Complex_Series_VecVecs.VecVec;
  function Series ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32; report : boolean )
                  return DoblDobl_Complex_Series_VecVecs.VecVec;
  function Series ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32; report : boolean )
                  return QuadDobl_Complex_Series_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Applies Newton's method to compute series solutions,
  --   starting at the initial solutions in the list sols.

  -- ON ENTRY :
  --   p        system with solutions for t = 0;
  --   sols     solutions of p for t = 0;
  --   maxdeg   the maximal degree of the series;
  --   nit      number of iterations of Newton's method;
  --   report   flag to write output to screen.

  -- ON RETURN :
  --   The sequence of series with leading coefficients in sols.

  procedure Concat ( sto : in out Standard_Complex_Series_VecVecs.VecVec;
                     idx : in out integer32;
                     sfrom : in Standard_Complex_Series_VecVecs.VecVec );
  procedure Concat ( sto : in out DoblDobl_Complex_Series_VecVecs.VecVec;
                     idx : in out integer32;
                     sfrom : in DoblDobl_Complex_Series_VecVecs.VecVec );
  procedure Concat ( sto : in out QuadDobl_Complex_Series_VecVecs.VecVec;
                     idx : in out integer32;
                     sfrom : in QuadDobl_Complex_Series_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Concatenates the series in sfrom to the series in sto,
  --   incrementing the index idx with each added series.

  -- REQUIRED : idx + sfrom'last <= sto'last.

  -- ON ENTRY :
  --   sto      vector of series vectors up to idx;
  --   idx      index of the last valid entry in sto;
  --   sfrom    series to be concatenated to sto.

  -- ON RETURN :
  --   sto      updated vector of series vectors up to idx;
  --   idx      updated index of the last valid entry in sto;

  function Series ( file : file_type;
                    p : Standard_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32 )
                  return Standard_Complex_Series_VecVecs.VecVec;
  function Series ( file : file_type;
                    p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32 )
                  return DoblDobl_Complex_Series_VecVecs.VecVec;
  function Series ( file : file_type;
                    p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32 )
                  return QuadDobl_Complex_Series_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the cells in mcc
  --   and computes series developments for the solution curves,
  --   writing output to file.

  -- ON ENTRY :
  --   file     to write output;
  --   p        a Laurent system of n equations in n+1 variables,
  --            where the last variable is the lifting;
  --   mcc      regular mixed cell configuration defined by the
  --            lower hull of the supports of p;
  --   mv       the mixed volume of the supports of p;
  --   maxdeg   the maximal degree of the series;
  --   nit      number of Newton steps per solution series.

  -- ON RETURN :
  --   A sequence of mv series in the range 1..mv.

  function Series ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32; report : boolean )
                  return Standard_Complex_Series_VecVecs.VecVec;
  function Series ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32; report : boolean )
                  return DoblDobl_Complex_Series_VecVecs.VecVec;
  function Series ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32; report : boolean )
                  return QuadDobl_Complex_Series_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the cells in mcc
  --   and computes series developments for the solution curves,
  --   writing output to screen if the report flag is true.

  -- ON ENTRY :
  --   p        a Laurent system of n equations in n+1 variables,
  --            where the last variable is the lifting;
  --   mcc      regular mixed cell configuration defined by the
  --            lower hull of the supports of p;
  --   mv       the mixed volume of the supports of p;
  --   maxdeg   the maximal degree of the series;
  --   nit      number of Newton steps per solution series;
  --   report   if true, then output is written to screen.

  -- ON RETURN :
  --   A sequence of mv series in the range 1..mv.

end Regular_Solution_Curves_Series;
