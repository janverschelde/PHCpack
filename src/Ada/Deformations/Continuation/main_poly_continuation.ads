with text_io;                            use text_io;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package Main_Poly_Continuation is

-- DESCRIPTION :
--   This package contains main procedures for two types of homotopies:
--   artificial and natural parameter.

  procedure Ask_Symbol;

  -- DESCRIPTION :
  --   This procedure asks for the symbol to display the additional unknown
  --   when continuation happens in projective space.

  procedure Driver_for_Process_io ( file : in file_type );
  procedure Driver_for_Process_io ( file : in file_type; oc : out natural32 );

  -- DESCRIPTION :
  --   Choice of kind of output information during continuation.

  -- ON ENTRY :
  --   file     must be opened for output.

  -- ON RETURN :
  --   oc       number between 0 and 8 indicating the output code:
  --              0 : no intermediate output information during continuation;
  --              1 : only the final solutions at the end of the paths;
  --              2 : intermediate solutions at each step along the paths;
  --              3 : information of the predictor: t and step length;
  --              4 : information of the corrector: corrections and residuals;
  --              5 : intermediate solutions and information of the predictor;
  --              6 : intermediate solutions and information of the corrector;
  --              7 : information of predictor and corrector;
  --              8 : intermediate solutions, info of predictor and corrector.

  procedure Driver_for_Continuation_Parameters;
  procedure Driver_for_Continuation_Parameters ( file : in file_type );
  procedure Driver_for_Continuation_Parameters ( precision : in natural32 );
  procedure Driver_for_Continuation_Parameters
                ( file : in file_type; precision : in natural32 );

  -- DESCRIPTION :
  --   This procedure allows to determine all continuation parameters.
  --   Writes the final settings on file when there is a file given.
  --   If the precision is given (number of decimal places),
  --   then the default settings for that value of the precision are set
  --   first before allowing the user to interactively change them.

  procedure Check_Continuation_Parameter
                ( sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Check_Continuation_Parameter
                ( sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Check_Continuation_Parameter
                ( sols : in out QuadDobl_Complex_Solutions.Solution_List );
  procedure Check_Continuation_Parameter
                ( sols : in out Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION ;
  --   Reads the value of the continuation parameter for the first solution.
  --   If different from zero, the user is given the opportunity to change it.

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type;
                  p : in Standard_Complex_Poly_Systems.Poly_Sys;
                  prclvl : in natural32;
                  ls : in String_Splitters.Link_to_Array_of_Strings;
                  sols : out Standard_Complex_Solutions.Solution_list;
                  ddsols : out DoblDobl_Complex_Solutions.Solution_list;
                  qdsols : out QuadDobl_Complex_Solutions.Solution_list;
                  mpsols : out Multprec_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 );
  procedure Driver_for_Laurent_Continuation
                ( file : in file_type;
                  p : in Standard_Complex_Laur_Systems.Laur_Sys;
                  prclvl : in natural32;
                  ls : in String_Splitters.Link_to_Array_of_Strings;
                  sols : out Standard_Complex_Solutions.Solution_list;
                  ddsols : out DoblDobl_Complex_Solutions.Solution_list;
                  qdsols : out QuadDobl_Complex_Solutions.Solution_list;
                  mpsols : out Multprec_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This is a driver for the polynomial continuation routine
  --   with an artificial parameter homotopy.
  --   It reads the start system and start solutions and enables the
  --   user to determine all relevant parameters.
  
  -- ON ENTRY :
  --   file       to write diagnostics and results on;
  --   p          a polynomial system;
  --   prclvl     preset precision level, is either 1, 2, or 4,
  --              for double, double double, or quad double precision;
  --   ls         string representations of the input polynomials;
  --   verbose    the verbose level.

  -- ON RETURN :
  --   sols       the computed solutions in standard double precision,
  --              if the precision is set to double;
  --   ddsols     double double precision representation of the solutions,
  --              if the precision is set to double double;
  --   qdsols     quad double precision representation of the solutions,
  --              if the precision is set to quad double;
  --   mpsols     multiprecision representation of the solutions,
  --              if the precision is set to an arbitrary number of decimals.

  procedure Driver_for_Parameter_Continuation
                ( file : in file_type;
                  p : in Standard_Complex_Poly_Systems.Poly_Sys;
                  k : in natural32; target : in Complex_Number;
                  sols : out Standard_Complex_Solutions.Solution_list;
                  verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This is a driver for the polynomial continuation routine
  --   with a natural parameter homotopy.
  --   The start solutions will be read from file.

  -- ON ENTRY :
  --   file       to write diagnostics and results on;
  --   p          a polynomial system, with n equations and n+1 unknowns;
  --   k          index of t = xk;
  --   target     target value for the continuation parameter;
  --   verbose    the verbose level.

  -- ON RETURN :
  --   sols       the computed solutions.

-- DRIVERS FOR DOUBLE DOUBLE, QUAD DOUBLE & MULTIPRECISION CONTINUATION :

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type;
                  p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                  sols : out DoblDobl_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 );
  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type;
                  p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                  sols : out QuadDobl_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 );
  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type; dp : in natural32;
                  p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : out Multprec_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This is a driver for the polynomial continuation routine
  --   with an artificial parameter homotopy.
  --   It reads the start system and start solutions and enables the
  --   user to determine all relevant parameters.
  
  -- ON ENTRY :
  --   file       to write diagnostics and results on;
  --   dp         decimal places in the working precision;
  --   p          a polynomial system;
  --   verbose    the verbose level.

  -- ON RETURN :
  --   sols       the computed solutions;
  --   target     end target value for the continuation parameter.

-- REDEFINING ARTIFIICIAL-PARAMETER HOMOTOPIES

  procedure Standard_Redefine_Homotopy;
  procedure DoblDobl_Redefine_Homotopy;
  procedure QuadDobl_Redefine_Homotopy;

  -- DESCRIPTION :
  --   In case the relaxation power (by default equal to 2) differs
  --   from one, the homotopy will be redefined when polyhedral end
  --   games are called in the path trackers.
  --   The redefinition is provided for artificial-parameter homotopies
  --   in standard double, double double, and quad double precision.

-- CALLING THE PATH TRACKERS :

  procedure Driver_for_Standard_Continuation
                ( file : in file_type;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  proj : in boolean; nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 );
  procedure Driver_for_Standard_Laurent_Continuation
                ( file : in file_type;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  proj : in boolean; nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 );
  procedure Driver_for_Multprec_Continuation
                ( file : in file_type;
                  sols : in out Multprec_Complex_Solutions.Solution_List;
                  proj : in boolean; deci : in natural32;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 );
  procedure Driver_for_DoblDobl_Continuation
                ( file : in file_type;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List;
                  nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 );
  procedure Driver_for_DoblDobl_Laurent_Continuation
                ( file : in file_type;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List;
                  nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 );
  procedure Driver_for_QuadDobl_Continuation
                ( file : in file_type;
                  sols : in out QuadDobl_Complex_Solutions.Solution_List;
                  nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 );
  procedure Driver_for_QuadDobl_Laurent_Continuation
                ( file : in file_type;
                  sols : in out QuadDobl_Complex_Solutions.Solution_List;
                  nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given a homotopy, contained in the Homotopy package,
  --   respectively Standard_Homotopy, Standard_Laurent_Homotopy,
  --   Multprec_Homotopy, DoblDobl_Homotopy, or QuadDobl_Homotopy.
  --   The continuation procedure will be invoked.
  --   The user may tune all continuation parameters.
   
  -- ON ENTRY :
  --   file       to write intermediate results and diagnostics on;
  --   sols       start solutions for the continuation;
  --   deci       number of decimal places;
  --   proj       true when a projective-perpendicular corrector will be used;
  --   nbq        number of equations to turn on the Gauss-Newton correctors;
  --   target     target value for the continuation parameter;
  --   verbose    the verbose level.

  -- ON RETURN :
  --   sols       the computed solutions.

end Main_Poly_Continuation;
