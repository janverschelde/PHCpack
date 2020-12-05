with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Poly_Systems;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_VecVecs;
with TripDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with TripDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;

package Run_Power_Series_Methods is

-- DESCRIPTION :
--   Newton's method on power series can start at a constant or at
--   power series, in double, double double, or quad double precision.

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in standard double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions;
  --   vrb     the verbose level.

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in double double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton or not;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions;
  --   vrb     the verbose level.

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in TripDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in triple double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton or not;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions;
  --   vrb     the verbose level.

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in quad double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions.

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in standard double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions;
  --   vrb     the verbose level.

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in double double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions;
  --   vrb     the verbose level.

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in TripDobl_Complex_Poly_Systems.Poly_Sys;
               s : in TripDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in triple double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions;
  --   vrb     the verbose level.

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in quad double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions;
  --   vrb     the verbose level.

  procedure Standard_Main_at_Constant
               ( infilename,outfilename : in string;
                 vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the infilename is not the empty string, then the polynomial
  --   system and solutions will be read from the file called infilename,
  --   else prompts for a system and solutions.
  --   Power series starting starting at a constant leading term 
  --   are computed in double precision.
  --   If the outfilename is not the empty string, then the output
  --   will be written to the file called outfilename,
  --   else prompts for a file name to write the output to.
  --   The vrb on input is the verbose level.

  procedure DoblDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the infilename is not the empty string, then the polynomial
  --   system and solutions will be read from the file called infilename,
  --   else prompts for a system and solutions.
  --   Power series starting starting at a constant leading term 
  --   are computed in double double precision.
  --   If the outfilename is not the empty string, then the output
  --   will be written to the file called outfilename,
  --   else prompts for a file name to write the output to.
  --   The vrb on input is the verbose level.

  procedure TripDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the infilename is not the empty string, then the polynomial
  --   system and solutions will be read from the file called infilename,
  --   else prompts for a system and solutions.
  --   Power series starting starting at a constant leading term 
  --   are computed in triple double precision.
  --   If the outfilename is not the empty string, then the output
  --   will be written to the file called outfilename,
  --   else prompts for a file name to write the output to.
  --   The vrb on input is the verbose level.

  procedure QuadDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the infilename is not the empty string, then the polynomial
  --   system and solutions will be read from the file called infilename,
  --   else prompts for a system and solutions.
  --   Power series starting starting at a constant leading term 
  --   are computed in quad double precision.
  --   If the outfilename is not the empty string, then the output
  --   will be written to the file called outfilename,
  --   else prompts for a file name to write the output to.
  --   The vrb on input is the verbose level.

  procedure Standard_Main_at_Series
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the infilename is not the empty string, then the polynomial
  --   system will be read from the file called infilename,
  --   else prompts for a system.
  --   Prompts for a series to run Newton's method in double precision.
  --   If the outfilename is not the empty string, then the output
  --   will be written to the file called outfilename,
  --   else prompts for a file name to write the output to.
  --   The vrb on input is the verbose level.

  procedure DoblDobl_Main_at_Series
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the infilename is not the empty string, then the polynomial
  --   system will be read from the file called infilename,
  --   else prompts for a system.  Prompts for a series to run
  --   Newton's method in double double precision.
  --   If the outfilename is not the empty string, then the output
  --   will be written to the file called outfilename,
  --   else prompts for a file name to write the output to.
  --   The vrb on input is the verbose level.

  procedure QuadDobl_Main_at_Series
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the infilename is not the empty string, then the polynomial
  --   system will be read from the file called infilename,
  --   else prompts for a system.  Prompts for a series to run
  --   Newton's method in quad double precision.
  --   If the outfilename is not the empty string, then the output
  --   will be written to the file called outfilename,
  --   else prompts for a file name to write the output to.
  --   The vrb on input is the verbose level.

end Run_Power_Series_Methods;
