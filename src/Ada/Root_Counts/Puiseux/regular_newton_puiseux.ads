with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;

package Regular_Newton_Puiseux is

-- DESCRIPTION :
--   Development of the Newton-Puiseux algorithm,
--   for regular solution curves, defined by complete intersections,
--   in Noether position, with sufficiently general coefficients.

  procedure Tropisms_by_Mixed_Cells
              ( file : in file_type; 
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32 );

  -- DESCRIPTION :
  --   Computes a regular mixed cell configuration for
  --   the supports in sup, with some writing to file.
  --   The mixed volume is in mv.

  procedure Tropisms_by_Mixed_Cells
              ( sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                report : in boolean );

  -- DESCRIPTION :
  --   Computes a regular mixed cell configuration for
  --   the supports in sup, with some writing to screen if report.
  --   The mixed volume is in mv.

  function Standard_Residual
              ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                s : Standard_Complex_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : double_float ) return double_float;
  function DoblDobl_Residual
              ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : DoblDobl_Complex_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : double_double ) return double_double;
  function QuadDobl_Residual
              ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : QuadDobl_Complex_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Returns the residual of the evaluation of s at t in p,
  --   or ||p(s(t),t)||, with weights for the exponents in w,
  --   in standard double, double double, or quad double precision.

  function Standard_Residuals
              ( file : file_type;
                p : Standard_Complex_Laur_Systems.Laur_Sys;
                s : Standard_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_float ) return double_float;
  function DoblDobl_Residuals
              ( file : file_type;
                p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : DoblDobl_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_double ) return double_double;
  function QuadDobl_Residuals
              ( file : file_type;
                p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : QuadDobl_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Computes the residuals of the series in s at t,
  --   evaluated in p, with weights for the exponents in w,
  --   in standard double, double double, or quad double precision.
  --   Output is written to file.  Returns the sum of the residuals.

  function Standard_Residuals
              ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                s : Standard_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_float; report : in boolean )
              return double_float;
  function DoblDobl_Residuals
              ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : DoblDobl_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_double; report : in boolean )
              return double_double;
  function QuadDobl_Residuals
              ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : QuadDobl_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : quad_double; report : in boolean )
              return quad_double;

  -- DESCRIPTION :
  --   Computes the residuals of the series in s at t,
  --   evaluated in p, with weights for the exponents in w,
  --   in standard double, double double, or quad double precision.
  --   Output is written to screen if report.
  --   Returns the sum of the residuals.

  procedure Standard_Test
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys );
  procedure DoblDobl_Test
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure QuadDobl_Test
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in standard double, double double, or quad double precision.
  --   The output is written to file.

  procedure Standard_Test 
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                report : in boolean );
  procedure DoblDobl_Test 
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean );
  procedure QuadDobl_Test 
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean );

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in standard double, double double, or quad double precision.
  --   The output is written to screen if report.

  procedure Standard_Read
              ( lp : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 );
  procedure DoblDobl_Read
              ( lp : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 );
  procedure QuadDobl_Read
              ( lp : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system and returns the system,
  --   with coefficients in standard double, double double, or quad
  --   double precision, in lp.
  --   The number of polynomials in returned in nq
  --   and the number of variables is returned in nv.
  --   For the algorithm to apply, we must have nv = nq + 1.

  procedure Standard_Test
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 );
  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 );
  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for the kind of output (to file or to screen)
  --   and then computes the power series for p in standard double,
  --   double double, or quad double precision.
  --   The number of polynomials is in nq and
  --   the number of variables is in nv.

  procedure Standard_Main;
  procedure DoblDobl_Main;
  procedure QuadDobl_Main;

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.
  --   Computations are done in standard double, double double,
  --   or quad double precision.

end Regular_Newton_Puiseux;
