with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;  
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Homotopy_Continuation_Parameters;

package Series_and_Trackers is

-- DESCRIPTION :
--   Path trackers with Newton power series predictors are provided
--   in standard double, double double, or quad double precision.
--   The versions may be silent or verbose.

  procedure Set_Step
              ( t,step : in out double_float;
                maxstep,target : in double_float );

  -- DESCRIPTION :
  --   Sets the step size and the new value for t,
  --   taking into account the maximal step size maxstep
  --   and the value for the target for t.

  procedure Update_Step_Sizes
              ( minsize,maxsize : in out double_float;
                step : in double_float );

  -- DESCRIPTION :
  --   Given in step is a size of a step,
  --   the smallest and largest sizes in minsize and maxsize are updated.

  function Residual_Prediction
              ( sol : Standard_Complex_Vectors.Vector;
                t : double_float ) return double_float;
  function Residual_Prediction
              ( sol : DoblDobl_Complex_Vectors.Vector;
                t : double_float ) return double_float;
  function Residual_Prediction
              ( sol : QuadDobl_Complex_Vectors.Vector;
                t : double_float ) return double_float;

  -- DESCRIPTION :
  --   Given a predicted solution and a value for t,
  --   returns the norm of the evaluated homotopy at the solution.

  -- REQUIRED :
  --   For double, double double, or quad double precision,
  --   the Standard_Homotopy, DoblDobl_Homotopy, or QuadDobl_Homotopy
  --   must have been properly initialized.

  procedure Track_One_Path
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float );
  procedure Track_One_Path
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float );
  procedure Track_One_Path
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom,
  --   in standard double, double double, or quad double precision.
  --   This version remains silent and does not write any output.

  -- ON ENTRY :
  --   hom      a homotopy with series coefficients;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances.

  -- ON RETURN :
  --   sol      solution at the end of the path;
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntfail  is the total number of corrector failures on the path;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path.

  procedure Track_One_Path
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false );
  procedure Track_One_Path
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false );
  procedure Track_One_Path
              ( file : in file_type;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom,
  --   in standard double, double double, or quad double precision.
  --   This version is verbose and writes extra diagnostic output to file.

  -- ON ENTRY :
  --   file     for writing output during the computations;
  --   hom      a homotopy with series coefficients;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances.

  -- ON RETURN :
  --   sol      solution at the end of the path;
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntfail  is the total number of corrector failes on the paths;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path.

-- VERSIONS WITH COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Track_One_Path
              ( file : in file_type;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom,
  --   in standard double precision.
  --   This version is verbose and writes extra diagnostic output to file.

  -- ON ENTRY :
  --   file     for writing output during the computations;
  --   fhm      coefficient-parameter homotopy for efficient evaluation,
  --            the series parameter is the continuation parameter;
  --   fcf      coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances.

  -- ON RETURN :
  --   sol      solution at the end of the path;
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntfail  is the total number of corrector failes on the paths;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path.

  procedure Update_Counters
              ( mincnt,maxcnt : in out natural32; cnt : in natural32 );

  -- DESCRIPTION :
  --   Given the value of a counter in cnt, updates the smallest and
  --   largest value for the counter in mincnt and maxcnt.

  procedure Update_MinMax
              ( smallest,largest : in out double_float;
                minsize,maxsize : in double_float );

  -- DESCRIPTION :
  --   Given the current smallest and largest step sizes,
  --   updates their values with the new minimum and maximum step sizes.

  procedure Track_Many_Paths
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure Track_Many_Paths
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );
  procedure Track_Many_Paths
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters );

  -- DESCRIPTION :
  --   Tracks the paths starting at the solutions in the list sols,
  --   as defined by the homotopy in hom,
  --   in double, double double, or quad double precision.
  --   The procedures are silent and do not write any output.

  -- ON ENTRY :
  --   hom      a homotopy with series coefficients;
  --   sols     start solutions in the homotopy;
  --   pars     values of the parameters and tolerances.

  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                monitor,verbose : in boolean := false );
  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                monitor,verbose : in boolean := false );
  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                monitor,verbose : in boolean := false );

  -- DESCRIPTION :
  --   Tracks the paths starting at the solutions in the list sols,
  --   as defined by the homotopy in hom,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     for output to file;
  --   hom      a homotopy with series coefficients;
  --   sols     start solutions in the homotopy;
  --   pars     values of the parameters and tolerances;
  --   verbose  if true, then extra output is written to file.

  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Write_Path_Statistics
              ( file : in file_type;
                nbrsteps,nbrcorrs,cntfail : in natural32;
                minsize,maxsize : in double_float );

  -- DESCRIPTION :
  --   Writes the path statistics to file.

  -- ON ENTRY :
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntfail  is the number of corrector failures on the path;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path.

  procedure Write_Total_Path_Statistics
              ( file : in file_type;
                minnbrsteps,maxnbrsteps : in natural32;
                minnbrcorrs,maxnbrcorrs : in natural32;
                smallestsize,largestsize : in double_float );

  -- DESCRIPTION :
  --   Writes the statistic for all paths.

  -- ON RETURN :
  --   minnbrsteps is the smallest number of steps on a path;
  --   maxnbrsteps is the largest number of steps on a path;
  --   minnbrcorrs is the smallest number of corrector iterations on a path;
  --   maxnbrcorrs is the largest number of corrector iterations on a path.
  --   smallestsize is the smallest step size on a path;
  --   largestsize is the largest step size on a path.

end Series_and_Trackers;
