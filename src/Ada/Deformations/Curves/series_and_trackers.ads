with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;  
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Hessians;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Homotopy_Continuation_Parameters;

package Series_and_Trackers is

-- DESCRIPTION :
--   Path trackers with Newton power series predictors are provided
--   in standard double, double double, or quad double precision.
--   The versions may be silent or verbose.

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

  procedure Update_Ratio_Sum
              ( ratsum : in out double_float; num,den : in natural32 );

  -- DESCRIPTION :
  --   Adds to ratsum the ratio num/den.

  procedure Update_Ratio_Sums
              ( ratsum1,ratsum2,ratsum3 : in out double_float;
                num1,num2,num3,den : in natural32 );

  -- DESCRIPTION :
  --   Adds to ratsum1 the ratio num1/den, to ratsum2 the ratio num2/den,
  --   and to ratsum3 the ratio num3/den.

  procedure Track_Many_Paths
              ( jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure Track_Many_Paths
              ( jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure Track_Many_Paths
              ( jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks the paths starting at the solutions in the list sols,
  --   as defined by the homotopy in hom,
  --   in double, double double, or quad double precision.
  --   The procedures are silent and do not write any output.

  -- ON ENTRY :
  --   jm       Jacobian matrix for the polynomial homotopy;
  --   hs       Hessians for all equations in the polynomial homotopy;
  --   hom      a homotopy with series coefficients;
  --   sols     start solutions in the homotopy;
  --   pars     values of the parameters and tolerances;
  --   mhom     0 for affine coordinates, 1 for 1-homogenization,
  --            and m for m-homogenization;
  --   idz      index representation of the partition of m-homogenization,
  --            in case mhom > 1;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Track_Many_Paths
              ( file : in file_type;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                monitor,verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Track_Many_Paths
              ( file : in file_type;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                monitor,verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Track_Many_Paths
              ( file : in file_type;
                jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                monitor,verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks the paths starting at the solutions in the list sols,
  --   as defined by the homotopy in hom,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     for output to file;
  --   jm       Jacobian matrix for the polynomial homotopy;
  --   hs       Hessians for all equations in the polynomial homotopy;
  --   hom      a homotopy with series coefficients;
  --   sols     start solutions in the homotopy;
  --   pars     values of the parameters and tolerances;
  --   mhom     0 for affine coordinates, 1 for 1-homogenization,
  --            and m for m-homogenization;
  --   idz      index representation of the partition of m-homogenization,
  --            in case mhom > 1;
  --   verbose  if true, then extra output is written to file;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Write_Path_Statistics
              ( file : in file_type;
                nbrsteps,nbrcorrs,cntcut,cntfail : in natural32;
                minsize,maxsize : in double_float;
                cntdstp,cntpstp : in natural32 );

  -- DESCRIPTION :
  --   Writes the path statistics to file.

  -- ON ENTRY :
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntcut   is the number of steps cut by predictor residual;
  --   cntfail  is the number of corrector failures on the path;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path;
  --   cntdstp  counts the number of times the Hessian step was minimal;
  --   cntpstp  counts the number of times the pole step was minimal.

  procedure Write_Total_Path_Statistics
              ( file : in file_type;
                minnbrsteps,maxnbrsteps : in natural32;
                minnbrcorrs,maxnbrcorrs : in natural32;
                smallestsize,largestsize : in double_float;
                ratdstp,ratpstp : in double_float );

  -- DESCRIPTION :
  --   Writes the statistic for all paths.

  -- ON ENTRY :
  --   minnbrsteps is the smallest number of steps on a path;
  --   maxnbrsteps is the largest number of steps on a path;
  --   minnbrcorrs is the smallest number of corrector iterations on a path;
  --   maxnbrcorrs is the largest number of corrector iterations on a path.
  --   smallestsize is the smallest step size on a path;
  --   largestsize is the largest step size on a path;
  --   ratdstp is the average ratio of times Hessian step was minimal;
  --   ratpstp is the average ratio of times pole step was minimal.

end Series_and_Trackers;
