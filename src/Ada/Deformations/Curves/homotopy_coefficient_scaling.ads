with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with Standard_CSeries_Poly_SysFun;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with DoblDobl_CSeries_Poly_SysFun;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with QuadDobl_CSeries_Poly_SysFun;

package Homotopy_Coefficient_Scaling is

-- DESCRIPTION :
--   The procedures in this package recompute the coefficients of the
--   coefficient and the series homotopies.

  procedure Last_Coefficients
              ( file : in file_type;
                fcf : in Standard_Complex_Series_Vectors.Link_to_Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number );
  procedure Last_Coefficients
              ( file : in file_type;
                fcf : in DoblDobl_Complex_Series_Vectors.Link_to_Vector;
                t : in double_double;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Last_Coefficients
              ( file : in file_type;
                fcf : in QuadDobl_Complex_Series_Vectors.Link_to_Vector;
                t : in quad_double;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Checks the last coefficients in the linear equation added
  --   by the projective transformation and writes diagnostics to file,
  --   in double, double double, and quad double precision.
  --   The input fcf are the current coefficient vectors of the homotopy.
  --   This is an exploratory testing procedure to verify that the last
  --   coefficients in fcf can be recomputed from the coefficients in the
  --   defined coefficient homotopy.
 
  procedure Scale_Solution_Coefficients
              ( hcf : in Standard_Complex_Series_VecVecs.VecVec;
                sol : in out Standard_Complex_Vectors.Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number );
  procedure Scale_Solution_Coefficients
              ( hcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                t : in double_double;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Scale_Solution_Coefficients
              ( hcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                t : in quad_double;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Scales the solution and the coefficients in the homotopy,
  --   in double, double double, and quad double precision.
  --   No output or extra verifications are computed.

  -- ON ENTRY :
  --   fhm      coefficient version of the series version for evaluation;
  --   hcf      coefficient vectors for fhm;
  --   sol      a solution;
  --   t        current value of the continuation parameter;
  --   gamma    the random gamma constant in the artificial-parameter homotopy.

  -- ON RETURN :
  --   sol      scaled solution, each component is divided by the magnitude
  --            of the largest coordinate in sol.

  procedure Scale_Solution_Coefficients
              ( file : in file_type;
                fhm : Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                hcf : in Standard_Complex_Series_VecVecs.VecVec;
                sol : in out Standard_Complex_Vectors.Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false );
  procedure Scale_Solution_Coefficients
              ( file : in file_type;
                fhm : DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                hcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                t : in double_double;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false );
  procedure Scale_Solution_Coefficients
              ( file : in file_type;
                fhm : QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                hcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                t : in quad_double;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Scales the solution and the coefficients in the homotopy,
  --   in double, double double, and quad double precision.
  --   Extra verifications are computed and diagnostics are written,
  --   mainly for debugging purposes.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   fhm      coefficient version of the series version for evaluation;
  --   hcf      coefficient vectors for fhm;
  --   sol      a solution;
  --   t        current value of the continuation parameter;
  --   gamma    the random gamma constant in the artificial-parameter homotopy;
  --   verbose  for extra verification output.

  -- ON RETURN :
  --   sol      scaled solution, each component is divided by the magnitude
  --            of the largest coordinate in sol.

-- MULTI-HOMOGENEOUS VERSIONS :

  procedure Last_Coefficients
              ( file : in file_type;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                m : in natural32 );
  procedure Last_Coefficients
              ( file : in file_type;
                fcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                t : in double_double;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number;
                m : in natural32 );
  procedure Last_Coefficients
              ( file : in file_type;
                fcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                t : in quad_double;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number;
                m : in natural32 );

  -- DESCRIPTION :
  --   Checks the last coefficients in the linear equations added
  --   by the m-projective transformation and writes diagnostics to file,
  --   in double, double double, and quad double precision.
  --   This is an exploratory testing procedure to verify that the last
  --   coefficients in fcf can be recomputed from the coefficients in the
  --   defined coefficient homotopy.

  -- ON ENTRY :
  --   file     file opened for output;
  --   fcf      current coefficient vectors of the homotopy;
  --   t        current value of the homotopy continuation parameter;
  --   gamma    gamma constant in the homotopy;
  --   m        number of sets in the partition of the variables,
  --            the variables added in the m-homogenization are always
  --            at the last m positions, in increasing order.

  procedure Scale_Solution_Coefficients
              ( hcf : in Standard_Complex_Series_VecVecs.VecVec;
                sol : in out Standard_Complex_Vectors.Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                m : in natural32;
                z : in Standard_Natural_Vectors.Vector );
  procedure Scale_Solution_Coefficients
              ( hcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                t : in double_double;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number;
                m : in natural32;
                z : in Standard_Natural_Vectors.Vector );
  procedure Scale_Solution_Coefficients
              ( hcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                t : in quad_double;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number;
                m : in natural32;
                z : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Scales the solution and the coefficients in the homotopy,
  --   after multi-projective coordinate transformations, for m > 1,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   hcf      coefficient vectors for fhm;
  --   sol      a solution;
  --   t        current value of the continuation parameter;
  --   gamma    the random gamma constant in the artificial-parameter homotopy;
  --   m        number of sets in the partition of the variables,
  --            the variables added in the m-homogenization are always
  --            at the last m positions, in increasing order;
  --   z        the index representation of the partition.

  -- ON RETURN :
  --   sol      scaled solution, each component is divided by the magnitude
  --            of the largest coordinate in sol.

  procedure Scale_Solution_Coefficients
              ( file : in file_type;
                fhm : Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                hcf : in Standard_Complex_Series_VecVecs.VecVec;
                sol : in out Standard_Complex_Vectors.Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                m : in natural32;
                z : in Standard_Natural_Vectors.Vector;
                verbose : in boolean := false );
  procedure Scale_Solution_Coefficients
              ( file : in file_type;
                fhm : DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                hcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                t : in double_double;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number;
                m : in natural32;
                z : in Standard_Natural_Vectors.Vector;
                verbose : in boolean := false );
  procedure Scale_Solution_Coefficients
              ( file : in file_type;
                fhm : QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                hcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                t : in quad_double;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number;
                m : in natural32;
                z : in Standard_Natural_Vectors.Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Scales the solution and the coefficients in the homotopy,
  --   after multi-projective coordinate transformations, for m > 1,
  --   in double, double double, and quad double precision.
  --   Extra verifications are computed and diagnostics are written,
  --   mainly for debugging purposes.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   fhm      coefficient version of the series version for evaluation;
  --   hcf      coefficient vectors for fhm;
  --   sol      a solution;
  --   t        current value of the continuation parameter;
  --   gamma    the random gamma constant in the artificial-parameter homotopy;
  --   m        number of sets in the partition of the variables,
  --            the variables added in the m-homogenization are always
  --            at the last m positions, in increasing order;
  --   z        the index representation of the partition;
  --   verbose  for extra verification output.

  -- ON RETURN :
  --   sol      scaled solution, each component is divided by the magnitude
  --            of the largest coordinate in sol.

end Homotopy_Coefficient_Scaling;
