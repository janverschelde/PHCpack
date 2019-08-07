with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with Standard_CSeries_Poly_SysFun;

package Homotopy_Coefficient_Scaling is

-- DESCRIPTION :
--   The procedures in this package recompute the coefficients of the
--   coefficient and the series homotopies.

  procedure Last_Coefficients
              ( file : in file_type;
                fcf : in Standard_Complex_Series_Vectors.Link_to_Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Checks the last coefficients in the linear equation added
  --   by the projective transformation and writes diagnostics to file.
  --   The input fcf are the current coefficient vectors of the homotopy.
  --   This is an exploratory testing procedure.
 
  procedure Scale_Solution_Coefficients
              ( hcf : in Standard_Complex_Series_VecVecs.VecVec;
                sol : in out Standard_Complex_Vectors.Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number );
  procedure Scale_Solution_Coefficients
              ( file : in file_type;
                fhm : Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                hcf : in Standard_Complex_Series_VecVecs.VecVec;
                sol : in out Standard_Complex_Vectors.Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Scales the solution and the coefficients in the homotopy.

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

end Homotopy_Coefficient_Scaling;
