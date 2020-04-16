with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Hyperplane_Convolution_Scaling is

-- DESCRIPTION :
--   After a projective coordinate transformation, the scaling of
--   the solution can be done with a simple adjustment of the last
--   constant in the added linear equation to the circuits.

  procedure Adjust ( cff : in Standard_Complex_VecVecs.VecVec;
                     cst : in Standard_Complex_Vectors.Link_to_Vector;
                     sol : in Standard_Complex_Vectors.Vector );
  procedure Adjust ( cff : in DoblDobl_Complex_VecVecs.VecVec;
                     cst : in DoblDobl_Complex_Vectors.Link_to_Vector;
                     sol : in DoblDobl_Complex_Vectors.Vector );
  procedure Adjust ( cff : in QuadDobl_Complex_VecVecs.VecVec;
                     cst : in QuadDobl_Complex_Vectors.Link_to_Vector;
                     sol : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Adjusts cst(0) so that the solution sol evaluated at the equation
  --   with series coefficients in cff and constant in cst at t = 0
  --   yields zero, relative to the standard double precision,
  --   the double double precision, or the quad double precision.

  -- REQUIRED : cst /= 0 and cff'range = sol'range.

  procedure Adjust_Last_Constant
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in Standard_Complex_Vectors.Vector );
  procedure Adjust_Last_Constant
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in DoblDobl_Complex_Vectors.Vector );
  procedure Adjust_Last_Constant
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Adjusts the last constant in the convolution circuits in hom
  --   so that the solution sol evaluated at the system at t = 0
  --   yields zero, relative to the standard double precision,
  --   the double double precision, or the quad double precision.

  -- REQUIRED : hom is in 1-homogeneous coordinates and the last
  --   equation in hom is linear and sol'range is appropriate.

  procedure Adjust_Last_Radius
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System );
  procedure Adjust_Last_Radius
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System );
  procedure Adjust_Last_Radius
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System );

  -- DESCRIPTION :
  --   Recomputes the last constant in the circuits of abh,
  --   from the corresponding last constant in the circuits of hom.

  procedure Scale_and_Adjust
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector );
  procedure Scale_and_Adjust
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector );
  procedure Scale_and_Adjust
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Scales the coefficients in sol by dividing by the largest component
  --   and the adjusts the constant in the last circuit of hom,
  --   in double, double double, or quad double precision.

end Hyperplane_Convolution_Scaling;
