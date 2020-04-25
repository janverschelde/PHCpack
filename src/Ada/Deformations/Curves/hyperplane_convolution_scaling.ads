with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
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

-- 1-HOMOGENIZATION :

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

-- MULTI-HOMOGENIZATION :

  procedure Adjust ( cff : in Standard_Complex_VecVecs.VecVec;
                     cst : in Standard_Complex_Vectors.Link_to_Vector;
                     sol : in Standard_Complex_Vectors.Vector;
                     idz : in Standard_Natural_Vectors.Link_to_Vector;
                     m,i : in integer32 );
  procedure Adjust ( cff : in DoblDobl_Complex_VecVecs.VecVec;
                     cst : in DoblDobl_Complex_Vectors.Link_to_Vector;
                     sol : in DoblDobl_Complex_Vectors.Vector;
                     idz : in Standard_Natural_Vectors.Link_to_Vector;
                     m,i : in integer32 );
  procedure Adjust ( cff : in QuadDobl_Complex_VecVecs.VecVec;
                     cst : in QuadDobl_Complex_Vectors.Link_to_Vector;
                     sol : in QuadDobl_Complex_Vectors.Vector;
                     idz : in Standard_Natural_Vectors.Link_to_Vector;
                     m,i : in integer32 );

  -- DESCRIPTION :
  --   Adjusts cst(0) so that the solution sol evaluated at the 
  --   last i-th added linear equation with series coefficients in cff
  --   and constant in cst at t = 0 yields zero,
  --   relative to double, double double, or quad double precision.

  -- REQUIRED : cst /= 0 and cff'range = sol'range.

  -- ON ENTRY :
  --   cff      coefficient vector of the series coefficients
  --            of a linear equation;
  --   sol      a solution;
  --   idz      index representation of the partition of the variables:
  --            idz(k) is the index of the set to which the k-th variable
  --            belongs, for k in 1..m, idz'range = 1..dim,
  --            where dim equals the number of variables.
  --   m        number of sets in the partition of the variables,
  --            the variables added in the m-homogenization are always
  --            at the last m positions, in increasing order;
  --   i        index of the added homogeneous variable,
  --            i must be in the range 1..m.

  -- ON RETURN :
  --   cst(0)   adjusted constant so the solution evaluates to zero.

  procedure Adjust_Last_Radii
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                m : in integer32 );
  procedure Adjust_Last_Radii
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                m : in integer32 );
  procedure Adjust_Last_Radii
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                m : in integer32 );

  -- DESCRIPTION :
  --   Recomputes the constants in the last m circuits of abh,
  --   from the corresponding constants in the last m circuits of hom.

  procedure Scale_and_Adjust
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                m : in integer32 );
  procedure Scale_and_Adjust
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                m : in integer32 );
  procedure Scale_and_Adjust
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                m : in integer32 );

  -- DESCRIPTION :
  --   Scales the solution in sol with respect to the multi-homogenization
  --   defined by idz and m and adjusts the constants in the last m circuits
  --   in the homotopy hom so the scaled solutions satisfy the system,
  --   relative to the double, double double, or quad double precision.

  -- ON ENTRY :
  --   hom      convolution circuits for a homotopy;
  --   sol      a solution;
  --   idz      index representation of the partition of the variables:
  --            idz(k) is the index of the set to which the k-th variable
  --            belongs, for k in 1..m, idz'range = 1..dim,
  --            where dim equals the number of variables.
  --   m        number of sets in the partition of the variables,
  --            the variables added in the m-homogenization are always
  --            at the last m positions, in increasing order;
  --   i        index of the added homogeneous variable,
  --            i must be in the range 1..m.

  -- ON RETURN :
  --   hom      adjusted constant coefficients of the last m circuits;
  --   sol      solution scaled with respect to the multi-homogenization.

end Hyperplane_Convolution_Scaling;
