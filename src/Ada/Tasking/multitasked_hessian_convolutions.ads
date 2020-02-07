with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_VecMats;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Multitasked_Hessian_Convolutions is

-- DESCRIPTION :
--   Provides a multitasked implementation of the Hessian criterion,
--   for systems given as convolution circuits,
--   in double, double double, and quad double precision.

  function Allocate ( nbr,dim : integer32 )
                    return Standard_Complex_VecMats.VecMat;
  function Allocate ( nbr,dim : integer32 )
                    return DoblDobl_Complex_VecMats.VecMat;
  function Allocate ( nbr,dim : integer32 )
                    return QuadDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Returns an array of range 1..nbr, with allocated dim-by-dim matrices,
  --   in double, double double, or quad double precision.

  function Allocate ( nbr,dim : integer32 )
                    return Standard_Complex_VecVecs.VecVec;
  function Allocate ( nbr,dim : integer32 )
                    return DoblDobl_Complex_VecVecs.VecVec;
  function Allocate ( nbr,dim : integer32 )
                    return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns an array of range 1..nbr, with allocated vectors of
  --   range 1..dim, in double, double double, or quad double preeision.

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_Vectors.Vector;
                jmsvls : out Standard_Complex_Vectors.Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := true );
  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Vector;
                jmsvls : out DoblDobl_Complex_Vectors.Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true );
  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Vector;
                jmsvls : out QuadDobl_Complex_Vectors.Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Evaluates all Hessians of the circuits in s at x
  --   and computes the singular values with nbt tasks.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        a system of convolution circuits;
  --   x        coordinates of the point to evaluate the Hessians;
  --   values   space allocated for all singular values,
  --            as a vector of range 1..s.dim;
  --   verbose  if verbose, then one line is written for each job,
  --            otherwise, the jobs are preformed without output.

  -- ON RETURN :
  --   jmsvls   singular values of the Jacobian matrix;
  --   values   values(k) contains the singular values of the k-th Hessian.

end Multitasked_Hessian_Convolutions;
