with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Circuits;
with DoblDobl_Complex_Circuits;
with QuadDobl_Complex_Circuits;

package Multitasked_Hessian_Circuits is

-- DESCRIPTION :
--   Provides a multitasked implementation of the Hessian criterion,
--   for systems given as sequences of complex circuits,
--   in double, double double, and quad double precision.

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return Standard_Complex_VecVecs.VecVec;
  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return DoblDobl_Complex_VecVecs.VecVec;
  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns an array of range neqstart..neq,
  --   with allocated vectors of range dimstart..dim,
  --   in double, double double, or quad double precision.

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false );
  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false );
  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Evaluates all Hessians of the circuits in s at x
  --   and computes the singular values with nbt tasks.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        a system of complex circuits;
  --   x        coordinates of the point to evaluate the Hessians;
  --   values   space allocated for all singular values,
  --            as a vector of range 1..s.dim;
  --   static   flag to to apply static load balancing,
  --            by default dynamic load balancing is applied;
  --   verbose  if verbose, then one line is written for each job,
  --            otherwise, the jobs are preformed without output.

  -- ON RETURN :
  --   values   values(0) contains the singular values of the Jacobian;
  --            values(k) contains the singular values of the k-th Hessian.

  function Standard_Distance
              ( values : Standard_Complex_VecVecs.VecVec )
              return double_float;
  function DoblDobl_Distance
              ( values : DoblDobl_Complex_VecVecs.VecVec )
              return double_double;
  function QuadDobl_Distance
              ( values : QuadDobl_Complex_VecVecs.VecVec )
              return quad_double;

  -- DESCRIPTION :
  --   Returns an estimate to the distance to the nearest solution
  --   based on the smallest singular value of the Jacobian matrix
  --   and the largest singular values of the Hessian matrices,
  --   in double, double double, or in quad double precision.

  -- ON ENTRY :
  --   values   values(0) contains the singular values of the Jacobian,
  --            values(k) contains the singular values of the k-th Hessian;

end Multitasked_Hessian_Circuits;
