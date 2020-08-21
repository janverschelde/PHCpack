with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;
with DoblDobl_Complex_Circuits;
with QuadDobl_Complex_Circuits;
with Multitasking;

package Multitasked_Hessian_Circuits is

-- DESCRIPTION :
--   Provides a multitasked implementation of the Hessian criterion,
--   for systems given as sequences of complex circuits,
--   in double, double double, and quad double precision.

  procedure Allocate_Hessian_Spaces
              ( dim : in integer32;
                hrp,hip : out Standard_Floating_VecVecVecs.VecVecVec );

  -- DESCRIPTION :
  --   Allocates space for the real and imaginary parts of Hessian
  --   matrices of dimension dim.
  --   The ranges of hrp and hip are typically 1..nbt,
  --   where nbt equals the number of tasks.

  -- REQUIRED : hrp'range = hip'range.

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Evaluates all Hessians of the circuits in s at x
  --   and computes the singular values with nbt tasks.
  --   This version uses coefficient circuits with complex
  --   coefficients splitted in real and imaginary parts.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        a system of complex circuits;
  --   xr       real parts of the coordinates of the point 
  --            to evaluate the Hessians;
  --   xi       imaginary parts of the coordinates of the point 
  --            to evaluate the Hessians;
  --   values   space allocated for all singular values,
  --            as a vector of range 1..s.dim;
  --   verbose  if verbose, then one line is written for each job,
  --            otherwise, the jobs are preformed without output.

  -- ON RETURN :
  --   values   values(0) contains the singular values of the Jacobian;
  --            values(k) contains the singular values of the k-th Hessian.

  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out Standard_Complex_VecMats.VecMat;
                e : in out Standard_Complex_VecVecs.VecVec;
                yd : in out Standard_Complex_VecVecs.VecVec;
		verbose : in boolean := false );
  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DoblDobl_Complex_VecMats.VecMat;
                e : in out DoblDobl_Complex_VecVecs.VecVec;
                yd : in out DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out QuadDobl_Complex_VecMats.VecMat;
                e : in out QuadDobl_Complex_VecVecs.VecVec;
                yd : in out QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Evaluates all Hessians and computes the singular values with static
  --   load balancing, in double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        a system of complex circuits;
  --   x        coordinates of the point to evaluate the Hessians;
  --   values   space allocated for all singular values,
  --            as a vector of range 1..s.dim;
  --   pwtdone  array of nbt flags to synchronize power table computation,
  --            must all be equal to false on entry;
  --   gradone  array of nbt flags to synchronize gradient computation,
  --            must all be equal to false on entry;
  --   A        vector of range 1..nbt of work space matrices;
  --   U        vector of range 1..nbt of work space matrices;
  --   V        vector of range 1..nbt of work space matrices;
  --   e        vector of range 1..nbt of vector space vectors;
  --   yd       vector of range 1..s.neq of works space for gradients;
  --   verbose  if verbose, then one line is written for each job,
  --            otherwise, the jobs are preformed without output.

  -- ON RETURN :
  --   values   values(0) contains the singular values of the Jacobian;
  --            values(k) contains the singular values of the k-th Hessian;
  --   pwtdone  all equal to true if the power table was needed;
  --   gradone  all equal to true.

  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DoblDobl_Complex_VecMats.VecMat;
                e : in out DoblDobl_Complex_VecVecs.VecVec;
                yd : in out DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out Standard_Complex_VecMats.VecMat;
                e : in out Standard_Complex_VecVecs.VecVec;
                yd : in out Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out QuadDobl_Complex_VecMats.VecMat;
                e : in out QuadDobl_Complex_VecVecs.VecVec;
                yd : in out QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Evaluates all Hessians ad computes the singular values with dynamic
  --   load balancing, in double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        a system of complex circuits;
  --   x        coordinates of the point to evaluate the Hessians;
  --   values   space allocated for all singular values,
  --            as a vector of range 1..s.dim;
  --   pwtdone  array of nbt flags to synchronize power table computation,
  --            must all be equal to false on entry;
  --   gradone  array of nbt flags to synchronize gradient computation,
  --            must all be equal to false on entry;
  --   A        vector of range 1..nbt of work space matrices;
  --   U        vector of range 1..nbt of work space matrices;
  --   V        vector of range 1..nbt of work space matrices;
  --   e        vector of range 1..nbt of vector space vectors;
  --   yd       vector of range 1..s.neq of works space for gradients;
  --   verbose  if verbose, then one line is written for each job,
  --            otherwise, the jobs are preformed without output.

  -- ON RETURN :
  --   values   values(0) contains the singular values of the Jacobian;
  --            values(k) contains the singular values of the k-th Hessian;
  --   pwtdone  all equal to true if the power table was needed;
  --   gradone  all equal to true.

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
  --   and computes the singular values with nbt tasks,
  --   in double, double double, or quad double precision.
  --   Wraps Static_Singular_Values and Dynamic_Singular_Values.

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
