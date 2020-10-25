with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with TripDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with PentDobl_Complex_Vectors;
with OctoDobl_Complex_Vectors;
with DecaDobl_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_VecMats;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_VecMats;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_VecMats;
with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;
with DoblDobl_Complex_Circuits;
with TripDobl_Complex_Circuits;
with QuadDobl_Complex_Circuits;
with PentDobl_Complex_Circuits;
with OctoDobl_Complex_Circuits;
with DecaDobl_Complex_Circuits;
with Multitasking;

package Multitasked_Hessian_Circuits is

-- DESCRIPTION :
--   Provides a multitasked implementation of the Hessian criterion,
--   for systems given as sequences of complex circuits,
--   in double, double double, triple double, quad double,
--   and penta double precision.

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
                s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                values : in out TripDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out TripDobl_Complex_VecMats.VecMat;
                e : in out TripDobl_Complex_VecVecs.VecVec;
                yd : in out TripDobl_Complex_VecVecs.VecVec;
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
  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                values : in out PentDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out PentDobl_Complex_VecMats.VecMat;
                e : in out PentDobl_Complex_VecVecs.VecVec;
                yd : in out PentDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                values : in out OctoDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out OctoDobl_Complex_VecMats.VecMat;
                e : in out OctoDobl_Complex_VecVecs.VecVec;
                yd : in out OctoDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                values : in out DecaDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DecaDobl_Complex_VecMats.VecMat;
                e : in out DecaDobl_Complex_VecVecs.VecVec;
                yd : in out DecaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Evaluates all Hessians and computes the singular values with static
  --   load balancing, in double, double double, triple double,
  --   quad double, penta double, octo double, or deca double precision.

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
                s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                values : in out TripDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out TripDobl_Complex_VecMats.VecMat;
                e : in out TripDobl_Complex_VecVecs.VecVec;
                yd : in out TripDobl_Complex_VecVecs.VecVec;
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
  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                values : in out PentDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out PentDobl_Complex_VecMats.VecMat;
                e : in out PentDobl_Complex_VecVecs.VecVec;
                yd : in out PentDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                values : in out OctoDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out OctoDobl_Complex_VecMats.VecMat;
                e : in out OctoDobl_Complex_VecVecs.VecVec;
                yd : in out OctoDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                values : in out DecaDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DecaDobl_Complex_VecMats.VecMat;
                e : in out DecaDobl_Complex_VecVecs.VecVec;
                yd : in out DecaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Evaluates all Hessians ad computes the singular values with dynamic
  --   load balancing, in double, double double, triple double,
  --   quad double, penta double, octo double, or deca double precision.

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
                s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                values : in out TripDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false );
  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false );
  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                values : in out PentDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false );
  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                values : in out OctoDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false );
  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                values : in out DecaDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Evaluates all Hessians of the circuits in s at x
  --   and computes the singular values with nbt tasks,
  --   in double, double double, triple double, quad double,
  --   or penta double precision.
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
  function TripDobl_Distance
              ( values : TripDobl_Complex_VecVecs.VecVec )
              return triple_double;
  function QuadDobl_Distance
              ( values : QuadDobl_Complex_VecVecs.VecVec )
              return quad_double;
  function PentDobl_Distance
              ( values : PentDobl_Complex_VecVecs.VecVec )
              return penta_double;
  function OctoDobl_Distance
              ( values : OctoDobl_Complex_VecVecs.VecVec )
              return octo_double;
  function DecaDobl_Distance
              ( values : DecaDobl_Complex_VecVecs.VecVec )
              return deca_double;

  -- DESCRIPTION :
  --   Returns an estimate to the distance to the nearest solution
  --   based on the smallest singular value of the Jacobian matrix
  --   and the largest singular values of the Hessian matrices,
  --   in double, double double, triple double, quad double,
  --   penta double, octo double, or in deca double precision.

  -- ON ENTRY :
  --   values   values(0) contains the singular values of the Jacobian,
  --            values(k) contains the singular values of the k-th Hessian;

end Multitasked_Hessian_Circuits;
