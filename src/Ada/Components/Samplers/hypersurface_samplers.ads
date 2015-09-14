with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Partitions_of_Sets_of_Unknowns;    use Partitions_of_Sets_of_Unknowns;

package Hypersurface_Samplers is

-- DESCRIPTION :
--   A generic point on a hypersurface is the intersection of the hypersurface
--   with a general line.  The routines in this package produce new generic
--   points by varying the general line.

  function Random_Multihomogeneous_Directions
                ( n : natural32; z : Partition )
                return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of random multi-homogeneous directions:
  --   the components of the i-th n-vector are nonzero if and only if
  --   the corresponding variable belongs to the i-th set in z.

  function Random_Set_Structure_Directions
                ( n,i : natural32 )
                return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of random directions according to the set structure
  --   for the i-th polynomial: the components of the j-th direction on
  --   return are nonzero if and only if the corresponding variable belongs
  --   to the j-th set of the set structure for the i-th polynomial.

  -- REQUIRED :
  --   The i-th component of the set structure has been built.

  procedure Generic_Points
                ( p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  z : in Partition;
                  b : in Standard_Complex_Vectors.Vector;
                  v : in Standard_Complex_VecVecs.VecVec;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_VecVecs.VecVec;
                  fail : out boolean );
  procedure Generic_Points
                ( file : in file_type;
                  p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  z : in Partition;
                  b : in Standard_Complex_Vectors.Vector;
                  v : in Standard_Complex_VecVecs.VecVec;
                  eps : in double_float; maxit : in natural32;
                  t : out Standard_Complex_VecVecs.VecVec;
                  fail : out boolean );

  -- DESCRIPTION :
  --   Computes generic points on the hypersurface p(x) = 0 intersected
  --   with lines which respect the multihomogeneous structure defined
  --   by the partition z.

  -- ON ENTRY :
  --   file       for intermediate output and diagnostics;
  --   p          polynomial in n variables;
  --   ep         evaluable Horner form of the polynomial p;
  --   z          partition of the set of unknowns of p;
  --   b          basis vector for all lines intersecting p(x) = 0;
  --   v          as many directions as z'length, the i-th vector has only
  --              nonzero entries corresponding to the variables in z(i),
  --              we intersect p(x) = 0 with z'length lines b + v(i);
  --   eps        accuracy requirement on the roots;
  --   maxit      maximal number of iterations allowed.

  -- ON RETURN :
  --   t          the i-th vector contains the roots on the i-th line;
  --   fail       true if one of the root computations failed,
  --              false if all roots are computed with the required accuracy.

end Hypersurface_Samplers;
