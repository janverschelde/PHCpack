with text_io;                           use text_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;    use Partitions_of_Sets_of_Unknowns;

package Intrinsic_Diagonal_Solvers is

-- DESCRIPTION :
--   This package offers tools to solve polynomial systems incrementally
--   by intrinsic diagonal homotopies.

  procedure Witness_Points
              ( n : in natural; p,q : in Poly; b2 : out Vector;
                w2 : out VecVec; sols : out Array_of_Solution_Lists;
                fail : out boolean );
  procedure Witness_Points 
              ( file : in file_type;
                n : in natural; p,q : in Poly; b2 : out Vector;
                w2 : out VecVec; sols : out Array_of_Solution_Lists;
                fail : out boolean );

  -- DESCRIPTION :
  --   Computes witness points on all components of p and q.

  -- ON ENTRY :
  --   file     for intermediate results or diagnostics;
  --   n        dimension of the ambient space = #variables;
  --   p        first polynomial in n variables;
  --   q        second polynomial in n variables.

  -- ON RETURN :
  --   b2       double of offset vector of affine 2-plane;
  --   w2       directions of affine 2-plane in double format;
  --   sols(1)  collects the common factors whose witness points
  --            lie on the line with direction in w2(1);
  --   sols(2)  linear combinations of vectors in w2, representing
  --            the witness points cut out by an affine 2-plane;
  --   fail     true if witness point calculation of the polynomials
  --            p and q failed, false otherwise.

  procedure Witness_Points 
              ( n : in natural; p,q : in Poly; z : in Partition;
                b2 : out Vector; w2 : out Array_of_VecVecs;
                sols : out Array_of_Solution_Lists );
  procedure Witness_Points 
              ( file : in file_type;
                n : in natural; p,q : in Poly; z : in Partition;
                b2 : out Vector; w2 : out Array_of_VecVecs;
                sols : out Array_of_Solution_Lists );

  -- DESCRIPTION :
  --   Computes witness points on the intersection of two hypersurfaces,
  --   respecting their multi-homogeneous structure.

  -- ON ENTRY :
  --   file     for intermediate results or diagnostics;
  --   n        dimension of the ambient space = #variables;
  --   p        first polynomial in n variables;
  --   q        second polynomial in n variables.

  -- ON RETURN :
  --   b2       double of offset vector of affine 2-plane;
  --   w2       direction vectors of affine 2-planes, the i-th plane is
  --            defined as b2 + t1*w2(1)(i) + t2*w2(2)(i);
  --   sols     the i-th solution list contains values for t1 and t2 
  --            to represent the witness points cut out by the i-th 
  --            multi-homogeneous affine 2-plane. 

  procedure Filter_by_Evaluation
              ( file : in file_type; n : in natural; p : in Poly;
                b : in Vector; v : in VecVec; tol : in double_float;
		sols : in Solution_List; on_p,off_p : out Solution_List );

  -- DESCRIPTION :
  --   Determines whether the polynomial p is superfluous with respect to
  --   the k-dimensional solution component given by its witness points.

  -- ON ENTRY :
  --   file     to write intermediate results and diagnostics;
  --   n        dimension of the ambient space;
  --   p        polynomial in n variables to be evaluated;
  --   b        offset vector of the general affine k-plane;
  --   v        k directions vectors of the general affine k-plane;
  --   sols     coefficients in the linear combinations of the directions;
  --   tol      tolerance to decide when a number is zero.

  -- ON RETURN :
  --   on_p     solutions of sols whose value in p is smaller than tol;
  --   off_p    solutions of sols whose value in p is larger than tol.

  procedure Shift_to_End ( p : in out Poly_Sys; ind : in natural );

  -- DESCRIPTION :
  --   The polynomial p(ind) is moved to p(p'last) and all polynomials
  --   at positions ind+1..p'last shift up one position.

  procedure Witness_Points
              ( file : in file_type; n : in natural;
                p : in out Poly_Sys; k,ind : in out natural;
                b2 : in out Vector; w : in out VecVec;
		sols : in out Solution_List; fail : out boolean );

  -- DESCRIPTION :
  --   Intersection of a set of co-dimension k with p(ind) = 0.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of variables in the equations of p;
  --   p        polynomial system in n variables;
  --   k        co-dimension of the current solution set;
  --   ind      index of new equation in p to intersect with;
  --   b2       doubled basis vector used in diagonal homotopy;
  --   w        contains directions spanning affine k-plane;
  --   sols     spanning vectors for current solution set.

  -- ON RETURN :
  --   p        if p(ind) is superfluous, it is shifted to the end;
  --   k        equals the previous value for k if p(ind) is superfluous,
  --            otherwise k has been augmented by one to equal co-dimension;
  --   ind      equals the previous value for ind if p(ind) is superfluous,
  --            otherwise ind has been augmented by one;
  --   b2       new doubled basis vector for new affine plane;
  --   w        directions extended with direction at position k+1;
  --   sols     if p(ind) is not superfluous, then sols contains spanning
  --            vectors of solution set of co-dimension k+1;
  --   fail     true if witness points on p(ind) = 0 could not be computed,
  --            false otherwise.

  procedure Total_Degree_Hypersurface_Solver
              ( file : in file_type; n : in natural; p : in out Poly_Sys;
                b : out Vector; w : out VecVec;
                sols : out Array_of_Solution_Lists );

  -- DESCRIPTION :
  --   Solves p incrementally, computing witness points on all components.

  -- ON ENTRY :
  --   file     to write diagnostics and intermediate results;
  --   n        dimension of the ambient space, number of variables in p;
  --   p        polynomial system in n variables;

  -- ON RETURN :
  --   p        superfluous polynomials are shifted towards the end;
  --   b        random offset vector for the affine planes;
  --   w        contains directions for the affine planes;
  --   sols     sols(k) contains witness points on solution components
  --            of co-dimension k on affine planes spanned by offset b
  --            and directions in w(1..k).

  procedure Multihomogeneous_Hypersurface_Solver
              ( file : in file_type; n : in natural; z : in Partition;
                p : in out Poly_Sys );

  -- DESCRIPTION :
  --   Computes witness points on pure dimenional components, solving
  --   the polynomial system p incrementally.  On return, superfluous
  --   polynomials in p have been shifted to the end of the system.
  --   Respects the multi-homogeneous structure defined by z.

end Intrinsic_Diagonal_Solvers;
