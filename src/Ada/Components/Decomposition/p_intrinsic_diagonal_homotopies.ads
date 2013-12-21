with text_io;                           use text_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;    use Partitions_of_Sets_of_Unknowns;

package P_Intrinsic_Diagonal_Homotopies is

-- DESCRIPTION :
--   Tools to implement the intrinsic version of diagonal homotopies.
--   With diagonal homotopies we compute witness points on positive
--   dimensional solution components.  Many routines have 4 versions:
--   with output or not, and with evaluable polynomials or not.
--   The versions without evaluable forms of the polynomials are older.
--   The purpose of the newer versions which have more parameters is to
--   avoid a repeated creation of those evaluable forms.

-- UTILITIES :

  function Product ( a,b : Vector ) return Solution_List;

  -- DESCRIPTION :
  --   Returns the product of the two vectors a and b as a list of solutions.

  function Product ( sols : Solution_List; a : Vector ) return Solution_List;

  -- DESCRIPTION :
  --   Returns the product of the solution list with the vector a.

  function On_Hypersurface
             ( p : Poly; b,v : Vector; t : Complex_Number;
               tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the point b + t*v evaluates in p to a number
  --   whose absolute value is less than tol.

  function On_Hypersurface
             ( p : Poly; b,v,t : Vector; tol : double_float )
             return Solution_List;

  -- DESCRIPTION :
  --   Returns a solution list with one vector:  all witness points
  --   on the line x(t) = bp + t*vp, which lie on p, within the given
  --   tolerance tol.

  function Remove_on_Hypersurface
             ( p : Poly; b,v,t : Vector; tol : double_float ) return Vector;

  -- DESCRIPTION :
  --   On return are those points in t on the line x(t) = b + t*v which
  --   do not evaluate in p to a number whose absolute value < tol.

-- WITNESS POINTS ON ONE HYPERSURFACE :

  procedure Hypersurface_Witness_Points
              ( n,dp : in natural; p : in Poly; ep : in Eval_Poly;
                b,v : in Vector; tp : out Vector; fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( n,dp : in natural; p : in Poly; b,v : in Vector;
                tp : out Vector; fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp : in natural; p : in Poly; ep : in Eval_Poly;
                b,v : in Vector; tp : out Vector; fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp : in natural; p : in Poly; b,v : in Vector;
                tp : out Vector; fail : out boolean );

  -- DESCRIPTION :
  --   Computes witness points on the surface defined by p(x) = 0
  --   and on the affine line x = b + t*v.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of unknowns in p;
  --   dp       degree of the polynomial p;
  --   p        a polynomial in n variables of degree dp;
  --   ep       evaluable form of the polynomial p;
  --   b        offset vector of an affine line;
  --   v        direction of an affine line.

  -- ON RETURN :
  --   tp       values representing dp points on the affine line b+t*v,
  --            which are witness points on the hypersurface p(x) = 0;
  --   fail     true if the desired accuracy was not reached for one of
  --            the hypersurfaces, false otherwise.

-- WITNESS POINTS ON TWO HYPERSURFACES :

  procedure Hypersurface_Witness_Points
              ( n,dp,dq : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                bp,vp,bq,vq : in Vector; tp,tq : out Vector;
                fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( n,dp,dq : in natural; p,q : in Poly;
                bp,vp,bq,vq : in Vector; tp,tq : out Vector;
                fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp,dq : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                bp,vp,bq,vq : in Vector; tp,tq : out Vector;
                fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp,dq : in natural; p,q : in Poly;
                bp,vp,bq,vq : in Vector; tp,tq : out Vector;
                fail : out boolean );

  -- DESCRIPTION :
  --   Computes a set of witness points on the hypersurfaces defined
  --   by the multivariate polynomials p and q.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of unknowns in p and q;
  --   dp       degree of the polynomial p;
  --   dq       degree of the polynomial q;
  --   p        a polynomial in n variables of degree dp;
  --   q        a polynomial in n variables of degree dq;
  --   ep       evaluable form of the polynomial p;
  --   eq       evaluable form of the polynomial q;
  --   bp       offset vector of an affine line bp + t*vp;
  --   vp       direction of an affine line bp + t*vp.
  --   bq       offset vector of an affine line bq + t*vq;
  --   vq       direction of an affine line bq + t*vq.

  -- ON RETURN :
  --   tp       values representing dp points on the affine line bp+t*vp,
  --            which are witness points on the hypersurface p(x) = 0;
  --   tq       values representing dq points on the affine line bq+t*vq,
  --            which are witness points on the hypersurface q(x) = 0;
  --   fail     true if the desired accuracy was not reached for one of
  --            the hypersurfaces, false otherwise.

-- WITNESS POINTS ON ONE MULTI-HOMOGENEOUS HYPERSURFACE :

  procedure Hypersurface_Witness_Points
              ( n,dp : in natural; p : in Poly; ep : in Eval_Poly;
                z : in Partition; b : in Vector; v : in VecVec;
                tp : out VecVec; fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( n : in natural; p : in Poly; z : in Partition;
                b : in Vector; v : in VecVec; tp : out VecVec;
                fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp : in natural; p : in Poly; ep : in Eval_Poly;
                z : in Partition; b : in Vector; v : in VecVec;
                tp : out VecVec; fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n : in natural; p : in Poly; z : in Partition;
                b : in Vector; v : in VecVec; tp : out VecVec;
                fail : out boolean );

  -- DESCRIPTION :
  --   Computes witness points on the surface defined by p(x) = 0, on
  --   lines x(t) = b + t*v(i), respecting a multi-homogeneous structure.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of unknowns in p;
  --   dp       degree of the polynomial p;
  --   p        a polynomial in n variables;
  --   ep       evaluable form of p;
  --   z        partitions of the sets of unknowns, 
  --            defining a multi-homogeneous structure;
  --   b        offset vector of an affine line;
  --   v(i)     directions of the affine lines x(t) = b + t*v(i).

  -- ON RETURN :
  --   tp       the values in tp(i) define witness points on the affine
  --            multi-homogeneous lines b + v(i) intersecting p(x) = 0;
  --   fail     true if the desired accuracy was not reached for one of
  --            the hypersurfaces, false otherwise.

-- WITNESS POINTS ON TWO MULTI-HOMOGENEOUS HYPERSURFACES :

  procedure Hypersurface_Witness_Points
              ( n,dp,dq : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                z : in Partition; bp,bq : in Vector; vp,vq : in VecVec;
                tp,tq : out VecVec; fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( n : in natural; p,q : in Poly;
                z : in Partition; bp,bq : in Vector; vp,vq : in VecVec;
                tp,tq : out VecVec; fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n,dp,dq : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                z : in Partition; bp,bq : in Vector; vp,vq : in VecVec;
                tp,tq : out VecVec; fail : out boolean );
  procedure Hypersurface_Witness_Points
              ( file : in file_type;
                n : in natural; p,q : in Poly;
                z : in Partition; bp,bq : in Vector; vp,vq : in VecVec;
                tp,tq : out VecVec; fail : out boolean );

  -- DESCRIPTION :
  --   Computation of witness points on the hypersurfaces p(x) = 0 = q(x),
  --   with respect to a multi-homogeneous structure.
  --   This means that the support of the directions in vp and vq
  --   respects the partition of the sets of unknowns.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of variables in the polynomials;
  --   dp,dq    degrees of p and q;
  --   p,q      polynomials in n variables with complex coefficients;
  --   ep,eq    evaluable forms of polynomials p and q;
  --   z        defines a multi-homogeneous structure;
  --   bp       offset vector for all the affine lines to cut p with;
  --   vp       direction vectors for the affine lines bp + t*vp(i);
  --   bq       offset vector for all the affine lines to cut q with;
  --   vq       direction vectors for the affine lines bq + t*vq(i).

  -- ON RETURN :
  --   tp       the values in tp(i) define witness points on the affine
  --            special lines bp + t*vp(i) intersecting p(x) = 0;
  --   tq       the values in tq(i) define witness points on the affine
  --            special lines bq + t*vq(i) intersecting q(x) = 0;
  --   fail     is true if the desired accuracy is not met for all lines,
  --            otherwise this value is false.  

-- WITNESS POINTS ON THE INTERSECTION OF TWO HYPERSURFACES :

  procedure Diagonal_Homotopy
              ( n : in natural; p,q : in Poly;
                bp,vp,bq,vq,tp,tq : in Vector;
                b2 : out Vector; w2 : out VecVec; sols : out Solution_List );
  procedure Diagonal_Homotopy
              ( n : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                bp,vp,bq,vq,tp,tq : in Vector;
                b2 : out Vector; w2 : out VecVec; sols : out Solution_List );
  procedure Diagonal_Homotopy
              ( file : in file_type;
                n : in natural; p,q : in Poly;
                bp,vp,bq,vq,tp,tq : in Vector;
                b2 : out Vector; w2 : out VecVec; sols : out Solution_List );
  procedure Diagonal_Homotopy
              ( file : in file_type;
                n : in natural; p,q : in Poly; ep,eq : in Eval_Poly;
                bp,vp,bq,vq,tp,tq : in Vector;
                b2 : out Vector; w2 : out VecVec; sols : out Solution_List );

  -- DESCRIPTION :
  --   Applies a diagonal homotopy to compute witness points at components
  --   of dimension n-2.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of unknowns in p and q, dimension of ambient space;
  --   p,q      two polynomials in n variables;
  --   ep,eq    evaluable forms of polynomials p and q;
  --   bp       offset vector of an affine line used to cut p with;
  --   vp       direction vector of an affine line used to cut p with;
  --   bq       offset vector of an affine line used to cut q with;
  --   vq       direction vector of an affine line used to cut q with;
  --   tp       witness points on p(x=bp+t*vp) = 0, on an affine line;
  --   tq       witness points on q(x=bq+t*vq) = 0, on an affine line.

  -- ON RETURN :
  --   b2       double of the offset vector of affine 2-plane;
  --   w2       two directions of an affine 2-plane in 2n-space;
  --   sols     witness points on the components of dimension n-2,
  --            on the affine 2-plane with offset b2 and directions w2.

-- GENERAL WITNESS POINTS :

  procedure Diagonal_Homotopy
              ( n : in natural; p : in Poly_Sys; q : in Poly;
                bp : in Vector; vp : in VecVec; tp : in Solution_List;
                bq : in Vector; vq : in Vector; tq : in Vector;
                b2 : out Vector; w2 : out VecVec; sols : out Solution_List );
  procedure Diagonal_Homotopy
              ( file : in file_type;
                n : in natural; p : in Poly_Sys; q : in Poly;
                bp : in Vector; vp : in VecVec; tp : in Solution_List;
                bq : in Vector; vq : in Vector; tq : in Vector;
                b2 : out Vector; w2 : out VecVec; sols : out Solution_List );

  -- DESCRIPTION :
  --   Applies a diagonal homotopy to compute witness points at the
  --   solution components of p intersected with the hypersurface q.

  -- ON ENTRY :
  --   file     for intermediate results and diagnostics;
  --   n        number of unknowns in p and q, dimension of ambient space;
  --   p        polynomial system in n variables;
  --   q        polynomial in n variables;
  --   bp       offset vector for affine plane which intersects p(x) = 0;
  --   vp       directions of the affine plane which intersects p(x) = 0;
  --   tp       witness points, solutions of p(x) = 0, on affine plane;
  --   bq       offset vector of affine line to cut q(x) = 0;
  --   vq       direction of affine line to cut q(x) = 0;
  --   tq       witness points on hypersurface q(x) = 0, on affine line.

  -- ON RETURN :
  --   b2       double of offset vector of affine plane used to
  --            intersect the solutions of p(x) = 0 with q(x) = 0;
  --   w2       directions of an affine plane in doubled format, used to
  --            intersect the solutions of p(x) = 0 with q(x) = 0;
  --   sols     witness points on the components of p and q on the 
  --            affine plane with offset b2 and directions w2.

-- CONVERTING FROM INTRINSIC TO EXTRINSIC REPRESENTATIONS :

  procedure Combine_Solutions
              ( file : in file_type;
                n : in natural; sols : in Solution_List;
                b2 : in Vector; w2 : in VecVec );
  procedure Combine_Solutions
              ( file : in file_type;
                n : in natural; sols : in Array_of_Solution_Lists;
                b2 : in Vector; w2 : in Array_of_VecVecs );

  -- DESCRIPTION :
  --   Writes to file the extrinsic format of (multi-homogeneous) solutions.

end P_Intrinsic_Diagonal_Homotopies;
