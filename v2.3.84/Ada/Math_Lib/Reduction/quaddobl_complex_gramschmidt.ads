with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Quad_Double_Numbers;              use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;

package QuadDobl_Complex_GramSchmidt is

-- DESCRIPTION :
--   Offers the QR decomposition with the modified version of
--   the Gram-Schmidt orthonormalization method with 
--   complex quad double arithmetic.
--   The routines do no pivoting.

  procedure QR ( n,m : in integer32;
                 q,r : in out QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Applies the modified Gram-Schmidt orthonormalization method
  --   on the m vectors of dimenion n in v.
 
  -- REQUIRED :
  --   q'range = 1..m, r'range = 1..m and
  --   sufficient memory has been allocated for q and r.
 
  -- ON ENTRY :
  --   n         the dimension of the vectors in v;
  --   m         the number of vectors in v;
  --   q         m vectors of range 1..n are the columns of
  --             an n-by-m matrix A.

  -- ON RETURN :
  --   q         orthonormal basis for the space spanned by
  --             the vectors in v, the columns of q span the
  --             n-by-m matrix Q;
  --   r         the columns of r define an m-by-m upper triangular 
  --             matrix R so that A = Q*R.

  procedure Select_Pivot_Column
               ( k,m : in integer32;
                 v : in QuadDobl_Complex_VecVecs.VecVec;
                 ind : out integer32; max : out quad_double );

  -- DESCRIPTION :
  --   Returns the index i in v(k..m) with the largest max norm.

  -- ON ENTRY :
  --   k         start index in v;
  --   m         index in v to end pivot computation;
  --   v         vector of vectors with range including k..m.
 
  -- ON RETURN :
  --   ind       index of vector in v with largest max norm;
  --   max       the max norm of the vector v(ind).

  procedure Swap ( v : in out QuadDobl_Complex_VecVecs.VecVec;
                   i,j : in integer32 );

  -- DESCRIPTION :
  --   Swaps v(i) with v(j).

  procedure Permute_Upper_Triangular
               ( n : in integer32;
                 r : in out QuadDobl_Complex_VecVecs.VecVec;
                 p : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the swaps done on the original vectors (of range 1..n)
  --   on the upper triangular matrix defined by the vectors in r.

  procedure QR ( n,m : in integer32; tol : in double_float;
                 q,r : in out QuadDobl_Complex_VecVecs.VecVec;
                 pivots : out Standard_Integer_Vectors.Vector;
                 rank : out integer32 );

  -- DESCRIPTION :
  --   Applies the modified Gram-Schmidt orthonormalization method
  --   on the m vectors of dimension n in v.
  --   Does pivoting and returns the rank.
 
  -- REQUIRED :
  --   q'range = 1..m, r'range = 1..m = pivots'range and
  --   sufficient memory has been allocated for q and r.
 
  -- ON ENTRY :
  --   n         the dimension of the vectors in v;
  --   m         the number of vectors in v;
  --   tol       tolerance on the max norm to decide whether
  --             a column in v is zero or not;
  --   q         m vectors of range 1..n are the columns of
  --             an n-by-m matrix A.

  -- ON RETURN :
  --   q         orthonormal basis for the space spanned by
  --             the vectors in v, the columns of q span the
  --             n-by-m matrix Q;
  --   r         the columns of r define an m-by-m upper triangular 
  --             matrix R so that A = Q*R, if pivots(k) /= k,
  --             then the A = Q*R will hold only after the permutation
  --             Permute_Upper_Triangular(n,r,pivots) has been applied;
  --   pivots    pivots(k) indicates which vector in v was picked
  --             in the k-th step, equals zero for k > rank;
  --   rank      number of vectors that are nonzero with respect
  --             to the tolerance.

  procedure Test_Orthonormality
               ( n,m : in integer32;
                 q : in QuadDobl_Complex_VecVecs.VecVec;
                 tol : in double_float; output : in boolean;
                 maxerr : out double_float; fail : out boolean );

  -- DESCRIPTION :
  --   Checks the orthonormality of the m vectors in q.
  --   Every vector in q has range 1..n.

  -- NOTE : if m is larger than n, the orthogonality will fail.

  -- ON ENTRY :
  --   n         dimension of the vectors in q;
  --   m         number of vectors in q;
  --   q         sequence of m vectors of range 1..n;
  --   tol       tolerance on inner products and difference from 1
  --             for the norm of the vectors in q;
  --   output    if true, then all norms and inner products are shown.
  --
  -- ON RETURN :
  --   maxerr    is the largest deviation from the orthornomality,
  --             maxerr > tol implies fail, otherwise, if not fail,
  --             maxerr indicates the accuracy level;
  --   fail      true if the norm of a vector or an inner product
  --             deviates by more than the given tolerance tol,
  --             false otherwise.

  procedure Test_Decomposition
               ( n,m : in integer32;
                 v,q,r : in QuadDobl_Complex_VecVecs.VecVec;
                 tol : in double_float; output : in boolean;
                 maxerr : out double_float; fail : out boolean );

  -- DESCRIPTION :
  --   Checks whether A = Q*R, for A spanned by the vectors in v,
  --   for Q and R defined by the vectors in q and r.
  --   Every vector in q has range 1..n.

  -- ON ENTRY :
  --   n         dimension of the vectors in q;
  --   m         number of vectors in q;
  --   v         original sequence of m vectors of range 1..n;
  --   q         sequence of m vectors of range 1..n,
  --             q is supposed to be orthonormal;
  --   r         sequence of m vectors of range 1..m,
  --             r is supposed to be upper triangular;
  --   tol       tolerance on inner products and difference from 1
  --             for the norm of the vectors in q;
  --   output    if true then the outcome of q*r is shown.

  -- ON RETURN :
  --   maxerr    measures the accuracy of the decomposition,
  --             maxerr > tol implies fail, otherwise, if not fail,
  --             maxerr indicates the accuracy level;
  --   fail      true if the norm of a vector or an inner product
  --             deviates by more than the given tolerance tol,
  --             false otherwise.
 
  function Matrix_Product
              ( n,m : integer32;
                v : QuadDobl_Complex_VecVecs.VecVec;
                x : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns A*x where A is an n-by-m matrix with columns in v.
 
  -- ON ENTRY :
  --   n        the range of the vectors in v is 1..n;
  --   m        there are m vectors in v;
  --   v        m vectors of range 1..n define the columns of A,
  --            so A is an n-by-m matrix;
  --   x        a vector of range 1..m.
  
  -- ON RETURN : y = A*x, where the range of y is 1..n.

  function Matrix_Projection
              ( n,m : integer32;
                q : QuadDobl_Complex_VecVecs.VecVec;
                b : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns Q^T*b where Q is the n-by-m matrix defined by 
  --   the m columns in q.  If the orthogonality of the Q is in doubt,
  --   then the Orthogonal_Projection below is more accurate.

  procedure Orthogonal_Projection
              ( n,m : in integer32;
                q : in QuadDobl_Complex_VecVecs.VecVec;
                b : in out QuadDobl_Complex_Vectors.Vector;
                qb : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   For an orthonormal basis q and some vector b, returns in b
  --   the projection of b onto q.

  -- REQUIRED : the vectors in q are orthonormal
  --   and b is of the same range as the vectors in q.

  -- ON ENTRY :
  --   n        dimension of the vectors in q;
  --   m        number of vectors in q;
  --   q        an orthonormal basis;
  --   b        some vector of range 1..n.

  -- ON RETURN :
  --   b        remainder of the orthogonal projection,
  --            if b lies in the column space of q, then b = 0;
  --   qb       contains in Q^T*b if Q is the matrix with columns in q,
  --            we have qb = Matrix_Project(n,m,q,b), but this procedure
  --            is supposed to be more accurate.

  function Solve ( m : integer32;
                   r : QuadDobl_Complex_VecVecs.VecVec;
                   b : QuadDobl_Complex_Vectors.Vector )
                 return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the solution of the m-by-m upper triangular system R*x = b,
  --   where the matrix R is defined by the columns in r.

  -- REQUIRED : r'range = 1..m = b'range and every vector in r
  --   has range 1..m.

  -- ON ENTRY :
  --   m        dimension of the linear system R*x = b;
  --   r        contains the columns of the upper triangular matrix R;
  --   b        right hand side vector of the linear system.

  -- ON RETURN : R^(-1)*b, a vector of range 1..m.
 
end QuadDobl_Complex_GramSchmidt;
