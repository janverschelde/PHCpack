with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Octo_Double_Numbers;              use Octo_Double_Numbers;
with OctoDobl_Complex_Vectors;         use OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Matrices;        use OctoDobl_Complex_Matrices;

package OctoDobl_Complex_Singular_Values is

-- DESCRIPTION :
--   This package provides a facility for computing a singular value
--   decomposition of matrices of octo double complex numbers.
--   Besides the SVD, the following operations are provided:
--     rank determination, condition number calculation,
--     pseudo inverse, and solving of a linear system.

-- ACKNOWLEDGMENT :
--   This implementation is a translation of the Linpack routine zsvd,
--   version dated 03/19/79, correction to shift calculation made 2/85.
--   G.W. Stewart, University of Maryland, Argonne National Lab.

  function Min0 ( a,b : integer32 ) return integer32;

  -- DESCRIPTION : returns the minimum of a and b.

  procedure SVD ( x : in out Matrix; n,p : in integer32;
                  s,e : out Vector; u : out Matrix; v : out Matrix;
                  job : in integer32; info : out integer32 );
  procedure SVD ( x : in out Matrix; n,p : in integer32;
                  s,e : out Vector; u : out Matrix; v : out Matrix;
                  job : in integer32; info : out integer32;
                  work : in out Vector );

  -- DESCRIPTION :
  --   Reduces a complex n-by-p matrix x by unitary transformations u
  --   and v to diagonal form.  The diagonal elements in s are the
  --   singular values of x.  The columns of u are the corresponding left
  --   singular vectors, and the columns of v are the right singular vectors.

  -- ON ENTRY :
  --   x        complex n-by-p matrix whose singular value decomposition
  --            is to be computed, x is destroyed by SVD;
  --   n        the number of rows of the matrix x;
  --   p        the number of columns of the matrix x;
  --   job      controls the computation of the singular vectors, it has
  --            the decimal expansion ab with the following meaning:
  --              a = 0 : do not compute the left singular vectors,
  --              a = 1 : return the n left singular vectors in u,
  --              a >=2 : returns the first min(n,p) left singular
  --                      vectors in u,
  --              b = 0 : do not compute the right singular vectors,
  --              b = 1 : return the right singular vectors in v;
  --   work     vector of range 1..n as work space;
  --            this is an optional argument, but without it,
  --            the procedure is not thread safe.

  -- ON RETURN :
  --   s        vector of range 1..mm, where mm = min(n+1,p),
  --            the first min(n,p) entries of s contain the singular values
  --            of x arranged in descending order of magnitude;
  --   e        vector of range 1..p, ordinarily containing zeros,
  --            however see the discussion of info for exceptions;
  --   u        matrix with n rows and k columns, 
  --            if joba = 1, then k = n, if joba >= 2 then k = min(n,p),
  --            u contains the matrix of left singular vectors,
  --            u is not referenced if joba = 0, if n <= p or if joba > 2,
  --            then u may be identified with x in the subroutine call;
  --   v        matrix with p rows and p columns,
  --            v contains the matrix of right singular vectors,
  --            v is not referenced if jobb = 0, if p <= n, then v may be
  --            identified with x in the subroutine call;
  --   info     the singular values (and their corresponding singular vectors)
  --            s(info+1),s(info+2),...,s(m) are correct (here m=min(n,p)),
  --            thus if info = 0, all the singular values and their vectors
  --            are correct, in any event, the matrix b = ctrans(u)*x*v is
  --            the bidiagonal matrix with the elements of s on its diagonal
  --            and the elements of e on its super diagonal (ctrans(u) is the
  --            conjugate-transpose of u), thus the singular values of x 
  --            and b are the same;
  --   work     updated work space vector, if provided on input.

  function Rank ( s : Vector ) return integer32;
  function Rank ( s : Vector; tol : double_float ) return integer32;

  -- DESCRIPTION :
  --   Given the singular values of a matrix, returns the number
  --   of singular values s(i) for which s(i) + 1.0 /= 1.0, or
  --   when tol is present, the number on return counts the singular
  --   values s(i) for which AbsVal(s(i)) >= tol.

  function Inverse_Condition_Number
             ( s : Vector ) return octo_double;
  function Inverse_Condition_Number 
             ( s : Vector; tol : double_float ) return octo_double;

  -- DESCRIPTION :
  --   Given the singular values of a matrix in s, the number on return
  --   is the smallest singular value divided by the lagest one.
  --   Singular values s(i) for which s(i) + 1.0 = 1.0  or for which
  --   AbsVal(s(i)) < tol are ignored.  In this sense the condition
  --   number is taken with respect to the numerical rank (see above)
  --   of the matrix and expresses how good the linear system can be
  --   solved by means of a pseudo inverse.

  function Conjugate_Transpose ( z : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns the conjugated transpose of the matrix z.

  function Inverse ( u,v : Matrix; s : Vector ) return Matrix;
  function Inverse ( u,v : Matrix; s : Vector; tol : double_float )
                   return Matrix;

  -- DESCRIPTION :
  --   Returns the pseudo inverse of the matrix whose singular value
  --   decomposition is given by u, v, and s.  Singular values s(i) for
  --   which AbsVal(s(i)) + 1.0 = 1.0 or AbsVal(s(i)) < tol are ignored.

  function Solve ( u,v : Matrix; s,b : Vector ) return Vector;
  function Solve ( u,v : Matrix; s,b : Vector; tol : double_float ) 
                 return Vector;

  -- DESCRIPTION :
  --   Given the singular value decomposition of the matrix a,
  --   i.e.: u'*a*v = s, the vector x = v*(1/s)*u'*b is returned,
  --   x is the solution of a*x = b.  This solver multiplies b with
  --   the pseudo inverse of the matrix a.
  --   Singular values s(i) for which s(i) + 1.0 = 1.0 are ignored,
  --   or in case tol is provided, singular values s(i) for which
  --   AbsVal(s(i)) < tol are ignored.

  procedure Solve ( ut,v : in Matrix; s,b : in Vector;
                    utb,sub : in out Vector; sol : out Vector );

  -- DESCRIPTION :
  --   Version of Solve with all local variables as vectors passed
  --   as workspace vectors.

  -- ON ENTRY :
  --   ut       the conjugated transpose of the U matrix in the SVD;
  --   v        the V matrix in the SVD;
  --   s        vector of singular values;
  --   b        the righthand side vector of the linear system;
  --   utb      work space for the product ut*b, of range u'range(2);
  --   sub      work space vector of range v(1)'range.

  -- ON RETURN :
  --   utb      equals ut*b;
  --   sub      used as work space;
  --   sol      the solution equals v*sub.

end OctoDobl_Complex_Singular_Values;
