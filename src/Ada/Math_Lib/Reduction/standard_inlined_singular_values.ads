with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;

package Standard_Inlined_Singular_Values is

-- DESCRIPTION :
--   The inlined implementation of the singular value decomposition works
--   on the real and imaginary parts of the columns of complex matrices.

-- ACKNOWLEDGMENT :
--   This implementation is based on the Linpack routine zsvd,
--   version dated 03/19/79, correction to shift calculation made 2/85.
--   G.W. Stewart, University of Maryland, Argonne National Lab.

  function Min0 ( a,b : integer32 ) return integer32;

  -- DESCRIPTION : returns the minimum of a and b.

  procedure SVD ( xrv,xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  n,p : in integer32;
                  sr,si : in Standard_Floating_Vectors.Link_to_Vector;
                  er,ei : in Standard_Floating_Vectors.Link_to_Vector;
                  urv,uiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  vrv,viv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  job : in integer32; info : out integer32;
                  rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Takes on input the real and imaginary parts of the columns of
  --   a complex n-by-p matrix and computes its singular values.
  --   As indicated by job, computes the left and right singular vectors.

  -- ON ENTRY :
  --   xrv      real parts of the columns of an n-by-p complex matrix;
  --   xiv      imaginary parts of the columns of an n-by-p complex matrix;
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
  --   sr       space allocated for min(n+1,p) doubles;  
  --   si       space allocated for min(n+1,p) doubles;  
  --   er       space allocated for p doubles;
  --   ei       space allocated for p doubles;
  --   urv      space allocated for n vectors of range 1..n;
  --   uiv      space allocated for n vectors of range 1..n;
  --   vrv      space allocated for p vectors of range 1..p;
  --   viv      space allocated for p vectors of range 1..p;
  --   rwrk     work space for a vector of range 1..n;
  --   iwrk     work space for a vector of range 1..n.

  -- ON RETURN :
  --   xrv      destroyed real parts of the columns;
  --   xiv      destroyed imaginary parts of the columns;
  --   sr       first min(n,p) entries contain the singular values,
  --            arranged in descending order of magnitude;
  --   er       ordinarily a vector of p zeros, see info for exceptions;
  --   urv      real parts of the n left singular vectors;
  --   uiv      imaginary parts of the n left singular vectors;
  --   vrv      real parts of the n right singular vectors;
  --   viv      imaginary parts of the n right singular vectors;
  --   info     if zero, then all singular values are correct,
  --            otherwise only s(info+1), .., s(min(n,p)) are correct.

  procedure SVD ( x : in out Standard_Complex_Matrices.Matrix;
                  n,p : in integer32;
                  s,e : out Standard_Complex_Vectors.Vector;
                  u : out Standard_Complex_Matrices.Matrix;
                  v : out Standard_Complex_Matrices.Matrix;
                  job : in integer32; info : out integer32 );

  -- DESCRIPTION :
  --   Reduces a complex n-by-p matrix x by unitary transformations u
  --   and v to diagonal form.  The diagonal elements in s are the
  --   singular values of x.  The columns of u are the corresponding left
  --   singular vectors, and the columns of v are the right singular vectors.
  --   This procedure has the same parameters as the complex version,
  --   but it wraps the inlined version, intended for testing purposes.

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
  --              b = 1 : return the right singular vectors in v.

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
  --            and b are the same.

end Standard_Inlined_Singular_Values;
