with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Standard_Integer_Vectors;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with Double_Double_Matrices;             use Double_Double_Matrices;

package Double_Double_Eigenvalues is

-- DESCRIPTION :
--   This package offers some simple routines to compute all eigenvalues
--   and eigenvectors of a matrix of double double numbers.
--   The code is a translation of the eispack routine rg.

-- TOP LEVEL WRAPPERS :

  procedure Eigenvalues ( A : in Matrix; ierr : out integer32;
                          L : out DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns the eigenvalues of the real matrix A.

  -- REQUIRED : A'range(2) = L'range.

  -- ON ENTRY :
  --   A        a square matrix.

  -- ON RETURN :
  --   ierr     zero when completion was normal, an abnormal termination
  --            occurs when number of iterations exhausts 30*A'last(1);
  --   L        eigenvalues of the matrix A.

  procedure Eigenvectors ( A : in Matrix; ierr : out integer32;
                           L : out DoblDobl_Complex_Vectors.Vector;
                           v : out DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Returns eigenvalues and eigenvectors of a real matrix A.

  -- REQUIRED : A'range(2) = L'range = v'range.

  -- ON ENTRY :
  --   A        a square matrix.

  -- ON RETURN :
  --   ierr     zero when completion was normal, an abnormal termination
  --            occurs when number of iterations exhausts 30*A'last(1);
  --   L        eigenvalues of the matrix A;
  --   v        v(i) is the eigenvector with L(i).

  function Create ( wr,wi : Double_Double_Vectors.Vector )
                  return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Joins the real and imaginary parts of the eigenvalues
  --   into one vector of complex numbers.

  -- REQUIRED : wr'range = wi'range equals range of vector on return.

  function Create ( z : in Double_Double_Matrices.Matrix;
                    wi : in Double_Double_Vectors.Vector ) 
                  return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the eigenvectors stored in the columns of the matrix z.
  --   Because eigenvectors of eigenvalues with nonzero imaginary parts
  --   are stored in consecutive columns of z, the imaginary parts of
  --   the eigenvalues are provided as second argument.

-- MAIN ROUTINE :

  procedure RG ( nm,n : in integer32; A : in Matrix; matz : in integer32;
                 wr,wi : out Double_Double_Vectors.Vector;
                 z : out Matrix; ierr : out integer32 );

  -- DESCRIPTION :
  --   Returns all eigenvalues and eigenvectors of a real general matrix.

  -- ON ENTRY :
  --   nm       number of rows of the matrix A;
  --   n        number of columns of the matrix A;
  --   A        contains the real general matrix;
  --   matz     set to zero if only eigenvalues are desired,
  --            if nonzero, the also eigenvectors will be computed.

  -- ON RETURN :
  --   wr       real parts of the eigenvalues of A;
  --   wi       imaginary parts of the eigenvalues of A,
  --            complex conjugate pairs of eigenvalues appear
  --            consecutively with the eigenvalue having the
  --            positive imaginary part first;
  --   z        contains real and imaginary parts of eigenvectors
  --            (for nonzero value of matz), of same ranges as A:
  --            if the j-th eigenvalue is real, the j-th column of z
  --            contains its eigenvector,
  --            if the j-th eigenvalue is complex with positive imaginary
  --            part, the j-th and (j+1)-th columns of z contain the real
  --            and imaginary parts of its eigenvector, the conjugate of
  --            this vector is the eigenvector for the conjugate eigenvalue;
  --   ierr     is zero for normal completion, otherwise an error
  --            occurred, see the documentation for hqr and hqr2.

-- SUBROUTINES :

  procedure balanc ( nm,n : in integer32; A : in out Matrix;
                     low,igh : out integer32;
                     scale : out Double_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Balances a real matrix and isolated eigenvalues whenever possible.

  -- ON ENTRY :
  --   nm       number of rows of the matrix A;
  --   n        number of columns of the matrix A;
  --   A        contains the real general matrix;

  -- ON RETURN :
  --   A        balanced matrix;
  --   low      integer such that A(i,j) is equal to zero if (1) i > j
  --            and (2) j = 1,2,..,low - 1, or
  --   igh      i = igh+1,..,n;
  --   scale    vector of length n with information determining the
  --            permutations and scaling factors used: suppose that the
  --            principal submatrix in rows low through igh has been 
  --            balanced, that p(j) denotes the index interchanged
  --            with j during the permutation step, and that the elements
  --            of the diagonal matrix used are denoted by d(i,j), then
  --               scale(j) = p(j),    for j = 1,...,low-1
  --                        = d(j,j),      j = low,...,igh
  --                        = p(j)         j = igh+1,...,n.
  --            the order in which the interchanges are made is n to igh+1,
  --            then 1 to low-1.
  --            Note that 1 is returned for igh if igh is zero formally.

  procedure elmhes ( nm,n,low,igh : in integer32; A : in out Matrix;
                     int : out Standard_Integer_Vectors.Vector);

  -- DESCRIPTION :
  --   Reduces a submatrix in rows and columns low through igh to upper 
  --   Hessenberg form by stabilized elementary similarity transformations.

  -- ON ENTRY :
  --   nm       number of rows of the matrix A;
  --   n        number of columns of the matrix A;
  --   low      obtained by balanc, or set to 1 if not used;
  --   igh      obtained by balanc, or set to n if not used;
  --   A        contains the real general matrix.

  -- ON RETURN :
  --   A        contains the hessenberg matrix, the multipliers
  --            which were used in the reduction are stored in the
  --            remaining triangle under the hessenberg matrix;
  --   int      contains information on the rows and columns
  --            interchanged in the reduction, only elements
  --            low through igh are used.

  procedure hqr ( nm,n,low,igh : in integer32; H : in out Matrix;
                  wr,wi : out Double_Double_Vectors.Vector;
                  ierr : out integer32 );

  -- DESCRIPTION :
  --   Finds the eigenvalues of a real upper Hessenberg matrix
  --   by the QR method.

  -- ON ENTRY :
  --   nm       number of rows of the matrix A;
  --   n        number of columns of the matrix A;
  --   low      obtained by balanc, or set to 1 if not used;
  --   igh      obtained by balanc, or set to n if not used;
  --   H        contains the upper Hessenberg matrix, information about
  --            the transformations used in the reduction to Hessenberg
  --            form by elmhes or orthes, if performed, is stored
  --            in the remaining triangle under the Hessenberg matrix.

  -- ON RETURN :
  --   H        has been destroyed, therefore, it must be saved
  --            before calling hqr if subsequent calculation and
  --            back transformation of eigenvectors is to be performed;
  --   wr       real parts of the eigenvalues;
  --   wi       imaginary parts of the eigenvalues,
  --            the eigenvalues are unordered except that complex conjugate 
  --            pairs of values appear consecutively with the eigenvalue
  --            having the positive imaginary part first, if an error exit
  --            is made, the eigenvalues should be correct for indices 
  --            ierr+1,...,n.
  --   ierr     is zero for normal return,
  --            is j if the limit of 30*n iterations is exhausted
  --            while the j-th eigenvalue is being sought.

  procedure eltran ( nm,n,low,igh : in integer32; A : in Matrix;
                     int : in Standard_Integer_Vectors.Vector;
                     z : out Matrix );

  -- DESCRIPTION :
  --   Accumulates the stabilized elementary similarity transformations
  --   used in the reduction of a real general matrix to upper Hessenberg
  --   form by elmhes.

  -- ON ENTRY :
  --   nm       number of rows of the matrix A;
  --   n        number of columns of the matrix A;
  --   low      obtained by balanc, or set to 1 if not used;
  --   igh      obtained by balanc, or set to n if not used;
  --   A        contains the multipliers which were used in the
  --            reduction by elmhes in its lower triangle
  --            below the subdiagonal;
  --   int      contains information on the rows and columns
  --            interchanged in the reduction by elmhes.
  --            only elements low through igh are used.

  -- ON RETURN :
  --   z        contains the transformation matrix produced in the
  --            reduction by elmhes.

  procedure hqr2 ( nm,n,low,igh : in integer32; H : in out Matrix;
                   wr,wi : out Double_Double_Vectors.Vector;
                   z : in out Matrix; ierr : out integer32 );

  -- DESCRIPTION :
  --   Finds the eigenvalues and eigenvectors of a real upper Hessenberg
  --   matrix by the QR method.  The eigenvectors of a real general matrix 
  --   can also be found if elmhes and eltran or orthes and ortran have
  --   been used to reduce this general matrix to Hessenberg form
  --   and to accumulate the similarity transformations.

  -- ON ENTRY :
  --   nm       number of rows of the matrix A;
  --   n        number of columns of the matrix A;
  --   low      obtained by balanc, or set to 1 if not used;
  --   igh      obtained by balanc, or set to n if not used;
  --   H        contains the upper Hessenberg matrix;
  --   z        contains the transformation matrix produced by eltran
  --            after the reduction by elmhes, or by ortran after the
  --            reduction by orthes, if performed; if the eigenvectors
  --            of the Hessenberg matrix are desired, z must contain the
  --            identity matrix;

  -- ON RETURN :
  --   H        has been destroyed;
  --   wr       contains the real parts of the eigenvalues;
  --   wi       contains the imaginary parts of te eigenvalues,
  --            the eigenvalues are unordered except that complex conjugate
  --            pairs of values appear consecutively with the eigenvalue
  --            having the positive imaginary part first;
  --            if an error exit is made, the eigenvalues should be correct
  --            for indices ierr+1,...,n;
  --   z        contains the real and imaginary parts of the eigenvectors,
  --            if the i-th eigenvalue is real, the i-th column of z
  --            contains its eigenvector; if the i-th eigenvalue is complex
  --            with positive imaginary part, the i-th and (i+1)-th
  --            columns of z contain the real and imaginary parts of its
  --            eigenvector; the eigenvectors are unnormalized; if an
  --            error exit is made, none of the eigenvectors has been found;
  --  ierr      is zero for normal return,
  --            is j if the limit of 30*n iterations is exhausted
  --            while the j-th eigenvalue is being sought.

  procedure balbak ( nm,n,low,igh : in integer32;
                     scale : in Double_Double_Vectors.Vector;
                     m : in integer32; z : in out Matrix );

  -- DESCRIPTION :
  --   Forms the eigenvectors of a real general matrix by back transforming
  --   those of the corresponding balanced matrix determined by balanc.

  -- ON ENTRY :
  --   nm       number of rows of the matrix A;
  --   n        number of columns of the matrix A;
  --   low      obtained by balanc;
  --   igh      obtained by balanc;
  --   scale    contains information determining the permutations
  --            and scaling factors used by balanc;
  --   m        is the number of columns of z to be back transformed;
  --   z        contains the real and imaginary parts of the eigenvectors
  --            to be back transformed in its first m columns.

  -- ON RETURN :
  --   z        contains the real and imaginary parts of the
  --            transformed eigenvectors in its first m columns.

-- AUXILIARY ROUTINE for hqr2 :

  procedure cdiv ( ar,ai,br,bi : in double_double;
                   cr,ci : out double_double );

  -- DESCRIPTION :
  --   Returns (cr,ci) = (ar,ai)/(br,bi), the result of the complex
  --   division with numbers given by real and imaginary parts.

end Double_Double_Eigenvalues;
