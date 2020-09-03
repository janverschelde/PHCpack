with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Deca_Double_Vectors;
with Deca_Double_Matrices;

package Deca_Double_QR_Least_Squares is

-- DESCRIPTION :
--   This package provides an implementation of QR-decomposition for
--   matrices in deca double precision.  With this decomposition,
--   we can orthogonalize matrices and solve linear systems in the
--   least squares sense.

  procedure QRD ( x : in out Deca_Double_Matrices.Matrix;
                  qraux : in out Deca_Double_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean );

  -- DESCRIPTION :
  --   Uses Householder transformations to compute the QR decomposition of x.
  --   Column pivoting based on 2-norms of reduced columns is optional.

  -- REQUIRED : jpvt'range = x'range(2) = qraux'range,
  --            x'length(1) >= x'length(2)

  -- ON ENTRY :
  --   x          matrix whose decomposition is to be computed;
  --   jpvt       controls the selection of the pivot columns.
  --              The k-th column x(k) of x is placed in one of three classes
  --              according to the value of jpvt(k):
  --                if jpvt(k) > 0, then x(k) is an initial column.
  --                if jpvt(k) = 0, then x(k) is a free column.
  --                if jpvt(k) < 0, then x(k) is a final column.
  --              before the decomposition is computed, initial columns are
  --              moved to the beginning of the array x and final columns to
  --              the end.  Both initial and final columns are frozen in place
  --              during the computation and only free columns are moved.
  --              At the k-th stage of the reduction, if x(k) is occupied by a
  --              free column it is interchanged with the free column of
  --              largest reduced norm.  jpvt is not referenced if not piv.
  --   piv        column pivoting is performed if and only if piv is true.

  -- ON RETURN :
  --   x          x contains in its upper triangle the upper triangular matrix
  --              R of the QR factorization.  Below its diagonal x contains
  --              information from which the orthogonal part of the
  --              decomposition can be recovered.  Note that if pivoting has
  --              been requested, the decomposition is not that of the original
  --              matrix x but that of x with its columns permuted as described
  --              by jpvt.
  --   qraux      qraux contains further information required to recover
  --              the orthogonal part of the decomposition.
  --   jpvt       jpvt(k) contains the index of the column of the original
  --              matrix that has been interchanged into the k-th column,
  --              if pivoting was requested.
  
  -- ACKNOWLEDGMENT :
  --   This Ada Version is a translation of the LINPACK version, dated 08/14/78,
  --   written by G.W. Stewart, University of Maryland, Argonne National Lab.

  procedure Permute_Columns ( x : in out Deca_Double_Matrices.Matrix;
                              jpvt : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Permutes the columns of the matrix x according to jpvt:
  --     x := (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k))), k = jpvt'last.
  --   This routine is useful to the Least Squares computation.

  procedure Permute ( x : in out Deca_Double_Vectors.Vector;
                      jpvt : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Permutes the columns of the vector x according to jpvt:
  --     x := (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k))), k = jpvt'last.
  --   This routine is useful to the Least Squares computation.

  procedure Basis ( qr : in out Deca_Double_Matrices.Matrix;
                    x : in Deca_Double_Matrices.Matrix );

  -- DESCRIPTION :
  --   Retrieves the orthogonal part of the decomposition.
  --   The columns of qr on output correspond to the column span of x.

  -- IMPORTANT :
  --   Note that this does not work when pivoting was requested.

  -- REQUIRED : qr'range(1) = qr'range(2) = x'range(1).

  -- ON ENTRY :
  --   qr         contains output of the routine qrd;
  --   x          original matrix as part of the input of qrd.

  -- ON RETURN : 
  --   qr         orthogonal part of the QR-decomposition.

  procedure QRLS ( x : in out Deca_Double_Matrices.Matrix;
                   ldx,n,k : in integer32;
                   qraux,y : in Deca_Double_Vectors.Vector;
                   qy,qty,b,rsd,xb : out Deca_Double_Vectors.Vector;
                   job : in integer32; info : out integer32 );

  -- DESCRIPTION :
  --   Applies the output of Deca_Double_QR_Decomposition.QRD to
  --   compute coordinate transformations, projections and least quares
  --   solutions.  For k <= min(n,p), let xk be the matrix
  --
  --             xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
  --
  --   formed from columnns jpvt(1), ... ,jpvt(k) of the original
  --   n x p matrix x that was input to dqrdc (if no pivoting was done,
  --   xk consists of the first k columns of x in their original order).
  --   QRD produces a factored orthogonal matrix q and an upper triangular
  --   matrix r such that
  --                       xk = q * (r)
  --                                (0)
  --   this information is contained in coded form in x and qraux.

  -- ON ENTRY :
  --   x         contains the output of QRD, of size ldx times p;
  --   ldx       leading dimension of x;
  --   n         number of rows in the matrix xk, must be same as in QRD
  --   k         number of columns of the matrix xk, k <= min(n,p)
  --   qraux     contains p entries, auxiliary output from QRD;
  --   y         n-vector to be manipulated by QRLS;
  --   job       specifies what is to be computed, job has the decimal
  --             expansion abcde, with the following meaning :
  --                  if a /= 0, compute qy,
  --                  if b,c,d, or e /= 0, compute qty,
  --                  if c /= 0, compute b,
  --                  if d /= 0, compute rsd,
  --                  if e /= 0, compute xb.
  --             Note that a request to compute b, rsd, or xb
  --             automatically triggers the computation of qty, for
  --             which an array must be provided in the calling sequence.

  -- ON RETURN :
  --   x         may be altered, used as work space;
  --   qy        qy'length = n,
  --             qy conntains q*y, if its computation has been requested;
  --   qty       qty'length = n,
  --             qty contains trans(q)*y, if its computation has been
  --             requested;  here trans(q) is the transpose of the matrix q;
  --   b         b'length = k,
  --             b contains the solution of the least squares problem
  --
  --                    minimize norm2(y - xk*b),
  --
  --             if its computation has been requested.  (Note that
  --             if pivoting was requested in dqrdc, the j-th
  --             component of b will be associated with column jpvt(j)
  --             of the original matrix x that was input into dqrdc.)
  --   rsd       rsd'length = n,
  --             rsd contains the least squares residual y - xk*b,
  --             if its computation has been requested;  rsd is also the
  --             orthogonal projection of y onto the orthogonal complement
  --             of the column space of xk;
  --   xb        x'length = n, xb contains the least squares approximation
  --             xk*b, if its computation has been requested;  xb is also
  --             the orthogonal projection of y onto the column space of x;
  --   info      is zero unless the computation of b has been requested
  --             and r is exactly singular.
  --             In this case, info is the index of the first zero
  --             diagonal element of r and b is left unaltered.

  -- The parameters qy, qty, b, rsd, and xb are not referenced
  -- if their computation is not requested and in this case can be
  -- replaced by dummy variables in the calling program.
  -- To save storage, the user may in some cases use the same array
  -- for different parameters in the calling sequence.
  -- A frequently occuring example is when one wishes to compute
  -- any of b, rsd, or xb and does not need y or qty.  In this
  -- case one may identify y, qty, and one of b, rsd, or xb, while
  -- providing separate arrays for anything else that is to be computed.
  -- Thus the calling sequence
  --
  --     QRLS(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
  --
  -- will result in the computation of b and rsd, with rsd overwriting y.
  -- More generally, each item in the following list contains groups of
  -- permissible identifications for a single callinng sequence.
  --
  --     1. (y,qty,b) (rsd) (xb) (qy)
  --     2. (y,qty,rsd) (b) (xb) (qy)
  --     3. (y,qty,xb) (b) (rsd) (qy)
  --     4. (y,qy) (qty,b) (rsd) (xb)
  --     5. (y,qy) (qty,rsd) (b) (xb)
  --     6. (y,qy) (qty,xb) (b) (rsd)
  --
  -- In any group the value returned in the array allocated to the group
  -- corresponds to the last member of the group.

  -- ACKNOWLEDGMENT :
  --   This Ada version is a translation of the LINPACK version, datet 08/14/78,
  --   written by G.W. Stewart, University of Maryland, Argonne National Lab.

end Deca_Double_QR_Least_Squares;
