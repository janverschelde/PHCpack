with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;

package Multitasked_Series_Linearization is

-- DESCRIPTION :
--   Linearization is applied to solve linear system of truncated series,
--   with multitasking.

  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return Standard_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return DoblDobl_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Allocates work space for nbt tasks, returns a vector of vectors
  --   of range 1..nbt, where each vector has range 1..dim.

  procedure MV_Multiply
             ( dim : in integer32;
               A : in Standard_Complex_Matrices.Link_to_Matrix;
               x,y : in Standard_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( dim : in integer32;
               A : in DoblDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( dim : in integer32;
               A : in QuadDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in QuadDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the matrix A with the vector x and stores the
  --   result in the vector y, with explicitly declared index variables,
  --   in double, double double, or quad double precision.

  -- REQUIRED : all vectors and matrices have the same dimension dim.

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in Standard_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in QuadDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Subtracts from the vector x the vector y,
  --   with explicitly declared index variables,
  --   in double, double double, and quad double precision.

  -- REQUIRED : both x and y have range 1..dim.

  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking for the backsubstitution (lusolve)
  --   to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, and quad double precision.

  -- REQUIRED :
  --   A'last = b'last >= 0.  Moreover, the system is square,
  --   and idx is in range 1..b'last.

  -- ON ENTRY :
  --   idx      index of the stage, all solutions in b(k),
  --            for k from 0 to idx-1 have been computed.
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series,
  --            A(0) contains the LU factorization of A(0);
  --   b        coefficients of vector series,
  --            for k < idx: b(k) is the k-th solution block,
  --            for k >= idx: b(k) is the (updated) right hand side vector;
  --   wrk      allocated work space for the nbt tasks;
  --   output   flag to indicate the extra output is needed.

  -- ON RETURN :
  --   b        all coefficients of the solution series up to b(idx)
  --            and updated right hand side vectors.

  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Allocates work space for every task and
  --   repeatedly calls the Multitasked_Solve_Next_by_lusolve
  --   to solve the linear system of power series,
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, or quad double precision.

  -- REQUIRED :
  --   A'last = b'last >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range 1..dim, dim = A(0)'last(1);
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   b        all coefficients of the solution series.

  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                qraux : in Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : in DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                qraux : in QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Allocates work space for every task and
  --   repeatedly calls the Multitasked_Solve_Next_by_QRLS
  --   to solve the linear system of power series,
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, or quad double precision.

  -- REQUIRED :
  --   A'last = b'last >= 0.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --   x        space allocated for the solution series;
  --   qraux    information to recover the orthogonal part;
  --   w1       work space vector of range 1..n, n = number of rows;
  --   w2       work space vector of range 1..n, n = number of rows;
  --   w3       work space vector of range 1..n, n = number of rows;
  --   w4       work space vector of range 1..n, n = number of rows;
  --   w5       work space vector of range 1..n, n = number of rows.
  --            with every vector of range 1..dim, dim = A(0)'last(1);
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   x        all coefficients of the solution series.

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, or quad double precision.

  -- REQUIRED :
  --   A'last = b'last >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range 1..dim, dim = A(0)'last(1);
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   info     if info /= 0, then the lead coefficient matrix of A
  --            was deemed singular and x is undefined,
  --            if info = 0, then the system is regular;
  --   b        all coefficients of the solution series,
  --            provided info = 0.

  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b,
  --   with an estimate for the inverse of the condition number.

  -- REQUIRED :
  --   A'last = b'last >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range 1..dim, dim = A(0)'last(1);
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   rcond    estimate for the inverse of the condition number,
  --            if rcond + 1.0 = 1.0, then the lead coefficient matrix of A
  --            was deemed singular and x is undefined,
  --            if rcond + 1.0 /= 1.0, then the system is regular;
  --   b        all coefficients of the solution series,
  --            provided rcond + 1.0 /= 1.0.

  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b,
  --   with QR factorization for least squares solving,
  --   in double, double double, or quad double precision.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A have the same dimension.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   x        space allocated for the solution series;
  --   w1       work space vector of range 1..n, n = number of rows;
  --   w2       work space vector of range 1..n, n = number of rows;
  --   w3       work space vector of range 1..n, n = number of rows;
  --   w4       work space vector of range 1..n, n = number of rows;
  --   w5       work space vector of range 1..n, n = number of rows.
  --   wrk      work vector, allocated of range at least A(0)'range(1).
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   A        A(0) contains the output of QRD on A(0);
  --   b        modified right hand side vectors after back substitution;
  --   x        contains the coefficients of the solution series x,
  --            provided info = 0;
  --   qraux    information to recover the orthogonal part;
  --   ipvt     pivoting information if that was requested;
  --   info     is zero of nonsingular, otherwise, a nonzero info
  --            indicates a singular matrix.

end Multitasked_Series_Linearization;
