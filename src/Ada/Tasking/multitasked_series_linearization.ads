with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_Matrices;
with PentDobl_Complex_VecMats;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_Matrices;
with OctoDobl_Complex_VecMats;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_Matrices;
with DecaDobl_Complex_VecMats;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_VecVecs;
with HexaDobl_Complex_Matrices;
with HexaDobl_Complex_VecMats;

package Multitasked_Series_Linearization is

-- DESCRIPTION :
--   Linearization is applied to solve linear system of truncated series,
--   with multitasking, in double, double double, triple double, quad double,
--   penta double, octo double, deca double, and hexa double precision.

  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return Standard_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return DoblDobl_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return TripDobl_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return QuadDobl_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return PentDobl_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return OctoDobl_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return DecaDobl_Complex_VecVecs.VecVec;
  function Allocate_Work_Space
             ( nbt,dim : integer32 ) return HexaDobl_Complex_VecVecs.VecVec;

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
               A : in TripDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in TripDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( dim : in integer32;
               A : in QuadDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in QuadDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( dim : in integer32;
               A : in PentDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in PentDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( dim : in integer32;
               A : in OctoDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in OctoDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( dim : in integer32;
               A : in DecaDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in DecaDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( dim : in integer32;
               A : in HexaDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in HexaDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the matrix A with the vector x and stores the
  --   result in the vector y, with explicitly declared index variables,
  --   in double, double double, triple double, quad double, penta double,
  --   octo double, deca double, or hexa double precision.

  -- REQUIRED : all vectors and matrices have the same dimension dim.

  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in Standard_Complex_Matrices.Link_to_Matrix;
               x,y : in Standard_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in DoblDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in TripDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in TripDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in QuadDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in QuadDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in PentDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in PentDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in OctoDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in OctoDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in DecaDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in DecaDobl_Complex_Vectors.Link_to_Vector );
  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in HexaDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in HexaDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the matrix A with the vector x and stores the
  --   result in the vector y, with explicitly declared index variables,
  --   in double, double double, triple double, quad double, penta double,
  --   octo double, deca double, or hexa double precision.

  -- REQUIRED : the number of rows of the matrix A is nrows and
  --   the number of columns of the matrix A is ncols, nrows >= ncols;
  --   the vector x has range 1..nbcols and the range of y is 1..nrows.

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in Standard_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in TripDobl_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in QuadDobl_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in PentDobl_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in OctoDobl_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in DecaDobl_Complex_Vectors.Link_to_Vector );
  procedure V_Subtract
              ( dim : in integer32;
                x,y : in HexaDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Subtracts from the vector x the vector y,
  --   with explicitly declared index variables,
  --   in double, double double, triple double, quad double, penta double,
  --   octo double, or deca double precision.

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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking for the backsubstitution (lusolve)
  --   to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

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

  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                qraux : in Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_VecVecs.VecVec;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : in DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_VecVecs.VecVec;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                qraux : in TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_VecVecs.VecVec;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                qraux : in QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_VecVecs.VecVec;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                qraux : in PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_VecVecs.VecVec;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                qraux : in OctoDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out OctoDobl_Complex_VecVecs.VecVec;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                qraux : in DecaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DecaDobl_Complex_VecVecs.VecVec;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                qraux : in HexaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out HexaDobl_Complex_VecVecs.VecVec;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to the least squares solving of the matrix series
  --   equation defined by the matrix series in A and right hand side in b,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

  -- REQUIRED :
  --   A'last = b'last >= 0 and all coefficients in x up to idx-1
  --   have been defined, for idx is in range 1..b'last.

  -- ON ENTRY :
  --   idx      index of the stage, all solutions in b(k),
  --            for k from 0 to idx-1 have been computed.
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series,
  --            A(0) contains the LU factorization of A(0);
  --   b        coefficients of vector series,
  --            for k < idx: b(k) is the k-th solution block,
  --            for k >= idx: b(k) is the (updated) right hand side vector;
  --   x        x(k) for k in 0..idx-1 contains the solutoins;
  --   qraux    information to recover the orthogonal part,
  --            as output of Solve_Lead_by_QRLS;
  --   w1       vector of work space vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows;
  --   w2       vector of work space vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows;
  --   w3       vector of work space vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows;
  --   w4       vector of work space vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows;
  --   w5       vector of work space vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows.
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range at least A(0)'range(1);
  --   output   flag to indicate the extra output is needed.

  -- ON RETURN :
  --   b        all coefficients of the solution series up to b(idx)
  --            and updated right hand side vectors;
  --   x        computed coefficient at idx with respect to input.

  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                S : in Standard_Complex_Vectors.Vector;
                Ut,V : in Standard_Complex_Matrices.Matrix;
                wrk,utb,sub : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                S : in DoblDobl_Complex_Vectors.Vector;
                Ut,V : in DoblDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                S : in TripDobl_Complex_Vectors.Vector;
                Ut,V : in TripDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                S : in QuadDobl_Complex_Vectors.Vector;
                Ut,V : in QuadDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                S : in PentDobl_Complex_Vectors.Vector;
                Ut,V : in PentDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                S : in OctoDobl_Complex_Vectors.Vector;
                Ut,V : in OctoDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                S : in DecaDobl_Complex_Vectors.Vector;
                Ut,V : in DecaDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                S : in HexaDobl_Complex_Vectors.Vector;
                Ut,V : in HexaDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to the least squares solving of the matrix series
  --   equation defined by the matrix series in A and right hand side in b,
  --   using the output of the singular value decomposition,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

  -- REQUIRED :
  --   A'last = b'last >= 0 and all coefficients in x up to idx-1
  --   have been defined, for idx is in range 1..b'last.

  -- ON ENTRY :
  --   idx      index of the stage, all solutions in b(k),
  --            for k from 0 to idx-1 have been computed.
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series,
  --            A(0) contains the LU factorization of A(0);
  --   b        coefficients of vector series,
  --            for k < idx: b(k) is the k-th solution block,
  --            for k >= idx: b(k) is the (updated) right hand side vector;
  --   x        x(k) for k in 0..idx-1 contains the solutoins;
  --   S        vector of singular values;
  --   Ut       the conjugate transpose of the U matrix in the SVD;
  --   V        matrix V in the SVD;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range at least A(0)'range(1);
  --   utb      work space as a vector of vectors of range 1..nbt,
  --            for the product of Ut with b, of range U'range(2);
  --   sub      work space as a vector of vectors of range 1..nbt,
  --            for the product with v, ofr ange v'range(1);
  --   output   flag to indicate the extra output is needed.

  -- ON RETURN :
  --   b        all coefficients of the solution series up to b(idx)
  --            and updated right hand side vectors;
  --   x        computed coefficient at idx with respect to input.

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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   With allocated work space for every task,
  --   repeatedly calls the Multitasked_Solve_Next_by_lusolve
  --   to solve the linear system of power series,
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

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
                w1,w2,w3,w4,w5 : in out Standard_Complex_VecVecs.VecVec;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : in DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_VecVecs.VecVec;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                qraux : in TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_VecVecs.VecVec;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                qraux : in QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_VecVecs.VecVec;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                qraux : in PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_VecVecs.VecVec;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                qraux : in OctoDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out OctoDobl_Complex_VecVecs.VecVec;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                qraux : in DecaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DecaDobl_Complex_VecVecs.VecVec;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                qraux : in HexaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out HexaDobl_Complex_VecVecs.VecVec;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   With allocated work space for every task,
  --   repeatedly calls the Multitasked_Solve_Next_by_QRLS
  --   to solve the linear system of power series,
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

  -- REQUIRED :
  --   A'last = b'last >= 0.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   x        space allocated for the solution series;
  --   qraux    information to recover the orthogonal part;
  --   w1       vector of work space vectors of range 1..nbt,
  --            each work space vector of range 1..n, n = number of rows;
  --   w2       vector of work space vectors of range 1..nbt,
  --            each work space vector of range 1..n, n = number of rows;
  --   w3       vector of work space vectors of range 1..nbt,
  --            each work space vector of range 1..n, n = number of rows;
  --   w4       vector of work space vectors of range 1..nbt,
  --            each work space vector of range 1..n, n = number of rows;
  --   w5       vector of work space vectors of range 1..nbt,
  --            each work space vector of range 1..n, n = number of rows.
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range at least A(0)'range(1);
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   x        all coefficients of the solution series.

  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                S : in Standard_Complex_Vectors.Vector;
                Ut,V : in Standard_Complex_Matrices.Matrix;
                wrk,utb,sub : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                S : in DoblDobl_Complex_Vectors.Vector;
                Ut,V : in DoblDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                S : in TripDobl_Complex_Vectors.Vector;
                Ut,V : in TripDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                S : in QuadDobl_Complex_Vectors.Vector;
                Ut,V : in QuadDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                S : in PentDobl_Complex_Vectors.Vector;
                Ut,V : in PentDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                S : in OctoDobl_Complex_Vectors.Vector;
                Ut,V : in OctoDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                S : in DecaDobl_Complex_Vectors.Vector;
                Ut,V : in DecaDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                S : in HexaDobl_Complex_Vectors.Vector;
                Ut,V : in HexaDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   With allocated work space for every task,
  --   repeatedly calls the Multitasked_Solve_Next_by_SVD
  --   to solve the linear system of power series,
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

  -- REQUIRED :
  --   A'last = b'last >= 0.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   x        lead x(0) of the solution has been computed;
  --   S,V      see the output of Solve_Lead_by_SVD;
  --   Ut       the conjugate transpose of the U matrix of the SVD;
  --   wrk      work space vector for the next coefficient computation;
  --   utb      work space as a vector of vectors of range 1..nbt,
  --            for the product of Ut with b, of range U'range(2);
  --   sub      work space as a vector of vectors of range 1..nbt,
  --            for the product with v, ofr ange v'range(1);
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   b        modified right hand side vectors after back substitution;
  --   x        coefficient vectors of the solution.

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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out triple_double;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out penta_double;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out octo_double;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out deca_double;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out hexa_double;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b,
  --   with an estimate for the inverse of the condition number,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

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
                w1,w2,w3,w4,w5 : in out Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                qraux : out TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                qraux : out PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                qraux : out OctoDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out OctoDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                qraux : out DecaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DecaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32; 
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                qraux : out HexaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out HexaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b,
  --   with QR factorization for least squares solving,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A have the same dimension.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   x        space allocated for the solution series;
  --   w1       work space vector of vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows;
  --   w2       work space vector of vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows;
  --   w3       work space vector of vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows;
  --   w4       work space vector of vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows;
  --   w5       work space vector of vectors of range 1..nbt,
  --            each work space vector has range 1..n, n = number of rows.
  --   wrk      work space as a vector has vectors of range 1..nbt,
  --            with every vector of range at least A(0)'range(1);
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

  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                S : out Standard_Complex_Vectors.Vector;
                U,Ut,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                S : out DoblDobl_Complex_Vectors.Vector;
                U,Ut,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                S : out TripDobl_Complex_Vectors.Vector;
                U,Ut,V : out TripDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out triple_double;
                ewrk : in TripDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                S : out QuadDobl_Complex_Vectors.Vector;
                U,Ut,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                ewrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                S : out PentDobl_Complex_Vectors.Vector;
                U,Ut,V : out PentDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out penta_double;
                ewrk : in PentDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                S : out OctoDobl_Complex_Vectors.Vector;
                U,Ut,V : out OctoDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out octo_double;
                ewrk : in OctoDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                S : out DecaDobl_Complex_Vectors.Vector;
                U,Ut,V : out DecaDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out deca_double;
                ewrk : in DecaDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                S : out HexaDobl_Complex_Vectors.Vector;
                U,Ut,V : out HexaDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out hexa_double;
                ewrk : in HexaDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using the SVD of the
  --   leading coefficient matrix of A for least squares solving,
  --   in double, double double, triple double, quad double, penta
  --   double, octo double, deca double, or hexa double precision,
  --   with multitasking.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A have the same dimension.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series;
  --   x        space allocated for the solution series;
  --   ewrk     work space allocated for the SVD of the lead A(0);
  --   wrkv     work space vectors for the next coefficient computation,
  --            of range 1..nbt, one for every task;
  --   utb      work space as a vector of vectors of range 1..nbt,
  --            for the product of Ut with b, of range U'range(2);
  --   sub      work space as a vector of vectors of range 1..nbt,
  --            for the product with v, ofr ange v'range(1);
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   A        A(0) modified as work space in SVD computation;
  --   b        modified right hand side vectors after back substitution;
  --   S,U,V    see the output of Solve_Lead_by_SVD;
  --   Ut       the conjugate transpose of U, used in the solver;
  --   info     see the output of Solve_Lead_by_SVD;
  --   rcond    inverse condition number computed from the singular
  --            values of the lead coefficient of A;
  --   x        all coefficients of the solution series up to b'last,
  --            provided rcond /= 0.0.

end Multitasked_Series_Linearization;
