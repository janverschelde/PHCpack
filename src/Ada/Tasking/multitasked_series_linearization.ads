with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
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

  procedure Multitasked_Solve_Next_by_lufac
              ( idx,nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lufac
              ( idx,nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );
  procedure Multitasked_Solve_Next_by_lufac
              ( idx,nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking for the backsubstitution
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
  --   b        all coefficients of the solution series up to b(idx),
  --            provided info = 0, and updated right hand side vectors.

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; output : in boolean := true );
  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to solve the matrix series equation
  --   defined by the matrix series in A and right hand side in b.

  -- REQUIRED :
  --   A'last >= 0 and b'last >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   output   if true, then intermediate output is written,
  --            otherwise, the multitasking remains silent.

  -- ON RETURN :
  --   info     if info /= 0, then the lead coefficient matrix of A
  --            was deemed singular and x is undefined,
  --            if info = 0, then the system is regular;
  --   x        all coefficients of the solution series up to b.deg,
  --            provided info = 0.

end Multitasked_Series_Linearization;
