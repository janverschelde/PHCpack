with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Boolean_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;

package Double_Puiseux_Operations is

  procedure Leading_Powers
              ( dim : in integer32; tol : in double_float;
                A : in Standard_Floating_Matrices.Matrix;
                b : in Standard_Floating_Vectors.Vector;
                wrkm : in Standard_Integer_VecVecs.VecVec;
                d : out Standard_Floating_Vectors.Vector;
                idxone,idxtwo : out Standard_Integer_Vectors.Vector;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns in d the leading powers via a tropical Cramer vector
  --   on A and b.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   tol      tolerance to decide whether a number is zero, or equivalently,
  --            whether two numbers are different from each other;
  --   A        matrix of real powers;
  --   b        leading powers of the right hand side vector;
  --   wrkm     pointers to work space for the indices in the
  --            computation of the tropical Cramer vector,
  --            the range of wrkm is 0..dim;
  --   vrblvl   is the verbose level, silent when zero.

  -- ON RETURN :
  --   d        leading powers of the solution,
  --            computed via the tropical Cramer vector.
  --   idxone   indices where the minimum is first obtained;
  --   idxtwo   indices where the minimum is obtained the second time;
  --   fail     true if the minimum is not everywhere exactly obtained twice,
  --            false if the minimum is obtained exactly twice, everywhere.

  procedure Check_Correctness
              ( dim : in integer32; tol : in double_float;
                x,d : in Standard_Floating_Vectors.Vector;
                idx1,idx2 : in Standard_Integer_Vectors.Vector;
                correct : out Boolean_Vectors.Vector;
                cd : in out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Determines which of the entries in d are correct,
  --   using the indices computed by the tropical Cramer vector,
  --   updating the Boolean vector correct.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   tol      tolerance to decide whether a number is zero, or equivalently,
  --            whether two numbers are different from each other;
  --   A        matrix with prescribed locations of the minima;
  --   x        original leading powers of the solution;
  --   d        leading powers of the solution,
  --            computed via the tropical Cramer vector;
  --   idx1     first set of indices in the tropical Cramer vector;
  --   idx2     second set of indices in the tropical Cramer vector;
  --   cd       current correct values;
  --   vrblvl   is the verbose level, silent if zero.

  -- ON RETURN :
  --   correct  updated vector of indices to correct values;
  --   cd       updated correct values.

end Double_Puiseux_Operations;
