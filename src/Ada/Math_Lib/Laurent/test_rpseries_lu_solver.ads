with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Double_Real_MatVecs;
with Double_Complex_MatVecs;

package Test_rpSeries_LU_Solver is

-- DESCRIPTION :
--   Tests the solving a linear system of random real powered series
--   via LU factorization.

  procedure Random_rpSeries_Vector
              ( nbt : in integer32;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Generates a random vector of real powered series,
  --   returned in (cff, pwr), with series of size nbt.

  procedure Random_rpSeries_Matrix
              ( nbt : in integer32;
                cff : out Double_Complex_MatVecs.MatVec;
                pwr : out Double_Real_MatVecs.MatVec );

  -- DESCRIPTION :
  --   Generates a random matrix returned by coefficients
  --   and powers in (cff, pwr), with series of size nbt.

  procedure Write ( cff : in Standard_Complex_VecVecs.VecVec;
                    pwr : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Writes the vector with coefficients and powers in (cff, pwr).

  procedure Write ( cff : in Double_Complex_MatVecs.MatVec;
                    pwr : in Double_Real_MatVecs.MatVec );

  -- DESCRIPTION :
  --   Writes the matrix with coefficients and powers in (cff, pwr).

  function Flatten ( v : Standard_Floating_VecVecs.VecVec ) 
                   return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Collects all numbers in v, assuming all vectors have the same size.

  function Index ( v : Standard_Floating_Vectors.Vector;
                   x : double_float ) return integer32;

  -- DESCRIPTION :
  --   Returns the index of x in v, or zero if x is not in v.

  procedure Matrix_Vector_Multiply
              ( Acff : in Double_Complex_MatVecs.MatVec;
                Xcff : in Standard_Complex_VecVecs.VecVec;
                Xpwr : in Standard_Floating_VecVecs.VecVec;
                Bcff : out Standard_Complex_VecVecs.VecVec;
                Bpwr : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Multiplies the matrix in (Acff, Apwr) with the vector in (Xcff, Xpwr)
  --   and returns the result in (Bcff, Bpwr).

  -- ASSUMED : all series have the same size.

  -- REQUIRED : the dimensions must be compatible, that is:
  --    Bcff'range = Acff'range(1), Bpwr'range = Apwr'range(1)
  --    Xcff'range = Acff'range(2), Xpwr'range = Apwr'range(2).

  function Constant_Coefficients
             ( A : Double_Complex_MatVecs.MatVec ) 
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the matrix of constant coefficients of A.

  procedure Solve_Linear_System
              ( Acf0 : in Standard_Complex_Matrices.Matrix;
                Bcff : in Standard_Complex_VecVecs.VecVec;
                Bpwr : in Standard_Floating_VecVecs.VecVec;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Solves the linear system defined by the matrix Acf0
  --   and the right hand side vector in (Bcff, Bpwr),
  --   returning coefficients and powers in (cff, pwr).

  function Difference_Sum
             ( x,y : Standard_Complex_VecVecs.VecVec ) return double_float;

  -- DESCRIPTION :
  --   Returns the sum of the differences in absolute value for all
  --   components of x and y.

  -- REQUIRED : dimensions of x and y must match.

  function Difference_Sum
             ( x,y : Standard_Floating_VecVecs.VecVec ) return double_float;

  -- DESCRIPTION :
  --   Returns the sum of the differences in absolute value for all
  --   components of x and y.

  -- REQUIRED : dimensions of x and y must match.

  procedure Test ( dim,nbt : in integer32 );

  -- DESCRIPTION :
  --   Tests a system of dimension dim with number of terms equal to nbt.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimensions and then generates a test problem.

end Test_rpSeries_LU_Solver;
