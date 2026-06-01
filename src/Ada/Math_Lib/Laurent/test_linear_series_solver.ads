with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Double_Real_Power_Series;
with Double_rpSeries_Vectors;
with Double_rpSeries_Matrices;

package Test_Linear_Series_Solver is

-- DESCRIPTION :
--   Tests solving a linear system of real power series.

  function Random_Series
             ( size : integer32 )
             return Double_Real_Power_Series.Link_to_Series;

  -- DESCRIPTION :
  --   Generates a random series with truncation index equal to size.

  function Random_Series_Vector
             ( dim,size : integer32 )
             return Double_rpSeries_Vectors.Vector;

  -- DESCRIPTION :
  --   Generates a vector of dimension dim of random series
  --   with truncation index equal to size.

  function Random_Series_Matrix
             ( dim,size : integer32 )
             return Double_rpSeries_Matrices.Matrix;

  -- DESCRIPTION :
  --   Generates a matrix of dimension dim of random series
  --   with truncation index equal to size.

  procedure Write ( x : in Double_rpSeries_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes all components of x.

  procedure Write ( A : in Double_rpSeries_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes all elements of A.

  function Right_Hand_Side
             ( A : Double_rpSeries_Matrices.Matrix;
               x : in Double_rpSeries_Vectors.Vector ) 
             return Double_rpSeries_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the product of A with x.

  function Extract_Constants
             ( A : Double_rpSeries_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the constant coefficients of the matrix A.

  function Extract_Leading_Powers
             ( A : Double_rpSeries_Matrices.Matrix )
             return Standard_Floating_Matrices.Matrix;

  -- DESCRIPION :
  --   Returns the matrix of leading powers of A.

  function Extract_Constants 
             ( v : Double_rpSeries_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the constant coefficients of the vector v.

  function Inverse ( A : Standard_Complex_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the inverse of the matrix A,
  --   computed via the singular value decomposition.

  function Matrix_Multiply
             ( A : Standard_Complex_Matrices.Matrix;
               x : Double_rpSeries_Vectors.Vector ) 
             return Double_rpSeries_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of A times x.

  procedure Test_Series_Solver 
              ( A : in Double_rpSeries_Matrices.Matrix;
                x,b : in Double_rpSeries_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Solves the system A*x = b.

  procedure Test ( dim,nbr : in integer32 );

  -- DESCRIPTION :
  --   Generates a dim-by-dim matrix A of series of size nbr,
  --   a dim-dimensional solution vector x of series of nbr size,
  --   and then computes the right hand side vector b.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimensions, generates a random problem,
  --   and then tests the solver.

end Test_Linear_Series_Solver;
