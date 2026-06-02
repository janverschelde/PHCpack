with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Double_rpSeries_Vectors;
with Double_rpSeries_Matrices;

package Double_Linear_rpSeries_Solver is

-- DESCRIPTION :
--   Solves a linear system of real power series.

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

  function Is_In ( A : Standard_Floating_Matrices.Matrix;
                   nbr : double_float; tol : double_float := 1.0E-12 )
                 return boolean;

  -- DESCRIPTION :
  --   Returns true if the number nbr occurs in A,
  --   with respect to the given tolerance.

  procedure Leading_Term
              ( invAb : in Double_rpSeries_Vectors.Vector;
                rA : in Standard_Floating_Matrices.Matrix;
                leadidx : out integer32; 
                leadpow : out double_float; leadcff : out complex_number;
                tol : in double_float := 1.0E-12 );

  -- DESCRIPTION :
  --   Returns in leadidx the index of the row in invAb where the
  --   smallest power not in the matrix rA and returns in leadpow
  --   the value of this smallest power and the corresponding coefficient.

  procedure Next_Term
              ( invAb : in Double_rpSeries_Vectors.Vector;
                rA : in Standard_Floating_Matrices.Matrix;
                powers : in Standard_Floating_Vectors.Vector;
                leadidx : out integer32; 
                leadpow : out double_float; leadcff : out complex_number;
                tol : in double_float := 1.0E-12 );

  -- DESCRIPTION :
  --   Returns in leadidx the index of the row in invAb where the
  --   smallest power not in the matrix rA, not among the already
  --   computed powers, and returns in leadpow the value of this 
  --   smallest power and the corresponding coefficient.

  procedure Real_Power_Series_Solver 
              ( A : in Double_rpSeries_Matrices.Matrix;
                b : in Double_rpSeries_Vectors.Vector;
                z0,z1c : out Standard_Complex_Vectors.Vector;
                z1p : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Solves the system A*x = b, for the special case when A has
  --   only one real power term in every element.

end Double_Linear_rpSeries_Solver;
