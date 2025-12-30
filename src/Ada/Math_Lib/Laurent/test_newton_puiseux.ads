with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;

package Test_Newton_Puiseux is

-- DESCRIPTION :
--   Tests the Newton-Puiseux algorithm on a random Laurent homotopy.

  procedure Evaluate_and_Differentiate
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                zt0 : in Standard_Complex_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                cjm : out Standard_Complex_Matrices.Matrix;
                ejm : out Standard_Floating_Matrices.Matrix;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates the Laurent homotopy at the leading constants
  --   and computes the Jacobian matrix.
 
  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   zt0      leading, constant coefficients of the series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      leading coefficients of the value;
  --   ydg      leading exponents of the value;
  --   cjm      leading coefficients of the Jacobian matrix;
  --   ejm      leading exponents of the Jacobian matrix.

  procedure Evaluate_All_Monomials
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                zt0 : in Standard_Complex_Vectors.Vector;
                ycf : in Standard_Complex_VecVecs.VecVec;
                ydg : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

   -- DESCRIPTION :
   --   Evaluates all monomials of the Laurent homotopy at a point.

  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   zt0      leading, constant coefficients of the series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ycf      ycf(i) is the coefficient value of monomial i;
  --   ydg      ydg(i) is the exponent value of monomial i.

  procedure Leading_Powers_by_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                lcf : in Standard_Complex_Vectors.Vector;
                psm : out Standard_Floating_Vectors.Vector;
                cfp : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the leading powers of a power series solution
  --   via evaluation of the Laurent homotopy at the constant
  --   of the power series.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   lcf      constant coefficients of a power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      powers computed as positive minima;
  --   cfp      coefficients corresponding to the powers in psm.

  procedure Second_Order_Derivatives
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a Laurent homotopy at a random t value
  --   and then computes the values of all second derivatives
  --   at the power series solution evaluated up to first order.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of the power series solution;
  --   pwr      exponents of the power series solution;
  --   vrblvl   is the verbose level.

  procedure Diagonal_Leading_Terms
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cA : out Standard_Complex_Matrices.Matrix;
                eA : out Standard_Floating_Matrices.Matrix;
                cf1 : out Standard_Complex_Vectors.Vector;
                pw1 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the leading terms of a power series solution,
  --   starting at the constant terms, exploiting the diagonal
  --   structure of the Jacobian matrix.

  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cf0      constant coefficients of a power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cA       coefficients of the Jacobian matrix;
  --   eA       corresponding leading exponents of the Jacobian matrix;
  --   cf1      coefficients corresponding to the exponents in pw1;
  --   pw1      leading exponents of a power series solution.

  procedure Diagonal_Second_Terms
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                tol : in double_float := 1.0E-12;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the second terms of a power series solution,
  --   starting at the constant and leading terms,
  --   exploiting the diagonal structure of the Jacobian matrix.

  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cf0      constant coefficients of a power series solution;
  --   cf1      coefficients corresponding to the exponents in pw1;
  --   pw1      leading exponents of a power series solution;
  --   cA       coefficients of the Jacobian matrix;
  --   tol      tolerance to decide if a number is zero;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cf2      coefficients corresponding to the exponents in pw2;
  --   pw2      second exponents of a power series solution.

  procedure Diagonal_Third_Terms
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cf3 : out Standard_Complex_Vectors.Vector;
                pw3 : out Standard_Floating_Vectors.Vector;
                tol : in double_float := 1.0E-12;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the third terms of a power series solution,
  --   starting at the constant, leading, and second terms,
  --   exploiting the diagonal structure of the Jacobian matrix.

  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cf0      constant coefficients of a power series solution;
  --   cf1      coefficients corresponding to the exponents in pw1;
  --   cf2      coefficients corresponding to the exponents in pw2;
  --   pw1      leading exponents of a power series solution;
  --   pw2      second exponents of a power series solution;
  --   cA       coefficients of the Jacobian matrix;
  --   tol      tolerance to decide if a number is zero;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cf3      coefficients corresponding to the exponents in pw3;
  --   pw3      second exponents of a power series solution.

  procedure Run_Newton_Step
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                tol : in double_float := 1.0E-12;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs one step with Newton's method on a Laurent homotopy,
  --   starting at the leading coefficients of a power series.

  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of a power series solution;
  --   pwr      exponents of a power series solution;
  --   tol      tolerance to decide if a number is zero;
  --   vrblvl   is the verbose level.

  procedure Define_Homotopy
              ( dim : in integer32;
                nbm,nbt : in Standard_Integer_Vectors.Vector;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                hct : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given the dimension, number of monomials and number of terms,
  --   defines a Laurent homotopy with a power series solution.

  -- ON ENTRY :
  --   dim      number of equations and variables in the homotopy;
  --   nbm      nbm(i) equals the number of monomials in the system;
  --   nbt      nbt(i) equals the number of terms in the i-th series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   cff      coefficients of the power series solution;
  --   pwr      real powers of the series solution;
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial.

  procedure Test ( dim : in integer32 );

  -- DESCRIPTION :
  --   Runs a test on a Laurent homotopy of dimension dim,
  --   on a random system with a random series solution.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the the dimension and then launches the test.

end Test_Newton_Puiseux;
