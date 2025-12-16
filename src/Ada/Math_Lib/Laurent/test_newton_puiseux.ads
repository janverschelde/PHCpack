with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
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

  procedure Sum ( cf : in out Standard_Complex_Vectors.Vector;
                  dg : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Exploiting that monomials with the same exponent in dg are consecutive,
  --   adds the coefficients with the same exponent in cf.

  procedure Sum ( cf : in Standard_Complex_VecVecs.VecVec;
                  dg : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Sums the coefficients with same exponents.

  function Positive_Minimum
             ( v : Standard_Floating_Vectors.Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the smallest positive number in v.

  function Positive_Minimum
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the smallest positive number in v,
  --   skipping the entries from which the corresponding c is zero.

  function Positive_Minimum_Index
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns index of the smallest positive number in v,
  --   skipping the entries from which the corresponding c is zero.

  function Coefficient ( c : Standard_Complex_Vectors.Vector;
                         e : Standard_Floating_Vectors.Vector;
                         p : double_float ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the coefficient c(i) for which e(i) = p.

  procedure Leading_Powers_by_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                lcf : in Standard_Complex_Vectors.Vector;
                lpw : in Standard_Floating_Vectors.Vector;
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
  --   lpw      exponents of the second term in the series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      powers computed as positive minima;
  --   cfp      coefficients corresponding to the powers in psm.

  procedure First_Order_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the first order.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of the power series solution;
  --   pwr      exponents of the power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      smallest positive powers in the evaluated series;
  --   csm      coefficients corresponding to the smallest positive powers.

  procedure Second_Order_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes a Taylor series expansion of the Laurent homotopy
  --   using the constant coefficients of the power series solution,
  --   truncated after the second order.
  --   Computes the smallest positive exponents of this evaluation.

  -- ON ENTRY :
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hdg      supports of the Laurent homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of the power series solution;
  --   pwr      exponents of the power series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   psm      smallest positive powers in the evaluated series;
  --   csm      coefficients corresponding to the smallest positive powers.

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

  procedure Run_Newton_Step
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
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
  --   vrblvl   is the verbose level.

  procedure Scale_Homotopy_Powers
              ( hct : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Subtracts the minimum power from all other powers in hct.

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
