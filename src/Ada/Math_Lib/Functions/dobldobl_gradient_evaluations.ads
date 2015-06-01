with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;

package DoblDobl_Gradient_Evaluations is

-- DESCRIPTION :
--   This package offers routines to evaluate a polynomial in several
--   variables with double double complex coefficients and to compute its
--   gradient, the vector of all partial derivatives.

  procedure Split_Common_Factor
              ( e : in Standard_Natural_Vectors.Vector;
                f,b : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the exponent vector e into a common factor f
  --   and a bit vector b of zeroes and ones.
  --   For all i in e'range: e(i) = f(i) + b(i), b(i) is 0 or 1.

  -- REQUIRED : f'range = b'range = e'range.

  procedure Split_Common_Factors
              ( e : in Standard_Natural_VecVecs.VecVec;
                f,b : out Standard_Natural_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Splits the exponents vectors in e into common factors f
  --   and bit vectors b.

  -- REQUIRED : f'range = b'range = e'range.

  function Reverse_Speel
              ( b : Standard_Natural_VecVecs.VecVec;
                x : DoblDobl_Complex_Vectors.Vector )
              return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Evaluates the monomials defined by the exponent vectors in b
  --   at x, returning the function value and the gradient of the
  --   monomials at x.  On returns, the vectors are of range 0..x'last,
  --   the 0-th position holds the function value and the i-th position
  --   contains the derivative with respect to the i-th variable at x.

  -- REQUIRED : the vectors in b are 0/1 vectors.

  procedure Reverse_Speel
              ( b : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                s : in out DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates the monomials defined by the exponent vectors in b
  --   at x, returning the function value and the gradient of the
  --   monomials at x.  On return in s, the vectors are of range 0..x'last,
  --   the 0-th position holds the function value and the i-th position
  --   contains the derivative with respect to the i-th variable at x.

  -- REQUIRED : the vectors in b are 0/1 vectors and sufficient space
  --   has been allocated for s.

  function Gradient_Monomials
             ( f,b : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Applies the reverse mode to compute the gradient of the Speelpenning
  --   products defined by the vectors in b and multiplies the results with
  --   the factors in f, all evaluated at x.

  procedure Gradient_Monomials
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               x : in DoblDobl_Complex_Vectors.Vector;
               s : in out DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Applies the reverse mode to compute the gradient of the Speelpenning
  --   products defined by the vectors in b and multiplies the results with
  --   the factors in f, all evaluated at x.

  -- REQUIRED : s'range = b'range and s(i)'range = 0..x'last.

  function Gradient_Sum_of_Monomials
             ( f,b : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the evaluation of the sum of the monomials defined by the
  --   Speelpenning products in b and common factors in f, evaluated at x,
  --   along with the gradient vectors.  The range of the vector on return
  --   is 0..x'last with its 0-th component the function value and the i-th
  --   component the i-th derivative of the sum at x.
  --   Note: this evaluation corresponds to a polynomial where all
  --   coefficients are equal to one.

  procedure Gradient_Sum_of_Monomials
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               x : in DoblDobl_Complex_Vectors.Vector;
               y : in out DoblDobl_Complex_VecVecs.VecVec;
               r : out DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in r the evaluation of the sum of the monomials defined by the
  --   Speelpenning products in b and common factors in f, evaluated at x,
  --   along with the gradient vectors.  The range of the vector r must 
  --   be 0..x'last with its 0-th component the function value and the i-th
  --   component the i-th derivative of the sum at x.
  --   Note: this evaluation corresponds to a polynomial where all
  --   coefficients are equal to one.

  -- REQUIRED : The y serves as work space and has been allocated,
  --   in particular y'range = b'range.

  function Gradient_of_Polynomial
             ( f,b : Standard_Natural_VecVecs.VecVec;
               c,x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of the polynomial and its gradient at x,
  --   with exponents in f, b, and coefficients in c.
  --   Compared to the procedures below, this function requires
  --   the allocation and deallocation of workspace memory, which
  --   is not efficient for repeated evaluations and differentiations.

  procedure Gradient_of_Polynomial
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.VecVec;
               ydx : out DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the value of the polynomial and its gradient at x,
  --   with exponents in f, b, and coefficients in c.
  --   This evaluation and differentiation is for general polynomials.
 
  -- REQUIRED :
  --   The range of the vector ydx must be 0..x'last with its 
  --   0-th component the function value and the i-th
  --   component the i-th derivative of the sum at x.
  --   The wrk serves as work space and has been allocated,
  --   in particular wrk'range = b'range.

  -- ON ENTRY :
  --   f       defines the common factors in each monomial,
  --           these are the monomials with higher powers that
  --           appear both in the function values and the derivatives;
  --   b       bit vectors that define products of variables,
  --           b(k)(i) = 1 if the i-th variable occurs in monomial k,
  --           b(k)(i) = 0 otherwise;
  --   c       c(k) is the coefficient of the k-th monomial;
  --   x       values for the variables: where to evaluate at;
  --   wrk     serves as workspace for all evaluated monomials.

  -- ON RETURN :
  --   wrk     used workspace, filled with values;
  --   ydx     ydx(0) is the value of the polynomial at x,
  --           ydx(k) is the k-th derivative of the polynomial at x.

  procedure Conditioned_Gradient_of_Polynomial
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.VecVec;
               ydx : out DoblDobl_Complex_Vectors.Vector;
               numcnd : out Double_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the value of the polynomial and its gradient at x,
  --   with exponents in f, b, and coefficients in c.
  --   This evaluation and differentiation is for general polynomials.
  --   In addition, the numerators of the condition numbers of the
  --   evaluation and differentiation problem are computed as well.
  --   These numerators are the sums of the absolute values of the
  --   evaluated monomials.  The denominators are the absolute values
  --   of the evaluated polynomials.
 
  -- REQUIRED :
  --   The range of the vector ydx must be 0..x'last with its 
  --   0-th component the function value and the i-th
  --   component the i-th derivative of the sum at x.
  --   Moreover: numrco'range = ydx'range.
  --   The y serves as work space and has been allocated,
  --   in particular wrk'range = b'range.

  -- ON ENTRY :
  --   f       defines the common factors in each monomial,
  --           these are the monomials with higher powers that
  --           appear both in the function values and the derivatives;
  --   b       bit vectors that define products of variables,
  --           b(k)(i) = 1 if the i-th variable occurs in monomial k,
  --           b(k)(i) = 0 otherwise;
  --   c       c(k) is the coefficient of the k-th monomial;
  --   x       values for the variables: where to evaluate at;
  --   wrk     serves as workspace for all evaluated monomials.

  -- ON RETURN :
  --   wrk     used workspace, filled with values;
  --   ydx     ydx(0) is the value of the polynomial at x,
  --           ydx(k) is the k-th derivative of the polynomial at x;
  --   numcnd  contains the numerators of the condition numbers of
  --           the polynomial and its derivatives at x.

end DoblDobl_Gradient_Evaluations;
