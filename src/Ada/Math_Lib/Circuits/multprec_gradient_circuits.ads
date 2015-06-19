with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Multprec_Complex_Numbers;            use Multprec_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Polynomials;        use Multprec_Complex_Polynomials;

package Multprec_Gradient_Circuits is

-- DESCRIPTION :
--   An arithmetic circuit to evaluate and differentate a polynomial
--   in several variables is defined as a tuple of coefficients,
--   common factors, and positions of the variables in the products.
--   The coefficients are complex, of arbitrary multiprecision.
--   A circuit to evaluate and differentiate a polynomial system is
--   a tuple of circuits to evaluate and differentiatie polynomials.
--   The purpose of this package is to bundle the different data
--   that defines a circuit into one single record.

  type Circuit is private;

-- CONSTRUCTORS :  

  function Create ( n : natural32;
                    c : Multprec_Complex_Vectors.Vector;
                    b : Standard_Natural_VecVecs.VecVec )
                  return Circuit;

  -- DESCRIPTION :
  --   Returns a circuit defined by the given coefficients in c
  --   and in b the corresponding variables in the products,
  --   to represent a polynomial in n variables.
  --   A deep copy is made of the data in b and c.

  -- REQUIRED : c'range = b'range = 1..m, for m = c'last = b'last.
  --   Moreover: every vector in b must be of range 1..n.

  function Create ( n : natural32;
                    c : Multprec_Complex_Vectors.Vector;
                    b,f : Standard_Natural_VecVecs.VecVec )
                  return Circuit;

  -- DESCRIPTION :
  --   Returns a circuit defined by the given coefficients in c,
  --   in b the corresponding variables in the products, and
  --   in f the common factors of the monomials,
  --   to represent a polynomial in n variables.
  --   A deep copy is made of the data in b, c, and f.

  -- REQUIRED : c'range = b'range = f'range = 1..m, 
  --   for m = c'last = b'last = f'last.
  --   Moreover: m > 0.

  function Create ( p : Poly ) return Circuit;

  -- DESCRIPTION :
  --   Returns the data structure that holds the circuit representation
  --   to evaluate and differentiate the polynomial p.

-- SELECTORS :

  function Number_of_Terms ( c : Circuit ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of coefficients in the circuit,
  --   or equivalently, the number of terms in the polynomial.

  function Number_of_Variables ( c : Circuit ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of variables in the circuit.

  function Coefficients
             ( c : Circuit ) return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of ceofficients of the polynomial.

  function Coefficient
             ( c : Circuit; k : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the k-th coefficient of the polynomial.

  -- REQUIRED : k is in range 1..c.m.

  function Positions
             ( c : Circuit ) return Standard_Natural_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the indices of the variables that occur in the
  --   products of variables in all terms of the polynomial.

  function Positions
             ( c : Circuit; k : integer32 )
             return Standard_Natural_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the indices of the variables that occur in the
  --   product of the k-th term in the polynomial.

  function Factors
             ( c : Circuit ) return Standard_Natural_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns all common factors in all terms of the polynomial.
  --   Note that if there are no common factors, then null is returned.

  function Factors
             ( c : Circuit; k : integer32 )
             return Standard_Natural_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the powers in the common factor of the k-th monomial.

-- EVALUATION AND DIFFERENTIATION :

  function WorkSpace ( c : Circuit )
                     return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns work space to evaluate and differentiate the circuit.
  --   The data structure on return should be used in the wrk below.

  procedure EvalDiff ( c : in Circuit;
                       x : in Multprec_Complex_Vectors.Vector;
                       wrk : in out Multprec_Complex_VecVecs.VecVec;
                       ydx : out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the value of the polynomial defined by a circuit,
  --   and evaluates its gradient at x.
 
  -- REQUIRED :
  --   The range of the vector ydx must be 0..x'last with its 
  --   0-th component the function value and the i-th
  --   component the i-th derivative of the sum at x.
  --   The wrk serves as work space and has been allocated,
  --   in particular wrk'range = 1..Numbers_of_Terms(c).

  -- ON ENTRY :
  --   c       a circuit for a polynomial in several variables;
  --   x       values for the variables: where to evaluate at;
  --   wrk     serves as workspace for all evaluated monomials,
  --           to be allocated with the WorkSpace function of above.

  -- ON RETURN :
  --   wrk     used workspace, filled with values;
  --   ydx     ydx(0) is the value of the polynomial at x,
  --           ydx(k) is the k-th derivative of the polynomial at x.

  function EvalDiff ( c : Circuit; x : Multprec_Complex_Vectors.Vector )
                    return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..x'last that contains at position 0
  --   the value of the polynomial defined by c at x and the other
  --   components of the vector on return contain the value of the
  --   gradient of the polynomial defined by c at x.
  --   The function wraps the procedure EvalDiff and may not be as
  --   efficient for repeated many point evaluations as the work space is
  --   allocated and deallocated with each evaluation and differentiation.

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the circuit.

private

  type Circuit_Rep;
  type Circuit is access Circuit_Rep;

end Multprec_Gradient_Circuits;
