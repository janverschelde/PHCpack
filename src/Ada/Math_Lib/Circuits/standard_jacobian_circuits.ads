with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Standard_Jacobian_Circuits is

-- DESCRIPTION :
--   An arithmetic circuit to evaluate a system of polynomials in several
--   variables and to compute its matrix of all partial derivatives
--   (the Jacobian matrix) consists of a tuple of coefficients and
--   tuples that define the exponents in the supports.

  type Circuit is private;

-- CONSTRUCTORS :

  function Create ( p : Poly_Sys ) return Circuit;

  -- DESCRIPTION :
  --   Returns the data structure that holds the circuit representation
  --   to evaluate and differentiate p.

-- SELECTORS :

  function Number_of_Polynomials ( c : Circuit ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of polynomials in the system.

  function Number_of_Variables ( c : Circuit ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of variables in each polynomial of the system.

  function Number_of_Monomials ( c : Circuit ) return natural32;

  -- DESCRIPTION :
  --   Returns the total number of distinct monomials
  --   as defined by the exponent vectors.

  function Number_of_Terms ( c : Circuit; i : integer32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of terms in the i-th polynomial.

  function Coefficients ( c : Circuit; i : integer32 )
                        return Standard_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of the i-th polynomial.

  function Coefficient ( c : Circuit; i,j : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the j-th coefficient of i-th polynomial.

  function Product ( c : Circuit; i,j : integer32 )
                   return Standard_Natural_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the indices of the variables in the product of variables
  --   for the j-th term in the i-th polynomials.

  function Factor ( c : Circuit; i,j : integer32 )
                  return Standard_Natural_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the common factor for the j-th term in the i-th polynomials.
  --   If there are no common factors, then the null pointer is returned.

-- EVALUATION AND DIFFERENTIATION :

  function WorkSpace ( c : Circuit )
                     return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns work space to evaluate and differentiate the circuit.
  --   The data structure on return should be used in the wrk below.

  procedure EvalDiff ( c : in Circuit;
                       x : in Standard_Complex_Vectors.Vector;
                       wrk : in out Standard_Complex_VecVecs.VecVec;
                       y : out Standard_Complex_Vectors.Vector;
                       A : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the value of the polynomial system defined by a circuit,
  --   and evaluates its Jacobian at x.
  --   This version returns the Jacobian as a regular two dimensional matrix.
 
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
  --   y       y(k) is the value of the k-th polynomial at x;
  --   A       Jacobian matrix of the system evaluated at x.

  procedure EvalDiff ( c : in Circuit;
                       x : in Standard_Complex_Vectors.Vector;
                       wrk : in out Standard_Complex_VecVecs.VecVec;
                       y : out Standard_Complex_Vectors.Vector;
                       A : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the value of the polynomial system defined by a circuit,
  --   and evaluates its Jacobian at x.
  --   This version returns the Jacobian as a vector of columns.
 
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
  --           to be allocated with the WorkSpace function of above;
  --   A       allocated space for a vector of columns.

  -- ON RETURN :
  --   wrk     used workspace, filled with values;
  --   y       y(k) is the value of the k-th polynomial at x;
  --   A       Jacobian matrix of the system evaluated at x.

-- DESTRUCTOR :

  procedure Clear ( c : in out Circuit );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the circuit.

private

  type Circuit_Rep;
  type Circuit is access Circuit_Rep;

end Standard_Jacobian_Circuits;
