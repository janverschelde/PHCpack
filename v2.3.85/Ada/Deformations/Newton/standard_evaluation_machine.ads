with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;

package Standard_Evaluation_Machine is

-- DESCRIPTION :
--   This package exports the evaluation of polynomials as general
--   functions on vectors of standard complex floating-point numbers.

  procedure Initialize ( p : in Poly );

  -- DESCRIPTION :
  --   Initializes the machine with a polynomial p.

  procedure Initialize ( p : in Poly_Sys );

  -- DESCRIPTION :
  --   Initializes the machine with a polynomial system p.

  procedure Deflate;

  -- DESCRIPTION :
  --   The last value in the value of a polynomial system at a point
  --   is the inverse condition number at that point.

  function Evaluate ( x : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the polynomial p given to Initialize
  --   at the vector x.

  function Evaluate ( x : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the value of the system p given to Initialize 
  --   at the vector x.

  procedure Clear;

  -- DESCRIPTION :
  --   Releases memory occupied by p.

end Standard_Evaluation_Machine;
