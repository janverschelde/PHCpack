with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;
with Generic_Matrices;
with Generic_Monomials;
with Generic_Monomial_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Matrices is new Generic_Matrices(Ring,Vectors);
  with package Monomials is new Generic_Monomials(Ring,Vectors,Matrices);
  with package Polynomials is
    new Generic_Monomial_Vectors(Ring,Vectors,Matrices,Monomials);

package Generic_Polynomial_Vectors is

-- DESCRIPTION :
--   A polynomial vector is a vector of monomial vectors.

  type Polynomial_Vector is 
    array ( integer32 range <> ) of Polynomials.Link_to_Polynomial;

  type Link_to_Polynomial_Vector is access Polynomial_Vector;

-- EVALUATORS :

  function Eval ( p : Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector;
  function Eval ( p : Link_to_Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate p at x.

-- DESTRUCTORS :

  procedure Clear ( p : in out Polynomial_Vector );
  procedure Clear ( p : in out Link_to_Polynomial_Vector );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial vector m.

end Generic_Polynomial_Vectors;
