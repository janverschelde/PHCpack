with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;
with Generic_VecVecs;
with Generic_Matrices;
with Generic_Monomials;
with Generic_Monomial_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);
  with package Matrices is new Generic_Matrices(Ring,Vectors);
  with package Monomials is new Generic_Monomials(Ring,Vectors,VecVecs);
  with package Polynomials is
    new Generic_Monomial_Vectors(Ring,Vectors,VecVecs,Monomials);

package Generic_Polynomial_Vectors is

-- DESCRIPTION :
--   A polynomial vector is a vector of monomial vectors.

  type Polynomial_Vector is 
    array ( integer32 range <> ) of Polynomials.Link_to_Polynomial;

  type Link_to_Polynomial_Vector is access Polynomial_Vector;

-- EVALUATION and DIFFERENTIATION :

  function Eval ( p : Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector;
  function Eval ( p : Link_to_Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate p at x.

  procedure Diff ( p : in Polynomial_Vector; x : in Vectors.Vector;
                   m : out Matrices.Matrix );
  procedure Diff ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                   m : out Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate all partial
  --   derivatives of p at x.  On return is the a matrix m.

  -- REQUIRED : m'range(1) = p'range and m'range(2) = x'range.

  procedure Speel ( p : in Polynomial_Vector; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector );
  procedure Speel ( p : in Polynomial_Vector; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix );
  procedure Speel ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector );
  procedure Speel ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the Speelpenning algorithm to evaluate all partial
  --   derivatives of p at x, returned in m and p(x) is in y.

  -- ON ENTRY :
  --   p        polynomials in several variables;
  --   x        a vector of values for the variables of p.

  -- ON RETURN :
  --   y        the values of the polynomials in p at x;
  --   m        all partial derivatives of p at x;
  --   yd       vector with the partial derivatives of the current polynomial;
  --   wrk      work space used in the computation of the yd vectors.

-- DESTRUCTORS :

  procedure Clear ( p : in out Polynomial_Vector );
  procedure Clear ( p : in out Link_to_Polynomial_Vector );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial vector m.

end Generic_Polynomial_Vectors;
