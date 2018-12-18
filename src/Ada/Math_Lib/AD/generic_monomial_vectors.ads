with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;
with Generic_Matrices;
with Generic_Monomials;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Matrices is new Generic_Matrices(Ring,Vectors);
  with package Monomials is new Generic_Monomials(Ring,Vectors,Matrices);

package Generic_Monomial_Vectors is

-- DESCRIPTION :
--   A monomial vector is a vector of monomials in several variables.

  type Monomial_Vector is 
    array ( integer32 range <> ) of Monomials.Link_to_Monomial;

  type Polynomial ( dim,nbr : integer32 ) is record
   -- dim is the ambient dimension, the total number of variables
   -- nbr is the number of monomials, excluding the constant term
    cff0 : Ring.number;             -- constant term of the polynomial
    mons : Monomial_Vector(1..nbr); -- monomials with at least one variable
  end record;

  type Link_to_Monomial_Vector is access Monomial_Vector;
  type Link_to_Polynomial is access Polynomial;

-- EVALUATION and DIFFERENTIATION :

  function Eval ( v : Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number;
  function Eval ( v : Link_to_Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate v at x.

  function Eval ( p : Polynomial;
                  x : Vectors.Vector ) return Ring.number;
  function Eval ( p : Link_to_Polynomial;
                  x : Vectors.Vector ) return Ring.number;

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to evaluate p at x.

  procedure Diff ( v : in Monomial_Vector; x : in Vectors.Vector;
                   yd,wrk : in out Vectors.Vector );
  procedure Diff ( v : in Monomial_Vector; x : in Vectors.Vector;
                   yd : in out Vectors.Vector );
  procedure Diff ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                   yd,wrk : in out Vectors.Vector );
  procedure Diff ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                   yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the straightforward algorithm to compute all derivatives
  --   of the polynomial defined by the sum of the monomials in v at x.

  -- REQUIRED : the ranges for yd and wrk are equal to x'range.

  -- ON ENTRY :
  --   v        a monomial vector;
  --   x        a vector of range 1..m.dim.

  -- ON RETURN :
  --   yd      all partial derivatives of v at x;
  --   wrk     work space to hold intermediate values of the derivatives
  --           of the monomials.

  procedure Speel ( v : in Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel ( v : in Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the monomials
  --   in v and all its derivatives at x.

  -- REQUIRED : m.nvr > 0, at least one variable has exponent 1 and
  --   m.n_base = 0, i.e.: no variables appear with exponent 2 or higher.
  --   The ranges of yd and wrk are equal to x'range.

  -- ON ENTRY :
  --   v       a monomial vector;
  --   x       a vector of range 1..v(i).dim.

  -- ON RETURN :
  --   y       the value of the monomial at x;
  --   yd      all partial derivatives of v at x;
  --   wrk     work space to hold intermediate values of the derivatives
  --           of the monomials.

  procedure Speel ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );
  procedure Speel ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector );
  procedure Speel ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the algorithm of Speelpenning to evaluate the polynomial p
  --   and all its derivatives at x.

  -- REQUIRED : m.nvr > 0, at least one variable has exponent 1 and
  --   m.n_base = 0, i.e.: no variables appear with exponent 2 or higher.
  --   The ranges of yd and wrk are equal to x'range.

  -- ON ENTRY :
  --   p       a polynomial;
  --   x       a vector of range 1..p.dim.

  -- ON RETURN :
  --   y       the value of the monomial at x;
  --   yd      all partial derivatives of p at x;
  --   wrk     work space to hold intermediate values of the derivatives
  --           of the monomials.

-- DESTRUCTORS :

  procedure Clear ( v : in out Monomial_Vector );
  procedure Clear ( v : in out Link_to_Monomial_Vector );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial vector v.

  procedure Clear ( p : in out Polynomial );
  procedure Clear ( p : in out Link_to_Polynomial );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the polynomial p.

end Generic_Monomial_Vectors;
