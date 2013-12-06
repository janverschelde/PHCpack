with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Abstract_Ring.Field;
with Generic_Vectors;
with Generic_Laurent_Polynomials;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Field is new Ring.Field(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Polynomials is new Generic_Laurent_Polynomials(Ring);

package Generic_Laur_Poly_Functions is

-- DESCRIPTION :
--   Besides the term by term evaluation, two special data structures are
--   provided for efficient evaluation of polynomials in several variables.
--   With negative exponents, numeric/constraint errors are raised when
--   zero is evaluated.

  use Ring,Field,Vectors,Polynomials;

-- FUNCTION TYPE :

  type Evaluator is access function ( x : Vector ) return number;

-- DATA STRUCTURES :

  type Eval_Poly is private;
  type Eval_Coeff_Poly is private;

-- CONSTRUCTORS :

  function Create ( p : Poly ) return Eval_Poly;
  function Create ( p : Poly ) return Eval_Coeff_Poly;

  procedure Diff ( p : in Poly; i : in integer32;
                   cp : out Eval_Coeff_Poly; m : out Vector );
    -- evaluable coefficient polynomial of the partial derivative,
    -- with m the multiplication factors of the coefficients of p

  function Coeff ( p : Poly ) return Vector;  -- returns coefficient vector

-- EVALUATORS :

  function Eval ( p : Poly; x : number; i : integer32 ) return Poly;
     -- return p(x1,..,xi=x,..,xn);
     -- Number_of_Unknowns(Eval(p,x,i)) = Number_of_Unknowns(p)-1

  function Eval ( d : Degrees; c : number; x : Vector ) return number; 
                                                            -- return c*x**d
  function Eval ( t : Term; c : number; x : Vector ) return number;
                                             -- return c*x**d, with d = t.dg
  function Eval ( t : Term; x : Vector ) return number;

  function Eval ( p : Poly; x : Vector ) return number;       -- return p(x)
  function Eval ( p : Poly; c,x : Vector ) return number;
                     -- return p(c,x), with c = vector of coefficients for p

  function Eval ( p : Eval_Poly; x : Vector ) return number;  -- return p(x)
  function Eval ( p : Eval_Coeff_Poly; c,x : Vector ) return number;
     -- return p(c,x), with c = vector of coefficients for p

-- DESTRUCTORS : deallocate memory.

  procedure Clear ( p : in out Eval_Poly );
  procedure Clear ( p : in out Eval_Coeff_Poly );

private

  type Eval_Poly_Rep;
  type Eval_Coeff_Poly_Rep;

  type Eval_Poly is access Eval_Poly_Rep;
  type Eval_Coeff_Poly is access Eval_Coeff_Poly_Rep;

  type kind is (coefficient,polynomial);
  type Poly_Rec ( k : kind := coefficient ) is record
    case k is
      when coefficient => c : number;
      when polynomial  => p : Eval_Poly;
    end case;
  end record;
  type Coeff_Poly_Rec ( k : kind := coefficient ) is record
    case k is
      when coefficient => c : integer32;
      when polynomial  => p : Eval_Coeff_Poly;
    end case;
  end record;

  type Eval_Poly_Rep is array ( integer32 range <> ) of Poly_Rec;
  type Eval_Coeff_Poly_Rep is array ( integer32 range <> ) of Coeff_Poly_Rec;

end Generic_Laur_Poly_Functions;
