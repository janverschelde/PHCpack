with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;
with Generic_Polynomials;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Polynomials is new Generic_Polynomials(Ring);

package Generic_Polynomial_Functions is

-- DESCRIPTION :
--   Besides the term by term evaluation, two special data structures are
--   provided for efficient evaluation of polynomials in several variables.

  use Ring,Vectors,Polynomials;

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

  function Coeff ( p : Poly ) return Vector;    -- returns coefficient vector

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
  type Poly_Rec is record
    k : kind;
    c : number;
    p : Eval_Poly;
  end record;
  type Coeff_Poly_Rec is record
    k : kind;
    c : integer32;
    p : Eval_Coeff_Poly;
  end record;

-- Note: the original variant records (with k as case) made "valgrind"
--   to report an error in the initialization to coefficient,
--   as c and p occupy the same memory location.

  type Eval_Poly_Rep is array ( integer32 range <> ) of Poly_Rec;
  type Eval_Coeff_Poly_Rep is array ( integer32 range <> ) of Coeff_Poly_Rec;

end Generic_Polynomial_Functions;
