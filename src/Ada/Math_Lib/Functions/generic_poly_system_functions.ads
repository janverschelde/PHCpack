with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors,Generic_VecVecs;
with Generic_Polynomials;
with Generic_Polynomial_Functions;
with Generic_Polynomial_Systems;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);
  with package Polynomials is new Generic_Polynomials(Ring);
  with package Poly_Functions is
         new Generic_Polynomial_Functions(Ring,Vectors,Polynomials);
  with package Poly_Systems is
         new Generic_Polynomial_Systems(Ring,Polynomials);

package Generic_Poly_System_Functions is 

-- DESCRIPTION :
--   This package provides data structures and evaluation functions for
--   systems of polynomials.

  use Ring,Vectors,VecVecs,Poly_Functions,Poly_Systems;

-- FUNCTION TYPE :

  type Evaluator is access function ( x : Vector ) return Vector;

-- DATA STRUCTURES :

  type Eval_Poly_Sys is array ( integer32 range <> ) of Eval_Poly;
  type Link_to_Eval_Poly_Sys is access Eval_Poly_Sys;
  type Eval_Coeff_Poly_Sys is array ( integer32 range <> ) of Eval_Coeff_Poly;
  type Link_to_Eval_Coeff_Poly_Sys is access Eval_Coeff_Poly_Sys;

  type Array_of_Eval_Poly_Sys is
    array ( integer32 range <> ) of Link_to_Eval_Poly_Sys;

-- CREATORS :

  function Create ( p : Poly_Sys ) return Eval_Poly_Sys;
  function Create ( p : Poly_Sys ) return Eval_Coeff_Poly_Sys;

  function Coeff ( p : Poly_Sys ) return VecVec;

  -- DESCRIPTION :
  --   Returns for each polynomial in p the coefficient vector.
  --   The coefficients on return are used as the "c" parameter in
  --   the evaluation of a coefficient parameter system below.

-- EVALUATORS :

  function Eval ( p : Poly_Sys; x : number; i : integer32 ) return Poly_Sys;
  function Eval ( p : Poly_Sys; x : Vector ) return Vector;
  function Eval ( p : Eval_Poly_Sys; x : Vector ) return Vector;
  function Eval ( p : Eval_Coeff_Poly_Sys; c : VecVec; x : Vector )
                return Vector;

-- DESTRUCTORS :

  procedure Clear ( p : in out Eval_Poly_Sys );
  procedure Clear ( p : in out Link_to_Eval_Poly_Sys );
  procedure Clear ( p : in out Eval_Coeff_Poly_Sys );
  procedure Clear ( p : in out Link_to_Eval_Coeff_Poly_Sys );

end Generic_Poly_System_Functions;
