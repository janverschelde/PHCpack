with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Abstract_Ring.Field;
with Generic_Vectors,Generic_VecVecs;
with Generic_Laurent_Polynomials;
with Generic_Laur_Poly_Functions;
with Generic_Laur_Poly_Systems;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Field is new Ring.Field(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);
  with package Polynomials is new Generic_Laurent_Polynomials(Ring);
  with package Poly_Functions is
         new Generic_Laur_Poly_Functions(Ring,Field,Vectors,Polynomials);
  with package Poly_Systems is
         new Generic_Laur_Poly_Systems(Ring,Polynomials);

package Generic_Laur_System_Functions is 

-- DESCRIPTION :
--   This package provides data structures and evaluation functions for
--   systems of Laurent polynomials.

  use Ring,Field,Vectors,VecVecs,Poly_Functions,Poly_Systems;

-- FUNCTION TYPE :

  type Evaluator is access function ( x : Vector ) return Vector;

-- DATA STRUCTURES :

  type Eval_Laur_Sys is array ( integer32 range <> ) of Eval_Poly;
  type Link_to_Eval_Laur_Sys is access Eval_Laur_Sys;
  type Eval_Coeff_Laur_Sys is array ( integer32 range <> ) of Eval_Coeff_Poly;
  type Link_to_Eval_Coeff_Laur_Sys is access Eval_Coeff_Laur_Sys;

-- CREATORS :

  function Create ( p : Laur_Sys ) return Eval_Laur_Sys;
  function Create ( p : Laur_Sys ) return Eval_Coeff_Laur_Sys;

  function Coeff ( p : Laur_Sys ) return VecVec;

  -- DESCRIPTION :
  --   Returns for each polynomial in p the coefficient vector.
  --   The coefficients on return are used as the "c" parameter in
  --   the evaluation of a coefficient parameter system below.

-- EVALUATORS :

  function Eval ( p : Laur_Sys; x : number; i : integer32 ) return Laur_Sys;
  function Eval ( p : Laur_Sys; x : Vector ) return Vector;
  function Eval ( p : Eval_Laur_Sys; x : Vector ) return Vector;
  function Eval ( p : Eval_Coeff_Laur_Sys; c : VecVec; x : Vector )
                return Vector;

-- DESTRUCTORS :

  procedure Clear ( p : in out Eval_Laur_Sys );
  procedure Clear ( p : in out Link_to_Eval_Laur_Sys );
  procedure Clear ( p : in out Eval_Coeff_Laur_Sys );
  procedure Clear ( p : in out Link_to_Eval_Coeff_Laur_Sys );

end Generic_Laur_System_Functions;
