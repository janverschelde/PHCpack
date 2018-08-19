with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Standard_Natural_Vectors;
with Generic_Laurent_Polynomials;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Polynomials is new Generic_Laurent_Polynomials(Ring);

package Generic_Laur_Poly_Systems is 

-- DESCRIPTION :
--   This package provides polynomial systems as an array of polynomials
--   with some arithmetic operations.

  use Ring,Polynomials;

-- DATA STRUCTURES :

  type Laur_Sys is array ( integer32 range <> ) of Poly;
  type Link_to_Laur_Sys is access Laur_Sys;
  type Array_of_Laur_Sys is array ( integer32 range <> ) of Link_to_Laur_Sys;
  type Link_to_Array_of_Laur_Sys is access Array_of_Laur_Sys;

-- COPYING :

  procedure Copy ( p : in Laur_Sys; q : in out Laur_Sys );

-- SELECTORS :

  function Variables_in_Support
             ( p : Laur_Sys ) return Standard_Natural_Vectors.Vector;
  -- returns 0/1 vector with 1 at place of a variable with nonzero exponent

  function Size_of_Support ( p : Laur_Sys ) return natural32;
  -- returns the number of variables with nonzero exponent

-- ARITHMETIC OPERATIONS :

  function "+" ( p,q : Laur_Sys ) return Laur_Sys;            -- return p+q
  function "-" ( p,q : Laur_Sys ) return Laur_Sys;            -- return p-q
  function "-" ( p : Laur_Sys ) return Laur_Sys;              -- return -p
  function "*" ( a : number; p : Laur_Sys ) return Laur_Sys;  -- return a*p
  function "*" ( p : Laur_Sys; a : number) return Laur_Sys;   -- return p*a

  procedure Add ( p : in out Laur_Sys; q : in Laur_Sys );     -- p := p+q
  procedure Sub ( p : in out Laur_Sys; q : in Laur_Sys );     -- p := p-q
  procedure Min ( p : in out Laur_Sys );                      -- p := -p
  procedure Mul ( p : in out Laur_Sys; a : in number );       -- p := a*p

-- DIFFERENTIATORS :

  function  Diff ( p : Laur_Sys; i : integer32 ) return Laur_Sys;
  procedure Diff ( p : in out Laur_Sys; i : in integer32 );

-- DESTRUCTORS :

  procedure Clear ( p : in out Laur_Sys );
  procedure Clear ( p : in out Link_to_Laur_Sys );
  procedure Shallow_Clear ( p : in out Link_to_Laur_Sys );
  procedure Clear ( p : in out Array_of_Laur_Sys );
  procedure Clear ( p : in out Link_to_Array_of_Laur_Sys );

  -- DESCRIPTION :
  --   Deallocation of memory, a shallow clear only releases the pointer
  --   to the system; the system is not destroyed.

end Generic_Laur_Poly_Systems;
