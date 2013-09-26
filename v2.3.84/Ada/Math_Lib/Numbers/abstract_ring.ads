generic

  type number is private;

  zero,one : number;          -- neutral element for + and * are zero and one

  with function Create ( i : integer ) return number;         -- i is number

  with function Equal ( a,b : number ) return boolean;          -- a = b
  with procedure Copy ( a : in number; b : in out number );     -- b := a

  with function "+" ( a,b : number ) return number;             -- return a+b;
  with function "+" ( a : number )   return number;             -- return +a;
  with function "-" ( a,b : number ) return number;             -- return a-b;
  with function "-" ( a : number )   return number;             -- return -a;
  with function "*" ( a,b : number ) return number;             -- return a*b;

  with procedure Add ( a : in out number; b : in number );      -- a := a+b;
  with procedure Sub ( a : in out number; b : in number );      -- a := a-b;
  with procedure Min ( a : in out number );                     -- a := -a;
  with procedure Mul ( a : in out number; b : in number );      -- a := a*b;

  with procedure Clear ( a : in out number );                   -- deallocation

package Abstract_Ring is end;

-- DESCRIPTION :
--   Abstract specification of any ring of numbers that has the standard
--   integer numbers as subrings.
