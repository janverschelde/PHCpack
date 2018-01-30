package Boolean_Numbers is

-- DESCRIPTION :
--   With the operations in this package the Boolean_Ring is built.
--   The arithmetic is defined modulo 2.

  function Create ( i : integer ) return boolean;     -- for instantiation

  function Equal ( a,b : boolean ) return boolean;       -- a = b
  procedure Copy ( a : in boolean; b : in out boolean ); -- b := a

  function "+" ( a,b : boolean ) return boolean;         -- return a+b;
  function "+" ( a : boolean )   return boolean;         -- return +a;
  function "-" ( a,b : boolean ) return boolean;         -- return a-b;
  function "-" ( a : boolean )   return boolean;         -- return -a;
  function "*" ( a,b : boolean ) return boolean;         -- return a*b;

  procedure Add ( a : in out boolean; b : in boolean );  -- a := a+b
  procedure Sub ( a : in out boolean; b : in boolean );  -- a := a-b
  procedure Min ( a : in out boolean );                  -- a := -a
  procedure Mul ( a : in out boolean; b : in boolean );  -- a := a*b

  function  Rmd ( a,b : boolean ) return boolean;        -- a mod b
  procedure Rmd ( a : in out boolean; b : in boolean );  -- a := a mod b

  procedure Div ( a : in out boolean; b : in boolean );  -- a := a/b
  procedure Div ( a,b : in boolean;                      -- a = b*q+r
                  q : out boolean; r : out boolean );    -- q := a/b
  procedure Div ( a : in out boolean; b : in boolean;    -- a := a/b
                  r : out boolean );                     -- r := a mod b

  procedure Clear ( a : in out boolean );                -- deallocation

end Boolean_Numbers;
