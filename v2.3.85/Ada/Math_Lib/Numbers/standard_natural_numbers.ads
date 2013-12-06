package Standard_Natural_Numbers is

-- DESCRIPTION :
--   With the operations in this package the Standard_Natural_Ring is build.

 -- subtype natural64 is long_long_integer range 0..long_long_integer'last;
  type natural32 is new long_integer range 0..long_integer'last;
  type natural64 is new long_long_integer range 0..long_long_integer'last;

  function Create ( i : integer ) return natural32;       -- instantiation
  function Create ( i : integer ) return natural64;       -- uniformity...

  function Equal ( a,b : natural32 ) return boolean;         -- a = b
  function Equal ( a,b : natural64 ) return boolean;
  procedure Copy ( a : in natural32; b : in out natural32 ); -- b := a
  procedure Copy ( a : in natural64; b : in out natural64 );

  procedure Add ( a : in out natural32; b : in natural32 );  -- a := a+b;
  procedure Add ( a : in out natural64; b : in natural64 );
  procedure Sub ( a : in out natural32; b : in natural32 );  -- a := a-b;
  procedure Sub ( a : in out natural64; b : in natural64 );
  procedure Min ( a : in out natural32 );                    -- a := -a;
  procedure Min ( a : in out natural64 );
  procedure Mul ( a : in out natural32; b : in natural32 );  -- a := a*b;
  procedure Mul ( a : in out natural64; b : in natural64 );

  procedure Clear ( a : in out natural32 );                  -- deallocation
  procedure Clear ( a : in out natural64 );                  -- a is set to 0

end Standard_Natural_Numbers;
