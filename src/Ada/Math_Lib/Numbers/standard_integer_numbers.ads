package Standard_Integer_Numbers is

-- DESCRIPTION :
--   With the operations in this package the Standard_Integer_Ring is build.

  type integer32 is new long_integer;
  type integer64 is new long_long_integer; -- up to 2^64 instead of 2^32

  function Create ( i : integer ) return integer32;     -- for instantiation
  function Create ( i : integer ) return integer64;

  function Equal ( a,b : integer32 ) return boolean;         -- a = b
  function Equal ( a,b : integer64 ) return boolean;
  procedure Copy ( a : in integer32; b : in out integer32 ); -- b := a
  procedure Copy ( a : in integer64; b : in out integer64 );

  procedure Add ( a : in out integer32; b : in integer32 );  -- a := a+b
  procedure Add ( a : in out integer64; b : in integer64 );
  procedure Sub ( a : in out integer32; b : in integer32 );  -- a := a-b
  procedure Sub ( a : in out integer64; b : in integer64 );
  procedure Min ( a : in out integer32 );                    -- a := -a
  procedure Min ( a : in out integer64 );
  procedure Mul ( a : in out integer32; b : in integer32 );  -- a := a*b
  procedure Mul ( a : in out integer64; b : in integer64 );

  function  Rmd ( a,b : integer32 ) return integer32;        -- a mod b
  function  Rmd ( a,b : integer64 ) return integer64;
  procedure Rmd ( a : in out integer32; b : in integer32 );  -- a := a mod b
  procedure Rmd ( a : in out integer64; b : in integer64 );

  procedure Div ( a : in out integer32; b : in integer32 );  -- a := a/b
  procedure Div ( a : in out integer64; b : in integer64 );
  procedure Div ( a,b : in integer32;                        -- a = b*q+r
                  q : out integer32; r : out integer32 );    -- q := a/b
  procedure Div ( a,b : in integer64;
                  q : out integer64; r : out integer64 );
  procedure Div ( a : in out integer32; b : in integer32;    -- a := a/b
                  r : out integer32 );                       -- r := a mod b
  procedure Div ( a : in out integer64; b : in integer64;
                  r : out integer64 );

  procedure Clear ( a : in out integer32 );                  -- deallocation
  procedure Clear ( a : in out integer64 );                  -- a is set to 0

end Standard_Integer_Numbers;
