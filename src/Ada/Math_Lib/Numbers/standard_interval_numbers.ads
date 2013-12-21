with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;

package Standard_Interval_Numbers is

-- DESCRIPTION :
--   Very basic definition of intervals of standard floating point numbers.

-- EXPORTED DATA TYPE :

  type Interval is private;

-- CREATORS :

  function Create ( i : integer ) return Interval;        -- returns [i,i]
  function Create ( n : natural32 ) return Interval;      -- returns [n,n]
  function Create ( i : integer32 ) return Interval;      -- returns [i,i]
  function Create ( x : double_float ) return Interval;   -- returns [x,x]
  function Create ( a,b : double_float ) return Interval; -- returns [a,b]

-- COPY and EQUAL :

  function Equal ( x,y : Interval ) return boolean;          -- x is y ?
  procedure Copy ( x : in Interval; y : in out Interval );   -- y := x

-- SELECTORS :

  function Left ( i : Interval ) return double_float;   -- returns a
  function Right ( i : Interval ) return double_float;  -- returns b

  function Width ( i : Interval ) return double_float;  -- returns b-a
  function Width ( i : Interval ) return Interval;  -- returns [b-a,b-a]
  function Middle ( i : Interval ) return double_float; -- return (a+b)/2

  function "<" ( x,y : Interval ) return boolean; -- compares Middle's
  function "<" ( x : Interval; y : double_float ) return boolean;
  function "<" ( x : double_float; y : Interval ) return boolean;
  function ">" ( x,y : Interval ) return boolean; -- compares Middle's
  function ">" ( x : Interval; y : double_float ) return boolean;
  function ">" ( x : double_float; y : Interval ) return boolean;

-- ARITHMETICAL OPERATIONS :

  function "+" ( x : Interval ) return Interval;    -- returns +x
  function "+" ( x,y : Interval ) return Interval;  -- returns x+y
  function "+" ( x : Interval; y : double_float ) return Interval;
  function "+" ( x : double_float; y : Interval ) return Interval;
  function "-" ( x : Interval ) return Interval;    -- returns -x
  function "-" ( x,y : Interval ) return Interval;  -- returns x-y
  function "-" ( x : Interval; y : double_float ) return Interval;
  function "-" ( x : double_float; y : Interval ) return Interval;
  function "*" ( x,y : Interval ) return Interval;  -- returns x*y
  function "*" ( x : Interval; y : double_float ) return Interval;
  function "*" ( x : double_float; y : Interval ) return Interval;
  function "/" ( x,y : Interval ) return Interval;  -- returns x/y
  function "/" ( x : Interval; y : double_float ) return Interval;
  function "/" ( x : double_float; y : Interval ) return Interval;

  procedure Add ( x : in out Interval; y : in Interval );   -- x := x + y
  procedure Add ( x : in out Interval; y : in double_float );
  procedure Sub ( x : in out Interval; y : in Interval );   -- x := x - y
  procedure Sub ( x : in out Interval; y : in double_float );
  procedure Min ( x : in out Interval );                    -- x := -x
  procedure Mul ( x : in out Interval; y : in Interval );   -- x := x * y
  procedure Mul ( x : in out Interval; y : in double_float );
  procedure Div ( x : in out Interval; y : in Interval );   -- x := x / y
  procedure Div ( x : in out Interval; y : in double_float );

-- DESTRUCTOR :

  procedure Clear ( x : in out Interval );

private

  type Interval is record
    a,b : double_float;
  end record;

end Standard_Interval_Numbers;
