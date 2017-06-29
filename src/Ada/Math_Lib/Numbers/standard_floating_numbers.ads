--with system;

with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Standard_Floating_Numbers is

-- DESCRIPTION :
--   This package sets floating-point types to be independent
--   of the compiler predefined floating-point declarations.

  type single_float is digits 7;                  -- single precision
  type double_float is digits 15;                 -- double precision
 -- type extra_float is digits system.Max_Digits; -- extra precision
 -- type double_float is digits system.Max_Digits;

  function Create ( i : integer ) return single_float;
  function Create ( i : integer ) return double_float;
  function Create ( i : integer32 ) return single_float;
  function Create ( i : integer32 ) return double_float;
  function Create ( n : natural32 ) return single_float;
  function Create ( n : natural32 ) return double_float;

  function Equal ( a,b : single_float ) return boolean;
  function Equal ( a,b : double_float ) return boolean;

  function AbsVal ( a : single_float ) return single_float;
  function AbsVal ( a : double_float ) return double_float;

  procedure Copy ( a : in single_float; b : in out single_float );
  procedure Copy ( a : in double_float; b : in out double_float );

  procedure Add ( a : in out single_float; b : in single_float ); -- a := a+b
  procedure Add ( a : in out double_float; b : in double_float ); 

  procedure Sub ( a : in out single_float; b : in single_float ); -- a := a-b
  procedure Sub ( a : in out double_float; b : in double_float );

  procedure Min ( a : in out single_float );                      -- a := -a
  procedure Min ( a : in out double_float );

  procedure Mul ( a : in out single_float; b : in single_float ); -- a := a*b
  procedure Mul ( a : in out double_float; b : in double_float );

  procedure Div ( a : in out single_float; b : in single_float ); -- a := a/b
  procedure Div ( a : in out double_float; b : in double_float );

  function Is_Valid ( x : single_float ) return boolean;
  function Is_Valid ( x : double_float ) return boolean;

  -- DESCRIPTION :
  --   A valid number is either positive, negative, or zero,
  --   in which case true is returned, otherwise false is returned.

  procedure Clear ( a : in out single_float );
  procedure Clear ( a : in out double_float );

end Standard_Floating_Numbers;
