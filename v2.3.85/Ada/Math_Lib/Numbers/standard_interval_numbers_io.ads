with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Interval_Numbers;         use Standard_Interval_Numbers;

package Standard_Interval_Numbers_io is

-- DESCRIPTION :
--   Provides basic input/output routines for intervals.

-- INPUT OPERATIONS :

  procedure get ( i : in out Interval );
  procedure get ( file : file_type; i : in out Interval);

  -- DESCRIPTION :
  --   Reads [a,b] from standard input or from file,
  --   where a and b are floating-point numbers.
  --   Spaces may separate the brackets from the numbers and the comma.

-- OUTPUT OPERATIONS :

  procedure put ( i : in Interval );
  procedure put ( i : in Interval; d : in natural32 );
  procedure put ( file : in file_type; i : in Interval );
  procedure put ( file : in file_type; i : in Interval; d : in natural32 );

  -- DESCRIPTION :
  --   Writes [a,b] to standard output or to file,
  --   using d digits for the floating point numbers.

end Standard_Interval_Numbers_io;
