with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Standard_Floating_Numbers_io is

-- DESCRIPTION :
--   This package defines output of standard floating numbers.

  procedure get ( f : in out single_float );
  procedure get ( file : in file_type; f : in out single_float );
  procedure get ( s : in string; f : in out single_float;
                  last : out positive );
  procedure get ( f : in out double_float );
  procedure get ( file : in file_type; f : in out double_float );
  procedure get ( s : in string; f : in out double_float;
                  last : out positive );

  -- DESCRIPTION :
  --   Reads a floating-point number from standard input, from file,
  --   or from a string.

  procedure put ( f : in single_float );
  procedure put ( s : out string; f : in single_float );
  procedure put ( file : in file_type; f : in single_float );
  procedure put ( f : in double_float );
  procedure put ( s : out string; f : in double_float );
  procedure put ( file : in file_type; f : in double_float );
  procedure put ( f : in single_float; fore,aft,exp : in natural32 );
  procedure put ( f : in double_float; fore,aft,exp : in natural32 );
  procedure put ( s : out string; f : in single_float; aft,exp : in natural32 );
  procedure put ( s : out string; f : in double_float; aft,exp : in natural32 );
  procedure put ( file : in file_type;
                  f : in single_float; fore,aft,exp : in natural32 );
  procedure put ( file : in file_type;
                  f : in double_float; fore,aft,exp : in natural32 );

  -- DESCRIPTION :
  --   Writes a floating-point number to standard output, to a string s,
  --   or to a file which must be opened for output.
  --   The triplet (fore,aft,exp) determines the output format, i.e.:
  --   they indicate the number of decimal places before and after
  --   the decimal point and the number of decimal places in the exponent.

  procedure put ( f : in double_float; dp : in natural32 );
  procedure put ( s : out string; f : in double_float; dp : in natural32 );
  procedure put ( file : in file_type; f : in double_float; dp : in natural32 );

  -- DESCRIPTION : put(f,dp) = put(f,dp,dp,dp).

-- FLOATS in C format :

  procedure printf ( s : in string; f : in double_float );
  pragma Import(C,printf,"printf");

  procedure scanf ( s : in string; f : in out double_float );
  pragma Import(C,scanf,"scanf");

  -- DESCRIPTION :
  --   scanf is similar to the C scanf, s is a control string;
  --   printf writes a float in scientific format, using the C convention.

end Standard_Floating_Numbers_io;
