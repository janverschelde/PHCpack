with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Complex64_Numbers;         use Multprec_Complex64_Numbers;

package Multprec_Complex64_Numbers_io is

-- DESCRIPTION :
--   This package provides input/output routines for complex numbers.
--   A complex number is displayed as two floating numbers, representing
--   respectively the real and imaginary part of the complex number.

  procedure get ( x : in out Complex_Number );
  procedure get ( file : in file_type; x : in out Complex_Number );

  procedure put ( x : in Complex_Number );
  procedure put ( file : in file_type; x : in Complex_Number );

  procedure put ( c : in Complex_Number; fore,aft,exp : in natural32 );
  procedure put ( file : in file_type;
                  c : in Complex_Number; fore,aft,exp : in natural32 );

  procedure put ( c : in Complex_Number; dp : in natural32 );
  procedure put ( file : in file_type;
                  c : in Complex_Number; dp : in natural32 ); 

end Multprec_Complex64_Numbers_io;
