with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Standard_Complex_Numbers_io is

-- DESCRIPTION :
--   This package provides io-routines for standard complex numbers.

  procedure get ( c : in out Complex_Number );
  procedure get ( file : in file_type; c : in out Complex_Number );

  procedure put ( c : in Complex_Number );
  procedure put ( file : in file_type; c : in Complex_Number );

  procedure put ( c : in Complex_Number; fore,aft,exp : in natural32 );
  procedure put ( file : in file_type;
                  c : in Complex_Number; fore,aft,exp : in natural32 );

  procedure put ( c : in Complex_Number; dp : in natural32 );
  procedure put ( file : in file_type;
                  c : in Complex_Number; dp : in natural32 );

end Standard_Complex_Numbers_io;
