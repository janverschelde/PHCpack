with text_io;                            use text_io;
with Abstract_Ring_io;
with Abstract_Ring.Field;
with Generic_Complex_Numbers;

generic

  with package Field is new Abstract_Ring.Field(<>);
  with package Ring_io is new Abstract_Ring_io(Field.Ring);
  with package Complex_Numbers is
         new Generic_Complex_Numbers(Ring_io.Ring,Field);

package Generic_Complex_Numbers_io is

-- DESCRIPTION :
--   This package provides io-routines for complex numbers.
--   A complex number is displayed as two floating numbers, representing
--   respectively the real and imaginary part of the complex number.

  use Complex_Numbers;

  procedure get ( x : in out Complex_Number );
  procedure get ( file : in file_type; x : in out Complex_Number );

  procedure put ( x : in Complex_Number );
  procedure put ( file : in file_type; x : in Complex_Number );

end Generic_Complex_Numbers_io;
