with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Abstract_Ring;

generic

  with package Ring is new Abstract_Ring(<>);
  with procedure get ( file : in file_type; n : in out Ring.number );
  with procedure put ( file : in file_type; n : in Ring.number );
  with procedure put ( file : in file_type;
                       n : in Ring.number; dp : in natural32 );

package Abstract_Ring_io is end;

-- DESCRIPTION :
--   Abstract specification of the input/output for any ring of numbers.
--   The second put operation allows to specify the number of decimal
--   places that have to be displayed.
