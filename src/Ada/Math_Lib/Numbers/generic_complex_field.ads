with Abstract_Ring;
with Abstract_Ring.Field;

generic

  with package The_Real_Ring is new Abstract_Ring(<>);
  with package The_Real_Field is new The_Real_Ring.Field(<>);
  with package The_Complex_Ring is new Abstract_Ring(<>);
  with package The_Complex_Field is new The_Complex_Ring.Field(<>);

package Generic_Complex_Field is 

-- DESCRIPTION :
--   Collects real and complex field, connected with the absolute value.

  package Real_Ring renames The_Real_Ring;
  package Real_Field renames The_Real_Field;
  package Complex_Ring renames The_Complex_Ring;

end Generic_Complex_Field;
