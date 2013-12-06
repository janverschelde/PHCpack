with Standard_Complex_Numbers_io;
with Standard_Complex_Ring;
with Abstract_Ring_io;

package Standard_Complex_Ring_io is
  new Abstract_Ring_io(Standard_Complex_Ring,
                       Standard_Complex_Numbers_io.get,
                       Standard_Complex_Numbers_io.put,
                       Standard_Complex_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for standard complex numbers.
