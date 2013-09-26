with Standard_Natural_Numbers_io;
with Standard_Natural64_Ring;
with Abstract_Ring_io;

package Standard_Natural64_Ring_io is
  new Abstract_Ring_io(Standard_Natural64_Ring,
                       Standard_Natural_Numbers_io.get,
                       Standard_Natural_Numbers_io.put,
                       Standard_Natural_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of standard natural numbers.
