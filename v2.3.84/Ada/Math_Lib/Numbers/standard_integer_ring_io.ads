with Standard_Integer_Numbers_io;
with Standard_Integer_Ring;
with Abstract_Ring_io;

package Standard_Integer_Ring_io is
  new Abstract_Ring_io(Standard_Integer_Ring,
                       Standard_Integer_Numbers_io.get,
                       Standard_Integer_Numbers_io.put,
                       Standard_Integer_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of standard integer numbers.
