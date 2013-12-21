with Standard_Integer_Numbers_io;
with Standard_Integer64_Ring;
with Abstract_Ring_io;

package Standard_Integer64_Ring_io is
  new Abstract_Ring_io(Standard_Integer64_Ring,
                       Standard_Integer_Numbers_io.get,
                       Standard_Integer_Numbers_io.put,
                       Standard_Integer_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of standard integer numbers,
--   of type long long integer.
