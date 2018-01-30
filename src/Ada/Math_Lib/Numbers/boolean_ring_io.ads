with Boolean_Numbers_io;
with Boolean_Ring;
with Abstract_Ring_io;

package Boolean_Ring_io is
  new Abstract_Ring_io(Boolean_Ring,
                       Boolean_Numbers_io.get,
                       Boolean_Numbers_io.put,
                       Boolean_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for the ring of Boolean numbers.
