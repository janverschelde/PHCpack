with Quad_Double_Numbers_io;
with Quad_Double_Ring;
with Abstract_Ring_io;

package Quad_Double_Ring_io is
  new Abstract_Ring_io(Quad_Double_Ring,
                       Quad_Double_Numbers_io.get,
                       Quad_Double_Numbers_io.put,
                       Quad_Double_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for quad double numbers.
