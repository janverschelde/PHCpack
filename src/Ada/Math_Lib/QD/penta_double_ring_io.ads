with Penta_Double_Numbers_io;
with Penta_Double_Ring;
with Abstract_Ring_io;

package Penta_Double_Ring_io is
  new Abstract_Ring_io(Penta_Double_Ring,
                       Penta_Double_Numbers_io.get,
                       Penta_Double_Numbers_io.put,
                       Penta_Double_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for penta double numbers.
