with Hexa_Double_Numbers_io;
with Hexa_Double_Ring;
with Abstract_Ring_io;

package Hexa_Double_Ring_io is
  new Abstract_Ring_io(Hexa_Double_Ring,
                       Hexa_Double_Numbers_io.get,
                       Hexa_Double_Numbers_io.put,
                       Hexa_Double_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for hexa double numbers.
