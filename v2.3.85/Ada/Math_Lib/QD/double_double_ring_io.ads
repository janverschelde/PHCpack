with Double_Double_Numbers_io;
with Double_Double_Ring;
with Abstract_Ring_io;

package Double_Double_Ring_io is
  new Abstract_Ring_io(Double_Double_Ring,
                       Double_Double_Numbers_io.get,
                       Double_Double_Numbers_io.put,
                       Double_Double_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for double double numbers.
