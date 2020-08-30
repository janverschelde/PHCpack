with Triple_Double_Numbers_io;
with Triple_Double_Ring;
with Abstract_Ring_io;

package Triple_Double_Ring_io is
  new Abstract_Ring_io(Triple_Double_Ring,
                       Triple_Double_Numbers_io.get,
                       Triple_Double_Numbers_io.put,
                       Triple_Double_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for triple double numbers.
