with Deca_Double_Numbers_io;
with Deca_Double_Ring;
with Abstract_Ring_io;

package Deca_Double_Ring_io is
  new Abstract_Ring_io(Deca_Double_Ring,
                       Deca_Double_Numbers_io.get,
                       Deca_Double_Numbers_io.put,
                       Deca_Double_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for deca double numbers.
