with Octo_Double_Numbers_io;
with Octo_Double_Ring;
with Abstract_Ring_io;

package Octo_Double_Ring_io is
  new Abstract_Ring_io(Octo_Double_Ring,
                       Octo_Double_Numbers_io.get,
                       Octo_Double_Numbers_io.put,
                       Octo_Double_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for octo double numbers.
