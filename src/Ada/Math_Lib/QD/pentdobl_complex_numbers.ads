with Penta_Double_Ring;
with Penta_Double_Ring.FField;
with Generic_Complex_Numbers;

package PentDobl_Complex_Numbers is 
  new Generic_Complex_Numbers(Penta_Double_Ring,
                              Penta_Double_Ring.FField);

-- DESCRIPTION :
--   Defines complex arithmetic with penta double numbers.
