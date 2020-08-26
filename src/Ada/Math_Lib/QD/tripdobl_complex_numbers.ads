with Triple_Double_Ring;
with Triple_Double_Ring.FField;
with Generic_Complex_Numbers;

package TripDobl_Complex_Numbers is 
  new Generic_Complex_Numbers(Triple_Double_Ring,
                              Triple_Double_Ring.FField);

-- DESCRIPTION :
--   Defines complex arithmetic with triple double numbers.
