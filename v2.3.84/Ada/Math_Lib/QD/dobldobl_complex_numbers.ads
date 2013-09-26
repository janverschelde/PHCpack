with Double_Double_Ring;
with Double_Double_Ring.FField;
with Generic_Complex_Numbers;

package DoblDobl_Complex_Numbers is 
  new Generic_Complex_Numbers(Double_Double_Ring,
                              Double_Double_Ring.FField);

-- DESCRIPTION :
--   Defines complex arithmetic with double double numbers.
