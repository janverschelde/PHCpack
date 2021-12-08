with Hexa_Double_Ring;
with Hexa_Double_Ring.FField;
with Generic_Complex_Numbers;

package HexaDobl_Complex_Numbers is 
  new Generic_Complex_Numbers(Hexa_Double_Ring,
                              Hexa_Double_Ring.FField);

-- DESCRIPTION :
--   Defines complex arithmetic with hexa double numbers.
