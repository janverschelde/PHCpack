with Deca_Double_Ring;
with Deca_Double_Ring.FField;
with Generic_Complex_Numbers;

package DecaDobl_Complex_Numbers is 
  new Generic_Complex_Numbers(Deca_Double_Ring,
                              Deca_Double_Ring.FField);

-- DESCRIPTION :
--   Defines complex arithmetic with deca double numbers.
