with Octo_Double_Ring;
with Octo_Double_Ring.FField;
with Generic_Complex_Numbers;

package OctoDobl_Complex_Numbers is 
  new Generic_Complex_Numbers(Octo_Double_Ring,
                              Octo_Double_Ring.FField);

-- DESCRIPTION :
--   Defines complex arithmetic with octo double numbers.
