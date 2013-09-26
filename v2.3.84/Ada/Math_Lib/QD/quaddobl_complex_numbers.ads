with Quad_Double_Ring;
with Quad_Double_Ring.FField;
with Generic_Complex_Numbers;

package QuadDobl_Complex_Numbers is 
  new Generic_Complex_Numbers(Quad_Double_Ring,
                              Quad_Double_Ring.FField);

-- DESCRIPTION :
--   Defines complex arithmetic with quad double numbers.
