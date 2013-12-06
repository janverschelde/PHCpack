with Double_Double_Ring;                  use Double_Double_Ring;
with Double_Double_Vectors;
with Generic_Matrices;

package Double_Double_Matrices is
  new Generic_Matrices(Double_Double_Ring,
                       Double_Double_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of double double numbers.
