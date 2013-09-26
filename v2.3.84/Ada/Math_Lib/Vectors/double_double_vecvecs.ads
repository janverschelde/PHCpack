with Generic_VecVecs;
with Double_Double_Ring;
with Double_Double_Vectors;

package Double_Double_VecVecs is 
  new Generic_VecVecs(Double_Double_Ring,Double_Double_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of double double numbers.
