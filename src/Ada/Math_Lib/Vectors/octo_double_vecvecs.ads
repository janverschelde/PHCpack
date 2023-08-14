with Generic_VecVecs;
with Octo_Double_Ring;
with Octo_Double_Vectors;

package Octo_Double_VecVecs is 
  new Generic_VecVecs(Octo_Double_Ring,Octo_Double_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of octo double numbers.
