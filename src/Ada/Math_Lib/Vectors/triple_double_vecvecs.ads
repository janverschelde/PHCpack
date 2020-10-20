with Generic_VecVecs;
with Triple_Double_Ring;
with Triple_Double_Vectors;

package Triple_Double_VecVecs is 
  new Generic_VecVecs(Triple_Double_Ring,Triple_Double_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of triple double numbers.
