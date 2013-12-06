with Generic_VecVecs;
with Quad_Double_Ring;
with Quad_Double_Vectors;

package Quad_Double_VecVecs is 
  new Generic_VecVecs(Quad_Double_Ring,Quad_Double_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of quad double numbers.
