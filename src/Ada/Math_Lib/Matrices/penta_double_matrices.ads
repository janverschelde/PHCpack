with Penta_Double_Ring;                   use Penta_Double_Ring;
with Penta_Double_Vectors;
with Generic_Matrices;

package Penta_Double_Matrices is
  new Generic_Matrices(Penta_Double_Ring,Penta_Double_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of penta double numbers.
