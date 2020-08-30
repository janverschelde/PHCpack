with Triple_Double_Ring;                  use Triple_Double_Ring;
with Triple_Double_Vectors;
with Generic_Matrices;

package Triple_Double_Matrices is
  new Generic_Matrices(Triple_Double_Ring,Triple_Double_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of triple double numbers.
