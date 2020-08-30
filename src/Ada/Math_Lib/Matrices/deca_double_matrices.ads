with Deca_Double_Ring;                    use Deca_Double_Ring;
with Deca_Double_Vectors;
with Generic_Matrices;

package Deca_Double_Matrices is
  new Generic_Matrices(Deca_Double_Ring,Deca_Double_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of deca double numbers.
