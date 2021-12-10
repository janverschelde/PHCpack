with Hexa_Double_Ring;                    use Hexa_Double_Ring;
with Hexa_Double_Vectors;
with Generic_Matrices;

package Hexa_Double_Matrices is
  new Generic_Matrices(Hexa_Double_Ring,Hexa_Double_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of hexa double numbers.
