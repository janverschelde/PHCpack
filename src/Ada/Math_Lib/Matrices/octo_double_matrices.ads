with Octo_Double_Ring;                    use Octo_Double_Ring;
with Octo_Double_Vectors;
with Generic_Matrices;

package Octo_Double_Matrices is
  new Generic_Matrices(Octo_Double_Ring,Octo_Double_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of octo double numbers.
