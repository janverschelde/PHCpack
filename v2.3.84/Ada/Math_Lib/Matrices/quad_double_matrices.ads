with Quad_Double_Ring;                    use Quad_Double_Ring;
with Quad_Double_Vectors;
with Generic_Matrices;

package Quad_Double_Matrices is
  new Generic_Matrices(Quad_Double_Ring,
                       Quad_Double_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of quad double numbers.
