with Quad_Double_Ring_io;
with Quad_Double_Vectors;
with Quad_Double_Matrices;
with Generic_Matrices_io;

package Quad_Double_Matrices_io is 
  new Generic_Matrices_io(Quad_Double_Ring_io,
                          Quad_Double_Vectors,
                          Quad_Double_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of quad double numbers.
