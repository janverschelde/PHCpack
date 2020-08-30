with Penta_Double_Ring_io;
with Penta_Double_Vectors;
with Penta_Double_Matrices;
with Generic_Matrices_io;

package Penta_Double_Matrices_io is 
  new Generic_Matrices_io(Penta_Double_Ring_io,
                          Penta_Double_Vectors,
                          Penta_Double_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of penta double numbers.
