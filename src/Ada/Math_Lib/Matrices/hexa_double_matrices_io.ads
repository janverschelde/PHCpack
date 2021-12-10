with Hexa_Double_Ring_io;
with Hexa_Double_Vectors;
with Hexa_Double_Matrices;
with Generic_Matrices_io;

package Hexa_Double_Matrices_io is 
  new Generic_Matrices_io(Hexa_Double_Ring_io,
                          Hexa_Double_Vectors,
                          Hexa_Double_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of hexa double numbers.
