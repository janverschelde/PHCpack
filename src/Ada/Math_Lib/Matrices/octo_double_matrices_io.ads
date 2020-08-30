with Octo_Double_Ring_io;
with Octo_Double_Vectors;
with Octo_Double_Matrices;
with Generic_Matrices_io;

package Octo_Double_Matrices_io is 
  new Generic_Matrices_io(Octo_Double_Ring_io,
                          Octo_Double_Vectors,
                          Octo_Double_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of octo double numbers.
