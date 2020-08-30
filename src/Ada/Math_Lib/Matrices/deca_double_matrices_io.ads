with Deca_Double_Ring_io;
with Deca_Double_Vectors;
with Deca_Double_Matrices;
with Generic_Matrices_io;

package Deca_Double_Matrices_io is 
  new Generic_Matrices_io(Deca_Double_Ring_io,
                          Deca_Double_Vectors,
                          Deca_Double_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of deca double numbers.
