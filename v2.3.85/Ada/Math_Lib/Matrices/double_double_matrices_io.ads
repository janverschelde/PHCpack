with Double_Double_Ring_io;
with Double_Double_Vectors;
with Double_Double_Matrices;
with Generic_Matrices_io;

package Double_Double_Matrices_io is 
  new Generic_Matrices_io(Double_Double_Ring_io,
                          Double_Double_Vectors,
                          Double_Double_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of double double numbers.
