with Triple_Double_Ring_io;
with Triple_Double_Vectors;
with Triple_Double_Matrices;
with Generic_Matrices_io;

package Triple_Double_Matrices_io is 
  new Generic_Matrices_io(Triple_Double_Ring_io,
                          Triple_Double_Vectors,
                          Triple_Double_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of triple double numbers.
