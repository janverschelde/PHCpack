with Penta_Double_Ring_io;
with Penta_Double_Vectors;
with Generic_Vectors_io;

package Penta_Double_Vectors_io is 
  new Generic_Vectors_io(Penta_Double_Ring_io,Penta_Double_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of penta double numbers.
