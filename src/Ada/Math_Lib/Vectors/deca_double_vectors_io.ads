with Deca_Double_Ring_io;
with Deca_Double_Vectors;
with Generic_Vectors_io;

package Deca_Double_Vectors_io is 
  new Generic_Vectors_io(Deca_Double_Ring_io,Deca_Double_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of deca double numbers.
