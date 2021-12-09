with Hexa_Double_Ring_io;
with Hexa_Double_Vectors;
with Generic_Vectors_io;

package Hexa_Double_Vectors_io is 
  new Generic_Vectors_io(Hexa_Double_Ring_io,Hexa_Double_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of hexa double numbers.
