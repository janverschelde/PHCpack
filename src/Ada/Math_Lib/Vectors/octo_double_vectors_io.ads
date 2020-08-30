with Octo_Double_Ring_io;
with Octo_Double_Vectors;
with Generic_Vectors_io;

package Octo_Double_Vectors_io is 
  new Generic_Vectors_io(Octo_Double_Ring_io,Octo_Double_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of octo double numbers.
