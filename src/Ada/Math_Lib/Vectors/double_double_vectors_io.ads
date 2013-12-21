with Double_Double_Ring_io;
with Double_Double_Vectors;
with Generic_Vectors_io;

package Double_Double_Vectors_io is 
  new Generic_Vectors_io(Double_Double_Ring_io,Double_Double_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of double double numbers.
