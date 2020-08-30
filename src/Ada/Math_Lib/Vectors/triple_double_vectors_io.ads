with Triple_Double_Ring_io;
with Triple_Double_Vectors;
with Generic_Vectors_io;

package Triple_Double_Vectors_io is 
  new Generic_Vectors_io(Triple_Double_Ring_io,Triple_Double_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of triple double numbers.
