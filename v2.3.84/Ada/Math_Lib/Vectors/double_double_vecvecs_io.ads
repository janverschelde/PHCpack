with Double_Double_Ring_io;
with Double_Double_Vectors;
with Double_Double_Vectors_io;
with Double_Double_VecVecs;
with Generic_VecVecs_io;

package Double_Double_VecVecs_io is 
  new Generic_VecVecs_io(Double_Double_Ring_io,
                         Double_Double_Vectors,
                         Double_Double_Vectors_io,
                         Double_Double_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of double double numbers.
