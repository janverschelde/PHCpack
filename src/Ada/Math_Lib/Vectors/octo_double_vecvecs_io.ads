with Octo_Double_Ring_io;
with Octo_Double_Vectors;
with Octo_Double_Vectors_io;
with Octo_Double_VecVecs;
with Generic_VecVecs_io;

package Octo_Double_VecVecs_io is 
  new Generic_VecVecs_io(Octo_Double_Ring_io,
                         Octo_Double_Vectors,
                         Octo_Double_Vectors_io,
                         Octo_Double_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of octo double numbers.
