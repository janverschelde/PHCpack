with Quad_Double_Ring_io;
with Quad_Double_Vectors;
with Quad_Double_Vectors_io;
with Quad_Double_VecVecs;
with Generic_VecVecs_io;

package Quad_Double_VecVecs_io is 
  new Generic_VecVecs_io(Quad_Double_Ring_io,
                         Quad_Double_Vectors,
                         Quad_Double_Vectors_io,
                         Quad_Double_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of quad double numbers.
