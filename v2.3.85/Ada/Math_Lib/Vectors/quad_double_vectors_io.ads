with Quad_Double_Ring_io;
with Quad_Double_Vectors;
with Generic_Vectors_io;

package Quad_Double_Vectors_io is 
  new Generic_Vectors_io(Quad_Double_Ring_io,Quad_Double_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of quad double numbers.
