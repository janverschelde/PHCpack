with Boolean_Ring_io;
with Boolean_Vectors;
with Generic_Vectors_io;

package Boolean_Vectors_io is 
  new Generic_Vectors_io(Boolean_Ring_io,Boolean_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors over the ring of Boolean numbers.
