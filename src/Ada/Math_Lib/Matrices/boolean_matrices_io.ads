with Boolean_Ring_io;
with Boolean_Vectors;
with Boolean_Matrices;
with Generic_Matrices_io;

package Boolean_Matrices_io is 
  new Generic_Matrices_io(Boolean_Ring_io,
                          Boolean_Vectors,
                          Boolean_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of Boolean numbers.
