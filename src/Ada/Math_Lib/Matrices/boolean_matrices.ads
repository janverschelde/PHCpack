with Boolean_Ring;                       use Boolean_Ring;
with Boolean_Vectors;
with Generic_Matrices;

package Boolean_Matrices is
  new Generic_Matrices(Boolean_Ring,Boolean_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of Boolean numbers.
