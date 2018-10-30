with Standard_Dense_Series2_Ring;
with Standard_Dense_Series2_Vectors;
with Generic_Matrices;

package Standard_Dense_Series2_Matrices is 
  new Generic_Matrices(Standard_Dense_Series2_Ring,
                       Standard_Dense_Series2_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard dense series.
