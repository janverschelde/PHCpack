with Standard_Dense_Series_Ring;
with Standard_Dense_Series_Vectors;
with Generic_Matrices;

package Standard_Dense_Series_Matrices is 
  new Generic_Matrices(Standard_Dense_Series_Ring,
                       Standard_Dense_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard dense series.
