with Standard_Complex_Series_Ring;
with Standard_Complex_Series_Vectors;
with Generic_Matrices;

package Standard_Complex_Series_Matrices is 
  new Generic_Matrices(Standard_Complex_Series_Ring,
                       Standard_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard complex series.
