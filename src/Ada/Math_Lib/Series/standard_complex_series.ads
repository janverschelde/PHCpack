with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Generic_Dense_Series;

package Standard_Complex_Series is
  new Generic_Dense_Series(Standard_Complex_Ring,
                           Standard_Complex_Vectors,"/",Div);

-- DESCRIPTION :
--   Defines series using vectors over the ring of standard complex numbers,
--   with the added division operator for the field.
