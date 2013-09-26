with Generic_VecMats;
with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;

package Standard_Complex_VecMats is 
  new Generic_VecMats(Standard_Complex_Ring,
                      Standard_Complex_Vectors,
                      Standard_Complex_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of standard complex numbers.
