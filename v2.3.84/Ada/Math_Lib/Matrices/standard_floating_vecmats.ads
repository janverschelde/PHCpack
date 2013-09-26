with Generic_VecMats;
with Standard_Floating_Ring;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;

package Standard_Floating_VecMats is 
  new Generic_VecMats(Standard_Floating_Ring,
                      Standard_Floating_Vectors,
                      Standard_Floating_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of standard floating numbers.
