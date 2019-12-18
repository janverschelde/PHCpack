with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Generic_Speelpenning_Convolutions;

package Standard_Speelpenning_Convolutions is
  new Generic_Speelpenning_Convolutions(Standard_Complex_Ring,
                                        Standard_Complex_Vectors,
                                        Standard_Complex_VecVecs);

-- DESCRIPTION :
--   Defines the evaluation and differentiation of products of power series
--   with convolutions computed over the ring of complex numbers,
--   in standard double precision.
