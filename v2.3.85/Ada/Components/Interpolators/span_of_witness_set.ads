with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;

package Span_of_Witness_Set is

-- DESCRIPTION :
--   Given a witness set, this package offers procedures to compute the
--   linear span starting at each point in the witness set.

  procedure Standard_Enumerate_Linear_Spans
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all linear spans for all points in the witness set,
  --   sampling the components using standard arithmetic.

  procedure Multprec_Enumerate_Linear_Spans
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all linear spans for all points in the witness set,
  --   sampling the components using multiprecision arithmetic.

end Span_of_Witness_Set;
