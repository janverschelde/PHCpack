with Standard_Integer64_Ring;            use Standard_Integer64_Ring;
with Standard_Integer64_Vectors;
with Generic_Matrices;

package Standard_Integer64_Matrices is
  new Generic_Matrices(Standard_Integer64_Ring,
                       Standard_Integer64_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard integer numbers,
--   of type long long integer.
