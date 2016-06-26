with QuadDobl_Dense_Series_Ring;
with QuadDobl_Dense_Series_Vectors;
with Generic_Matrices;

package QuadDobl_Dense_Series_Matrices is 
  new Generic_Matrices(QuadDobl_Dense_Series_Ring,
                       QuadDobl_Dense_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of quad double dense series.
