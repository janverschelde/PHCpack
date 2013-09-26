with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;

function mv ( fn : in string; nVar,nPts : in integer32;
              ind,cnt,sup : in Standard_Integer_Vectors.Vector )
            return natural32;
--
-- DESCRITPION :
--   Computes the mixed volume of the polytopes spanned by given supports.
--   Writes the resulting mixed-cell configuration to file.
--
-- ON ENTRY :
--   fn      name of the output file for the mixed-cell configuration;
--   nVar    ambient dimension, length of the vectors in supports;
--   nPts    total number of points in the supports;
--   ind     ind(i) is the start of the i-th support;
--   cnt     cnt(i) counts the length of the i-th support;
--   sup     coordinates of the points in the supports.
