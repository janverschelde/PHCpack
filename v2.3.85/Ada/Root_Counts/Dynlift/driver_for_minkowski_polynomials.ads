with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

procedure Driver_for_Minkowski_Polynomials
                ( file : in file_type;
                  n : in integer32; mix : in Vector; t : in Triangulation;
                  alltri : in boolean; mixsub : out Mixed_Subdivision );

-- DESCRIPTION :
--   Driver for the computation of the Minkowski-polynomial.

-- ON ENTRY :
--   file         to write all results on;
--   n            dimension before lifting and embedding;
--   mix          type of mixture;
--   t            triangulation of the Cayley polytope;
--   alltri       true when all triangulations are wanted, false otherwise.

-- ON OUTPUT :
--   mixed        mixed subdivision, corresponding the type of mixture.
