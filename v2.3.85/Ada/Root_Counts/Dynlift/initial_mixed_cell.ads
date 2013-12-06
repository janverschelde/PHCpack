with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

procedure Initial_Mixed_Cell
                ( n : in integer32; mix : in Vector; pts : in Array_of_Lists;
                  mic : out Mixed_Cell; rest : in out Array_of_Lists );

-- DESCRIPTION :
--   Computes an initial mixed cell for the supports in pts.
--   The lifting for this initial cell equals zero.

-- RECOMMENDED :
--   The list pts consists solely out of vertex points.

-- ON ENTRY :
--   n            the length of the points in pts;
--   mix          type of mixture;
--   pts          the supports, pts'range = mixture'range.
   
-- ON RETURN :
--   mic          an initial mixed cell, with lifting zero,
--                if Mixed_Volume(s) = 0, then Mixed_Volume(pts) = 0;
--   rest         the rest of the supports: pts - s.
