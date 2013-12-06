with text_io;                            use text_io;
with Normal_Cone_Intersections;          use Normal_Cone_Intersections;

package Normal_Cone_Intersections_io is

-- DESCRIPTION :
--   This package provides output routines for intersection matrices.

  procedure put ( ima : in Intersection_Matrix );
  procedure put ( file : in file_type; ima : in Intersection_Matrix );

end Normal_Cone_Intersections_io;
