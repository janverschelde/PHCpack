with text_io;                            use text_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Main_Vertex_Points is

-- DESCRIPTION :
--   This package provides two drivers for extracting the vertex
--   point out of a tuple of point lists.

  procedure Vertex_Points
                ( file : in file_type; L : in out List );
  procedure Vertex_Points 
                ( file : in file_type; L : in out Array_of_Lists );
  procedure Vertex_Points 
                ( file : in file_type; mix : in Link_to_Vector;
                  L : in out Array_of_Lists );

  -- DESCRIPTION :
  --   Reduces the lists to the lists of vertex points.

  -- REQUIRED :
  --   If the type of mixture (mix) is provided, then the tuple of lists
  --   must be sorted according to this vector mix.

  -- ON ENTRY :
  --   file       for writing diagnostics and statistics;
  --   mix        number of different lists in the tuple l,
  --              if not provided, then it will be assumed that all lists
  --              are different from each other;
  --   L          (tuple of) list(s).

  -- ON RETURN :
  --   L          (tuple of) list(s) with nothing except for vertex points.

  procedure Vertex_Points
                ( file : in file_type; p : in out Poly_Sys );
  procedure Vertex_Points
                ( file : in file_type; mix : in Link_to_Vector;
                  p : in out Poly_Sys );

  -- DESCRIPTION :
  --   Reduces the supports of the polynomials to their vertex points.
  --   Merely a driver to the procedures listed above.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a polynomial system, an output file, and then
  --   extracts the vertex points of the support sets of the system.

end Main_Vertex_Points;
