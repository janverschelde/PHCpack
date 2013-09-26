with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Triangulations;                     use Triangulations;

package Face_Structures is

-- DESCRIPTION :
--   This package offers a data abstraction to work with face
--   structures of tuples of polytopes.

-- DATA STRUCTURES :

  type Face_Structure is record
    l,last : List;      -- l contains the lifted points of the faces,
                        -- last is a pointer to the last element of l
    t : Triangulation;  -- triangulation of the points in l
    f : Faces;          -- k-faces of conv(l)
  end record;
  type Array_of_Face_Structures is array(integer range <>) of Face_Structure;

-- CONSTRUCTOR :

  procedure Compute_New_Faces
                 ( fs : in out Face_Structure; k,n : in natural;
                   point : in out Link_to_Vector; newfs : out Faces );

  -- DESCRIPTION :
  --   Given a point, the new faces will be computed and returned.

  -- ON ENTRY :
  --   fs          a face structure;
  --   k           dimension of the faces to generate;
  --   n           dimension without the lifting;
  --   point       point that has to be in all the faces.

  -- ON RETURN :
  --   fs          face structure with updated triangulation,
  --               the new faces are not added, also the list is not updated;
  --   point       lifted conservatively w.r.t. fs.t;
  --   newfs       k-faces, which all contain the given point.

-- FLATTENING OPERATIONS :

  procedure Flatten ( l : in out List );
  procedure Flatten ( f : in out Face );
  procedure Flatten ( fs : in out Faces );
  procedure Flatten ( fs : in out Face_Structure );
  procedure Flatten ( fs : in out Array_of_Face_Structures );

  -- DESCRIPTION :
  --   These operations are use to flatten face structures.

-- SELECTORS :

  function Is_In ( fs : Face_Structure; point : vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the point belongs to the face structure.

-- DESTRUCTORS :

  procedure Deep_Clear    ( fs : in out Face_Structure );
  procedure Shallow_Clear ( fs : in out Face_Structure );
  procedure Deep_Clear    ( fs : in out Array_of_Face_Structures );
  procedure Shallow_Clear ( fs : in out Array_of_Face_Structures );

end Face_Structures;
