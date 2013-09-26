with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Common_Faces_of_Polytope is

-- DESCRIPTION :
--   This package provides functions to implement the neighborship
--   relations of cells in a mixed subdivision, w.r.t. their faces.

  function Is_Neighbor ( L : List; fc : Face ) return boolean;

  -- DESCRIPTION :
  --   Defines the neighborship relation: returns true
  --   if #intersection(list,fc.points) >= Length_Of(fc.points)-1.

  function Neighboring_Faces
             ( mic : Mixed_Cell; fs : Faces; i : integer32 ) return Faces;

  -- DESCRIPTION :
  --   Returns the neighboring faces of fs to the ith component
  --   of the mixed cell mic.

  function Neighboring_Faces
             ( mic : Mixed_Cell; afs : Array_of_Faces ) return Array_of_Faces;

  -- DESCRIPTION :
  --   Returns the neighboring faces of afs to the mixed cell mic.

end Common_Faces_of_Polytope;
