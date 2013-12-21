with text_io;                           use text_io;
with Integer_Faces_of_Polytope;         use Integer_Faces_of_Polytope;

package Integer_Faces_of_Polytope_io is

-- DESCRIPTION :
--   Output of faces of integer polytopes.

  procedure put ( f : in Face );
  procedure put ( file : in file_type; f : in Face );

  procedure put ( f : in Faces );
  procedure put ( file : in file_type; f : in Faces );

  procedure put ( f : in Array_of_Faces );
  procedure put ( file : in file_type; f : in Array_of_Faces );

end Integer_Faces_of_Polytope_io;
