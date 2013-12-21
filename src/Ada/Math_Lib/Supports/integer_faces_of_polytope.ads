with Generic_Lists;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Integer_Faces_of_Polytope is

-- DESCRIPTION :
--   This package offers an abstraction for manipulating the k-faces
--   of an n-dimensional polytope.

-- DATA STRUCTURES :

  type Face is new Link_to_VecVec;
    -- points that span the face, range 0..k, with k = dimension of face

  type Face_Array is array ( integer32 range <> ) of Face;

  package Lists_of_Faces is new Generic_Lists(Face);
  type Faces is new Lists_of_Faces.List;

  type Array_of_Faces is array ( integer32 range <> ) of Faces;

-- CONSTRUCTORS :

  function Create ( k,n : integer32; p : List ) return Faces;
  function Create ( k,n : integer32; p : List; x : Vector ) return Faces;

  -- DESCRIPTION :
  --   The input for this routine is a list of points, defining a polytope 
  --   in n-dimensional space.  On return, a list of lists of points will be 
  --   given, defining k-faces of the polytope, each spanned by k+1 points.
  --   When x is provided, then only those k-faces will be returned that
  --   contain the vector x.  In the latter case, x must belong to p.

  function Create_Lower ( k,n : integer32; p : List ) return Faces;
  function Create_Lower ( k,n : integer32; p : List; x : Vector ) return Faces;

  -- DESCRIPTION :
  --   Only the k-faces of the lower hull will be generated.
  --   When x is added as parameter, only the k-faces that contain x
  --   will be returned.  In the latter case, x must belong to p.

  procedure Construct ( first : in out Faces; fs : in Faces );

  -- DESCRIPTION :
  --   All the faces in fs will be constructed to the front of first.

  procedure Copy ( f1 : in Face; f2 : in out Face );

  -- DESCRIPTION :
  --   Copies the face f1 to the face f2, after clearing f2.

  procedure Deep_Copy ( f1 : in Faces; f2 : in out Faces );

  -- DESCRIPTION :
  --   Makes a deep copy from the list f1 to the list f2.

-- SELECTORS :

  function Is_Equal ( f1,f2 : Face ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both faces are spanned by the same vertices.
  --   Note that the order of the vertices in the vector do not matter.

  function Is_In ( f : Face; x : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the face contains the vector x.

  function Is_In ( fs : Faces; f : Face ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the face already belongs to the list of faces.

  function Extract_Faces ( fs : Faces; x : Vector ) return Faces;

  -- DESCRIPTION :
  --   Returns a list of faces that contain x.
  --   The list on return shares the faces with the original list fs.

-- DESTRUCTORS :

  procedure Deep_Clear ( f : in out Face );
  procedure Deep_Clear ( fa : in out Face_Array );
  procedure Deep_Clear ( fs : in out Faces );
  procedure Deep_Clear ( afs : in out Array_of_Faces );
  procedure Shallow_Clear ( f : in out Face );
  procedure Shallow_Clear ( fa : in out Face_Array );
  procedure Shallow_Clear ( fs : in out Faces );
  procedure Shallow_Clear ( afs : in out Array_of_Faces );

  -- DESCRIPTION :
  --   All allocated memory space will be freed.  A deep clear destroys the 
  --   whole structure, while a shallow clear only destroys the structures 
  --   especially set up for the faces.

end Integer_Faces_of_Polytope;
