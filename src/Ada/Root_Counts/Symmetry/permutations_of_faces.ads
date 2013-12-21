with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Permutations;                       use Permutations;

package Permutations_of_Faces is

-- DESCRIPTION :
--   Provides elementary permutation operations on faces of polytopes.
--   Each time, the second operation concerns a lifted face and leaves
--   the lifting invariant.

  function Invariant ( f : Face; p : Permutation ) return boolean;
  function Invariant_Lifted ( f : Face; p : Permutation ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the face is invariant under the permutation.

  function Permute ( f : Face; p : Permutation ) return Face;
  function Permute_Lifted ( f : Face; p : Permutation ) return Face;

  -- DESCRIPTION :
  --   Permutations the coordinates of the points which span the face.

  function Permutable ( f1,f2 : Face ) return boolean;
  function Permutable_Lifted ( f1,f2 : Face ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there exists a permutation of face f1 into face f2.

  function Permutable ( f1 : Face; f2 : Faces ) return boolean;
  function Permutable_Lifted ( f1 : Face; f2 : Faces ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the face f1 is permutable with one of the faces in f2.

end Permutations_of_Faces;
