with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Symmetry_Group;                     use Symmetry_Group;

package Faces_of_Symmetric_Polytopes is

-- DESCRIPTION :
--   This package contains some routines to construct the tuple of
--   generating faces of symmetric polytopes.
--   When the group representation is not given as parameter,
--   full permutation symmetry is assumed.
--   If `Lifted'-routines are called, then the permutations will leave
--   the lifting value invariant.

-- ON A FACE : group * faces -> invariant subgroup

  function Stabilizer ( v : List_of_Permutations; f : Face )
                      return List_of_Permutations;
  function Stabilizer_Lifted ( v : List_of_Permutations; f : Face )
                             return List_of_Permutations;

  -- DESCRIPTION :
  --   Returns those permutations that leave the face invariant.

-- ON FACES : group * faces -> invariant faces

  function Invariant_Faces ( v : List_of_Permutations;
                             f : Faces ) return Faces;

  function Invariant_Lifted_Faces ( v : List_of_Permutations;
                                    f : Faces ) return Faces;

  -- DESCRIPTION :
  --   Returns those faces which are invariant under the permutations.
  --   To check this for the full permutation group, the list of
  --   generators of the group should be supplied.

-- ON FACES : group * faces -> generated faces

  function Generated_Faces ( v : List_of_Permutations; f : Faces )
                           return Faces;

  function Generated_Lifted_Faces
                           ( v : List_of_Permutations; f : Faces )
                           return Faces;

  -- DESCRIPTION :
  --   Returns those faces that generate all faces in f.
  --   For the full permutation group, supply the generators of the group.

-- ON FACES : group * faces -> generators of faces

  function Generating_Faces ( f : Faces ) return Faces;
  function Generating_Faces ( v : List_of_Permutations; f : Faces )
                            return Faces;

  function Generating_Lifted_Faces ( f : Faces ) return Faces;
  function Generating_Lifted_Faces
                            ( v : List_of_Permutations; f : Faces )
                            return Faces;

  -- DESCRIPTION :
  --   Returns those faces that generate all faces in f.

-- ON TUPLES OF FACES : group * faces -> invariant faces

  function Invariant_Faces ( v : List_of_Permutations;
                             af : Array_of_Faces ) return Array_of_Faces;

  function Invariant_Lifted_Faces ( v : List_of_Permutations;
                                    af : Array_of_Faces ) return Array_of_Faces;

  -- DESCRIPTION :
  --   Returns for each component those faces which are invariant under the
  --   permutations. To check this for the full permutation group, the list
  --   of generators of the group should be supplied.
  --   It is assumed that the tuple is invariant under v.

-- ON TUPLES OF FACES : group * faces -> generators of faces

  function Generating_Faces ( af : Array_of_Faces ) return Array_of_Faces;
  function Generating_Faces ( v : List_of_Permutations; af : Array_of_Faces )
                            return Array_of_Faces;
  function Generating_Faces ( v,w : List_of_Permutations; af : Array_of_Faces )
                            return Array_of_Faces;

  function Generating_Lifted_Faces
                            ( af : Array_of_Faces ) return Array_of_Faces;
  function Generating_Lifted_Faces
                            ( v : List_of_Permutations; af : Array_of_Faces )
                            return Array_of_Faces;
  function Generating_Lifted_Faces
                            ( v,w : List_of_Permutations; af : Array_of_Faces )
                            return Array_of_Faces;

  -- DESCRIPTION :
  --   Returns the generating faces of the tuple.  When w is left out,
  --   it is assumed that the tuple is invariant.  When both v and w are not
  --   supplied, then it is assumed that the tuple is equi-invariant w.r.t.
  --   the full permutation group.

end Faces_of_Symmetric_Polytopes;
