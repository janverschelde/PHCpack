with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Symmetry_Group;                     use Symmetry_Group;

package Symmetric_Lifting_Functions is

-- DESCRIPTION :
--   This package provides symmetric lifting functions.

  procedure Classify_Orbits
               ( supports : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mix : in Standard_Integer_Vectors.Vector;
                 v,w : in List_of_Permutations; norb : out natural32;
                 orbits : out Arrays_of_Integer_Vector_Lists.Array_of_Lists );

  -- DESCRIPTION :
  --   This procedure provides a symmetric lifting for the given supports,
  --   w.r.t. the symmetry group, represented by the lists of permutations
  --   v and w.

  -- ON ENTRY :
  --   supports  supports of a polynomial system;
  --   mix       type of mixture;
  --   v,w       linear representations of the group actions:
  --               wk P(x) = P(vk x).

  -- ON RETURN :
  --   norb      the number of different orbits;
  --   orbits    the lifted supports, w.r.t. the symmetry group,
  --             as lifting value the orbit number has been used.

  procedure Float_Lift_Orbits
               ( orbits : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 lifting : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Given supports with classified orbits, each orbit will receive
  --   the corresponding lifting value in the vector lifting.

  -- REQUIRED :
  --   lifting has range 1..#orbits

  procedure Integer_Lift_Orbits
               ( orbits : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 lifting : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Given supports with classified orbits, each orbit will receive
  --   the corresponding lifting value in the vector lifting.

  -- REQUIRED :
  --   lifting has range 1..#orbits

  procedure Float_Random_Lift_Orbits
               ( orbits : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 norb : in natural32; lower,upper : in double_float );

  -- DESCRIPTION :
  --   Orbits are given random floating-point lifting values in [lower,upper].
  --   Every point in the same orbit receives the same lifting.

  -- ON ENTRY :
  --   orbits    classified orbits, last entry of vectors is orbit number;
  --   norb      the number of orbits;
  --   lower     lower bound for random number generator;
  --   upper     upper bound for random number generator.

  -- ON RETURN :
  --   orbits    each point in the same orbit has the same random number.

  procedure Integer_Random_Lift_Orbits
               ( orbits : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 norb : in natural32; lower,upper : in integer32 );

  -- DESCRIPTION :
  --   Orbits are given random integer lifting values in [lower,upper].
  --   Every point in the same orbit receives the same lifting.

  -- ON ENTRY :
  --   orbits    classified orbits, last entry of vectors is orbit number;
  --   norb      the number of orbits;
  --   lower     lower bound for random number generator;
  --   upper     upper bound for random number generator.

  -- ON RETURN :
  --   orbits    each point in the same orbit has the same random number.

end Symmetric_Lifting_Functions;
