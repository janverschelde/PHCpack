with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;

package Dynamic_Mixed_Subdivisions is

-- DESCRIPTION :
--   This package provides implementations of the dynamic lifting
--   algorithm for the construction of a mixed subdivision.

-- OPTIONS :
--   1. Choice of the order of the points;
--   2. For vertices or not, allows to add interior points;
--   3. With flattening or not, allows maximum lifting value.

-- DATA STRUCTURES :

  type Face_Structure is record
    l,last : List;      -- l contains the lifted points of the faces,
                        -- last is a pointer to the last element of l
    t : Triangulation;  -- triangulation of the points in l
    f : Faces;          -- k-faces of conv(l)
  end record;
  type Face_Structures is array ( integer32 range <> ) of Face_Structure;

-- BASIC VERSIONS : WITHOUT OUTPUT GENERICS :

  procedure Dynamic_Lifting
                ( n : in integer32; mix : in Vector;
                  points : in Array_of_Lists;
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision; 
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the points.

  -- ON ENTRY :
  --   n          dimension before lifting;
  --   mix        type of mixture;
  --   points     tuple of lists of points of length n;
  --   order      if true, then the points are already ordered and will
  --              be processed like they occur in the list,
  --              if false, then a random order will be chosen;
  --   inter      if true, then the list may contain interior points,
  --              that are points x in conv(l\{x}),
  --              if false, no interior points occur in the list;
  --   conmv      if true, then online checks on zero contributions will 
  --              be made, otherwise this will not happen;
  --   maxli      maximum value of the lifting function,
  --              if = 0, then no flattening will be applied,
  --              if > 0, then no points will be given a lifting value
  --              larger than maxli, eventually flattening will be applied;
  --   mixsub     regular mixed subdivision of the lifted points;
  --   fs         the face structures of the lower hulls of the lifted
  --              points, necessary to re-start the algorithm efficiently;
  --   nbsucc     number of successful face-face combinations;
  --   nbfail     number of unsuccessful face-face combinations.

  -- ON RETURN :
  --   mixsub     regular mixed subdivision of the lifted points;
  --   fs         face structures of the updated lower hulls;
  --   nbsucc     updated number of successful face-face combinations;
  --   nbfail     updated number of unsuccessful face-face combinations.

-- EXTENDED VERSIONS : WITH OUTPUT GENERICS

  generic
    with procedure Before_Flattening
                ( mixsub : in out Mixed_Subdivision; fs : in Face_Structures );

    -- DESCRIPTION :
    --   Before flattening, the current mixed subdivision and face
    --   structures are given as parameter to the procedure above.

  procedure Dynamic_Lifting_with_Flat
                ( n : in integer32; mix : in Vector;
                  points : in Array_of_Lists;
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the points.
  --   Before flattening, the generic procedure will be invoked.
  --   The parameters have the same meaning as in the basic version.

  generic
    with procedure Process_New_Cells
                ( mixsub : in out Mixed_Subdivision;
                  i : in integer32; point : in vector );

    -- DESCRIPTION :
    --   After the addition of a new point to the ith component,
    --   this lifted point together with the new mixed cells are returned.
    --   After the computation of the initial cell, i=0.

  procedure Dynamic_Lifting_with_New
                ( n : in integer32; mix : in Vector;
                  points : in Array_of_Lists;
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the points.
  --   After each addition of a point, the generic procedure will
  --   be invoked.
  --   The parameters have the same meaning as in the basic version.

  generic
    with procedure Before_Flattening
                ( mixsub : in out Mixed_Subdivision; fs : in Face_Structures );

    -- DESCRIPTION :
    --   Before flattening, the current mixed subdivision with
    --   the current lists of lifted points are given.

    with procedure Process_New_Cells
                ( mixsub : in out Mixed_Subdivision;
                  i : in integer32; point : in vector );

    -- DESCRIPTION :
    --   After the addition of a new point to the ith component,
    --   this lifted point together with the new simplices are returned.
    --   After the computation of the initial cell, i=0.

  procedure Dynamic_Lifting_with_Flat_and_New
                ( n : in integer32; mix : in Vector;
                  points : in Array_of_Lists;
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the list l.
  --   Before flattening, the first generic procedure will be invoked.
  --   After each addition of a point, the second generic procedure
  --   will be invoked.
  --   The parameters have the same meaning as in the basic version.

end Dynamic_Mixed_Subdivisions;
