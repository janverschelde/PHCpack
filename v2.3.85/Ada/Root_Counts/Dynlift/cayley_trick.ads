with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Cayley_Trick is

-- DESCRIPTION :
--   This procedure provides some routines for constructing a mixed
--   subdivision by applicaton of the dynamic lifting algorithm to
--   the Cayley trick.

-- OPTIONS :
--   1. choice of the order of the points
--   2. for vertices or not: allows to add interior points
--   3. with maximum value on the lifting function or not

-- VERSIONS :
--   1. with or without output generics, before flattening or after new cell
--   2. with as result only the mixed subdivision or the whole triangulation

-- BASIC VERSION : WITHOUT OUTPUT GENERICS :

  procedure Dynamic_Cayley 
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  mixsub : out Mixed_Subdivision; numtri : out natural32 );

  procedure Dynamic_Cayley
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  t : in out Triangulation );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the Cayley trick.

  -- ON ENTRY :
  --   n          length of the vectors in the supports;
  --   mix        type of mixture;
  --   supports   lists of integer vector of length n;
  --   order      if true, then the points are already ordered and will
  --              be processed like they occur in the list,
  --              if false, then a random order will be chosen;
  --   inter      if true, then the list may contain interior points,
  --              i.e. points x in conv(l\{x}),
  --              if false, no interior points occur in the list;
  --   maxli      maximum value of the lifting function,
  --              if = 0, then no flattening will be applied,
  --              i.e. there is no maximum lifting value,
  --              if > 0, then no points will be given a lifting value
  --              greater than maxli.

  -- ON RETURN :
  --   lifted     the lifted supports;
  --   mixsub     the mixed subdivision;
  --   numtri     number of simplices in the triangulation of the
  --              auxiliary polytope;
  --   t          triangulation of the auxiliary polytope.

-- EXTENDED VERSIONS : WITH OUTPUT GENERICS

  generic
    with procedure Before_Flattening
                ( mixsub : in out Mixed_Subdivision;
                  lifted : in Array_of_Lists );

    -- DESCRIPTION :
    --   Before flattening, the current mixed subdivision with
    --   the current lists of lifted points are given.

  procedure Dynamic_Cayley_with_Flat
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  mixsub : out Mixed_Subdivision; numtri : out natural32 );

  generic
    with procedure Before_Flattening
                ( mixsub : in out Mixed_Subdivision;
                  lifted : in Array_of_Lists );
  procedure Dynamic_Cayley_with_Flatt
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  t : in out Triangulation );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the Cayley trick.
  --   Before flattening, the generic procedure will be invoked.
  --   The parameters have the same meaning as in the basic version.

  generic
    with procedure Process_New_Cells
                ( mixsub : in out Mixed_Subdivision;
                  i : in integer32; point : in vector );

    -- DESCRIPTION :
    --   After the addition of a new point to the ith component,
    --   this point together with the new mixed cells are returned.
    --   For the initial cell, i=0.

  procedure Dynamic_Cayley_with_New
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  mixsub : out Mixed_Subdivision; numtri : out natural32 );

  generic
    with procedure Process_New_Cells
                ( mixsub : in out Mixed_Subdivision;
                  i : in integer32; point : in vector );
  procedure Dynamic_Cayley_with_Newt
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  t : in out Triangulation );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the Cayley trick.
  --   After each addition of a point, the generic procedure will
  --   be invoked.
  --   The parameters have the same meaning as in the basic version.

  generic
    with procedure Before_Flattening
                ( mixsub : in out Mixed_Subdivision;
                  lifted : in Array_of_Lists );

    -- DESCRIPTION :
    --   Before flattening, the current mixed subdivision with
    --   the current lists of lifted points are given.

    with procedure Process_New_Cells
                ( mixsub : in out Mixed_Subdivision; 
                  i : in integer32; point : in vector );

    -- DESCRIPTION :
    --   After the addition of a new point to the ith component,
    --   this point together with the new simplices are returned.
    --   For the initial cell, i=0.

  procedure Dynamic_Cayley_with_Flat_and_New
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  mixsub : out Mixed_Subdivision; numtri : out natural32 );

  generic
    with procedure Before_Flattening
                ( mixsub : in out Mixed_Subdivision;
                  lifted : in Array_of_Lists );
    with procedure Process_New_Cells
                ( mixsub : in out Mixed_Subdivision;
                  i : in integer32; point : in vector );
  procedure Dynamic_Cayley_with_Flat_and_Newt
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  t : in out Triangulation );


  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the list l.
  --   Before flattening, the first generic procedure will be invoked.
  --   After each addition of a point, the second generic procedure
  --   will be invoked.
  --   The parameters have the same meaning as in the basic version.

end Cayley_Trick;
