with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;

package Standard_Dynamic32_Triangulations is

-- DESCRIPTION :
--   This package offers some implementations of the dynamic lifting
--   algorithm applied to the unmixed case.

-- OPTIONS :
--   1. choice of the order of the points
--   2. for vertices or not: allows to add interior points
--   3. with maximum value on the lifting function or not

-- BASIC VERSION : WITHOUT OUTPUT GENERICS :

  procedure Dynamic_Lifting
                ( L : in List; order,inter : in boolean;
                  maxli : in integer32; lifted,lifted_last : in out List;
                  t : in out Triangulation );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the points in l.

  -- ON ENTRY :
  --   L          list of integer vectors;
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
  --              greater than maxli;
  --   lifted     points that already have been lifted;
  --   lifted_last is pointer to the last element of the list lifted;
  --   t          contains a triangulation of the lifted points.

  -- ON RETURN :
  --   lifted     the lifted points;
  --   t          a regular triangulation of the points in L.

-- EXTENDED VERSIONS : WITH OUTPUT GENERICS

  generic
    with procedure Before_Flattening
                ( t : in Triangulation; lifted : in List );

    -- DESCRIPTION :
    --   Before flattening, the current triangulation with
    --   the current list of lifted points is given.

  procedure Dynamic_Lifting_with_Flat
                ( L : in List; order,inter : in boolean;
                  maxli : in integer32; lifted,lifted_last : in out List;
                  t : in out Triangulation );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the list L.
  --   Before flattening, the generic procedure will be invoked.
  --   The parameters have the same meaning as in the basic version.

  generic
    with procedure Process_New_Simplices
                ( t : in Triangulation; point : in vector );

    -- DESCRIPTION :
    --   After the addition of a new point, this lifted point together
    --   with the new simplices are returned.

  procedure Dynamic_Lifting_with_New
                ( L : in List; order,inter : in boolean;
                  maxli : in integer32; lifted,lifted_last : in out List;
                  t : in out Triangulation );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the list L.
  --   After each addition of a point, the generic procedure will be invoked.
  --   The parameters have the same meaning as in the basic version.

  generic
    with procedure Before_Flattening
                ( t : in Triangulation; lifted : in List );

    -- DESCRIPTION :
    --   Before flattening, the current triangulation with
    --   the current list of lifted points is given.

    with procedure Process_New_Simplices
                ( t : in Triangulation; point : in vector );

    -- DESCRIPTION :
    --   After the addition of a new point, this lifted point together
    --   with the new simplices are returned.

  procedure Dynamic_Lifting_with_Flat_and_New
                ( L : in List; order,inter : in boolean;
                  maxli : in integer32; lifted,lifted_last : in out List;
                  t : in out Triangulation );
  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the list L.
  --   Before flattening, the first generic procedure will be invoked.
  --   After each addition of a point, the second generic procedure
  --   will be invoked.
  --   The parameters have the same meaning as in the basic version.

end Standard_Dynamic32_Triangulations;
