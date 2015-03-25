with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Integer_Matrices;          use Multprec_Integer_Matrices;
with Multprec_Lattice_3d_Facets;
with Multprec_Lattice_4d_Facets;

package Multprec_Giftwrap_Container is

-- DESCRIPTION :
--   This package stores the list of facets computed with the
--   giftwrapping method for convex hulls as persistent storage
--   for interfacing to C and Python.

-- CONSTRUCTOR :

  procedure Create ( A : in Matrix );

  -- DESCRIPTION :
  --   Given in the columns of A, the coordinates of points,
  --   computes and stores the list of facets of the convex hull of A,
  --   either in 3d or in 4d, depending on the number of rows in A.
  --   If the span of the points in the matrix A do not have the
  --   expected dimension equal to A'last(1), then the list of facets 
  --   will be empty and Number_of_{3,4}_Facets will return zero.

  -- REQUIRED :
  --   A'range(1) = 1..3 or A'range(1) = 1..4.

  procedure Store_String ( s : in string );

  -- DESCRIPTION :
  --   Stores the string s in the container.

-- SELECTORS :

  function Point_Configuration_in_3d return Link_to_Matrix;

  -- DESCRIPTION :
  --   Returns a pointer to the matrix that contains in its columns
  --   the coordinates of the points in 3-space.

  function Point_Configuration_in_4d return Link_to_Matrix;

  -- DESCRIPTION :
  --   Returns a pointer to the matrix that contains in its columns
  --   the coordinates of the points in 4-space.

  function Number_of_3d_Facets return natural32;

  -- DESCRIPTION :
  --   Returns the number of facets in the list of the 3-dimensional hull.

  function Number_of_4d_Facets return natural32;

  -- DESCRIPTION :
  --   Returns the number of facets in the list of the 4-dimensional hull.

  function Facet_3d_Data
             ( k : natural32 )
             return Multprec_Lattice_3d_Facets.Link_to_3d_Facet;

  -- DESCRIPTION :
  --   Returns a pointer to the data representation to a facet,
  --   assumed to be available as the k-th facet in the list of
  --   computed facets for the 3-dimensional convex hull.
  --   The k corresponds to the identification number of the facet.
  --   If k is out of the range 0..Number_of_3d_Facets-1,
  --   then the null pointer will be returned.

  function Facet_4d_Data
             ( k : natural32 )
             return Multprec_Lattice_4d_Facets.Link_to_4d_Facet;

  -- DESCRIPTION :
  --   Returns a pointer to the data representation to a facet,
  --   assumed to be available as the k-th facet in the list of
  --   computed facets for the 4-dimensional convex hull.
  --   The k corresponds to the identification number of the facet.
  --   If k is out of the range 0..Number_of_4d_Facets-1,
  --   then the null pointer will be returned.

  function Retrieve_String return string;

  -- DESCRIPTION :
  --   Returns the string stored in the container;

-- DESTRUCTORS : 

  procedure Clear_3d;

  -- DESCRIPTION :
  --   Clears the data structures stored for the 3-dimensional convex hull.

  procedure Clear_4d;

  -- DESCRIPTION :
  --   Clears the data structures stored for the 3-dimensional convex hull.

  procedure Clear_String;

  -- DESCRIPTION :
  --   Clears the string stored in the container.

  procedure Clear;

  -- DESCRIPTION :
  --    Clears all internal data structures.

end Multprec_Giftwrap_Container;
