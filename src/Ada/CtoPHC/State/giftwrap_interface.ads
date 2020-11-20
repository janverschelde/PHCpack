with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Giftwrap_Interface is

  function Giftwrap_Planar_Hull
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Parses the input data of a[0] integers in b as a string,
  --   constructs the planar convex hull and returns the vertices
  --   and inner normals in their string representation.

  -- ON ENTRY :
  --   a       the number of characters in the string representation
  --           of the point configuration;
  --   b       as many integers as there are characters in the string
  --           representation of the point configuration;
  --   vrblvl  is the verbose level.

  function Giftwrap_Spatial_Hull
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Parses the string representation for a matrix and initializes the 
  --   gift wrapping container with the list of facets.

  -- ON ENTRY :
  --   a       the number of characters in the string representation
  --           of the point configuration;
  --   b       as many integers as there are characters in the string
  --           representation of the point configuration;
  --   vrblvl  is the verbose level.

  function Giftwrap_Number_of_Facets
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
 
  -- DESCRIPTION :
  --   Returns in b the number of facets of a convex hull in 3-space 
  --   or 4-space, depending on the value of a.

  -- ON ENTRY :
  --   a       either 3 or 4, for the dimension;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of facets in the spatial hull.

  function Giftwrap_String_of_Facet
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Extracts from a and b the dimension of the space and the
  --   number of a facet, computes the string representation of
  --   the facet and returns its number of characters.

  -- ON ENTRY :
  --   a       either 3 or 4, for the dimension;
  --   b       index of the facet;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of characters in the string
  --           representation of the facet.

  function Giftwrap_String_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the length of a string of a support set.
  --   The string is stored internally and is returned
  --   by the Giftwrap_String_of_Support.

  -- ON ENTRY :
  --   a       is the index of the Laurent polynomial stored in the
  --           container for Laurent systems in double precision,
  --           must be between 1 and the number of polynomials;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of characters in the string representation
  --           computed in Giftwrap_String_of_Facet.

  function Giftwrap_String_of_Support
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string stored in the giftwrap container.

  -- REQUIRED :
  --   The Giftwrap_String_Size has been called,
  --   or the Giftwrap_String_of_Facet.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       characters in the string representation of a support set,
  --           as many as the value returned by Giftwrap_String_Size,
  --           or Giftwrap_String_of_Facet.

  function Giftwrap_String_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates the string computed by Giftwrap_String_of_Facet.

  function Giftwrap_Laurent_Initial_Form
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes initial form of a Laurent system in double precision,
  --   and replaces the system in the container by its initial form.

  -- ON ENTRY :
  --   a       the number of variables;
  --   b       an inner normal to define the initial form with;
  --   vrblvl  is the verbose level.

  function Giftwrap_3d_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocation of the storage for the list of facets of a convex hull
  --   in 3-space.  The verbose level is given in vrblvl.

  function Giftwrap_4d_Clear 
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocation of the storage for the list of facets of a convex hull
  --   in 4-space.  The verbose level is given in vrblvl.

end Giftwrap_Interface;
