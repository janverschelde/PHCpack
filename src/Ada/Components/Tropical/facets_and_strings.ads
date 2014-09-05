with Standard_Integer_Vectors;
with Multprec_Integer_Vectors;
with Multprec_Integer_Matrices;
with Multprec_Lattice_3d_Facets;
with Multprec_Lattice_4d_Facets;

package Facets_and_Strings is

-- DESCRIPTION :
--   This package provides functions to convert the data of a facet
--   into a string representation.

  function write ( v : Standard_Integer_Vectors.Vector ) return string;
  function write ( v : Multprec_Integer_Vectors.Vector ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the vector as a tuple,
  --   a sequence of numbers separated by commas, enclosed by square brackets.
  --   For example, as "[-9, 3, 4]", suitable for parsing by Python.

  function write ( A : Multprec_Integer_Matrices.Matrix;
                   f : Multprec_Lattice_3d_Facets.Facet_in_3d ) return string;
  function write ( A : Multprec_Integer_Matrices.Matrix;
                   f : Multprec_Lattice_4d_Facets.Facet_in_4d ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the facet f, spanned by points
  --   with coordinates in the matrix A.

end Facets_and_Strings;
