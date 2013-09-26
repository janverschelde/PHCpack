with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;          use Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;

package Filtered_Points is

-- DESCRIPTION :
--   This package provides some utilities to do the bookkeeping for
--   filtering the points at the end of the paths with respect to
--   the hypersurface equations for solution components of a system.

  function Create ( n,L,d : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..n where the L-th entry is d.
  --   The n-th entry will be L to mark the dimension of the component.

  procedure Update ( fp : in out List; d,i,L : in integer32 );

  -- DESCRIPTION :
  --   One path at level l has converged to the i-th d-dimensional component.

  procedure Write ( file : in file_type; fp : in List );

  -- DESCRIPTION :
  --   Writes a formatted list of filtered points onto the file.

end Filtered_Points;
