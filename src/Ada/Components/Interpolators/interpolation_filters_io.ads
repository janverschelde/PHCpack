with text_io;                           use text_io;
with Interpolation_Filters;             use Interpolation_Filters;

package Interpolation_Filters_io is

-- DESCRIPTION :
--   This package offers input-output facilities for filters.

  procedure get ( f : in out Standard_Filter );
  procedure get ( file : in file_type; f : in out Standard_Filter );
  procedure get ( f : in out Multprec_Filter );
  procedure get ( file : in file_type; f : in out Multprec_Filter );

  -- DESCRIPTION :
  --   Reads a filter from standard input or from file.
  --   The expected format on input is the same as the one on output.

  procedure put ( f : in Standard_Filter );
  procedure put ( file : in file_type; f : in Standard_Filter );
  procedure put ( f : in Multprec_Filter );
  procedure put ( file : in file_type; f : in Multprec_Filter );

  -- DESCRIPTION :
  --   Writes on standard output or on file the dimension, degree 
  --   and centrality, followed by interpolator and projector data.

end Interpolation_Filters_io;
