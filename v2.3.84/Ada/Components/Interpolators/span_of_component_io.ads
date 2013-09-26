with text_io;                            use text_io;
with Span_of_Component;                  use Span_of_Component;

package Span_of_Component_io is

-- DESCRIPTION :
--   This package offers input/output routines for spans of components.

  procedure get ( sp : in out Standard_Span );
  procedure get ( file : file_type; sp : in out Standard_Span );
  procedure get ( sp : in out Multprec_Span );
  procedure get ( file : file_type; sp : in out Multprec_Span );

  -- DESCRIPTION :
  --   Reads in the representation for the span of a component from
  --   standard input or from file, following the output format below.

  procedure put ( sp : in Standard_Span );
  procedure put ( file : file_type; sp : in Standard_Span );
  procedure put ( sp : in Multprec_Span );
  procedure put ( file : file_type; sp : in Multprec_Span );

  -- DESCRIPTION :
  --   Writes the span of a component on standard output or on file,
  --   first writing Ambient_Dimension(sp) and Dimension(sp),
  --   followed by Free_Variables(sp) and Equations(sp).

end Span_of_Component_io;
