with text_io;                              use text_io;
with Irreducible_Components;               use Irreducible_Components;

package Irreducible_Components_io is

-- DESCRIPTION :
--   This package offers input-output facilities for irreducible 
--   solution components.  The formats for input and output are matched.
--   See the formats for interpolation filters and the span of component
--   for details.

  procedure get ( c : out Standard_Irreducible_Component );
  procedure get ( c : out Multprec_Irreducible_Component );
  procedure get ( file : in file_type;
                  c : out Standard_Irreducible_Component );
  procedure get ( file : in file_type;
                  c : out Multprec_Irreducible_Component );

  procedure put ( c : in Standard_Irreducible_Component );
  procedure put ( c : in Multprec_Irreducible_Component );
  procedure put ( file : in file_type;
                  c : in Standard_Irreducible_Component );
  procedure put ( file : in file_type; 
                  c : in Multprec_Irreducible_Component );

-- MINIMAL DATA :
--   The minimal representation of an irredubible component consists of
--   the labels of the generic points that lie on the component.

  procedure get_labels ( c : out Standard_Irreducible_Component );
  procedure get_labels ( file : in file_type;
                         c : out Standard_Irreducible_Component );
  procedure get_labels ( c : out Multprec_Irreducible_Component );
  procedure get_labels ( file : in file_type;
                         c : out Multprec_Irreducible_Component );

  -- DESCRIPTION :
  --   A new component is created from the labels read from standard input
  --   or from file.  The format of the label is "d : l(1) .. l(d)", where
  --   d stands for the degree of the component and l(i) is the label of
  --   the i-th generic point.

  procedure put_labels ( c : in Standard_Irreducible_Component );
  procedure put_labels ( file : in file_type;
                         c : in Standard_Irreducible_Component );

  -- DESCRIPTION :
  --   Writes the label of c on standard output or on file, according
  --   to the same input format.  An empty label is "0 :".
  --   The writing of the label concludes with the end of line symbol.
 
end Irreducible_Components_io;
