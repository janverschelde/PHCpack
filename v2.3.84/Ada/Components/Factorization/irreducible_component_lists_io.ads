with text_io;                            use text_io;
with Irreducible_Component_Lists;        use Irreducible_Component_Lists;

package Irreducible_Component_Lists_io is

-- DESCRIPTION :
--   This package offers input-output facilities for lists of
--   irreducible components.  The first element in both input
--   and output format is the number of elements in the list.

  procedure get ( first,last : in out Standard_Irreducible_Component_List );
  procedure get ( file : in file_type;
                  first,last : in out Standard_Irreducible_Component_List );
  procedure get ( first,last : in out Multprec_Irreducible_Component_List );
  procedure get ( file : in file_type;
                  first,last : in out Multprec_Irreducible_Component_List );

  -- DESCRIPTION :
  --   Reads first the number of elements, then the components one after
  --   the other.  See the format in Irreducible_Components_io.

  procedure put ( L : in Standard_Irreducible_Component_List );
  procedure put ( file : file_type;
                  L : in Standard_Irreducible_Component_List );
  procedure put ( L : in Multprec_Irreducible_Component_List );
  procedure put ( file : file_type;
                  L : in Multprec_Irreducible_Component_List );

  -- DESCRIPTION :
  --   Writes first the number of elements, followed by the components.

-- MINIMAL DATA :
--   The minimal representation of a list of irreducible components consists
--   of the labels to the generic points on each irreducible component.

  procedure get_labels
               ( first,last : in out Standard_Irreducible_Component_List );
  procedure get_labels
               ( file : in file_type;
                 first,last : in out Standard_Irreducible_Component_List );
  procedure get_labels
               ( first,last : in out Multprec_Irreducible_Component_List );
  procedure get_labels
               ( file : in file_type;
                 first,last : in out Multprec_Irreducible_Component_List );

  -- DESCRIPTION :
  --   Reads first the number of components from standard input or from file.
  --   Then reads the labels for each component.  For the input format of the
  --   labels, see the specifications of Irreducible_Components_io.

  procedure put_labels ( L : in Standard_Irreducible_Component_List );
  procedure put_labels ( file : in file_type;
                         L : in Standard_Irreducible_Component_List );

  -- DESCRIPTION :
  --   Writes first the number of components on standard output or on file.
  --   Then writes the labels for each component, according to the format in
  --   the specifications of Irreducible_Components_io.

end Irreducible_Component_Lists_io;
