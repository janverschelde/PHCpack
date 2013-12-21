with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Monodromy_Group_Actions;           use Monodromy_Group_Actions;

package Monodromy_Group_Actions_io is

-- DESCRIPTION :
--   This package provides output facilities for irreducible components
--   of the i-th dimensional solution set of a polynomial system.

  procedure put ( ic : in Irreducible_Components; i : in integer32 );
  procedure put ( file : in file_type;
                  ic : in Irreducible_Components; i : in integer32 );

  -- DESCRIPTION :
  --   Writes the i-th irreducible component, if it is nonempty.
  --   Nothing happens when Empty(ic,i).

  procedure put ( ic : in Irreducible_Components );
  procedure put ( file : in file_type; ic : in Irreducible_Components );

  -- DESCRIPTION :
  --   Writes the irreducible components on standard output or on file.

  procedure put_labels ( ic : in Irreducible_Components );
  procedure put_labels ( file : in file_type;
                         ic : in Irreducible_Components );

  -- DESCRIPTION :
  --   Writes the irreducible components on standard output or on file,
  --   in a format also shared by the rest of the program...

end Monodromy_Group_Actions_io;
