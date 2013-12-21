with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Brackets;                          use Brackets;
with Localization_Posets;               use Localization_Posets;

package Localization_Poset_Strings is

-- DESCRIPTION :
--   This package provides operations to convert a localization poset
--   to compute the Pieri root count into a string.

  function Bracket_to_String ( b : Bracket ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the bracket b.

  function Node_to_String
            ( top,bottom : Bracket; roco : natural32 ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the top, bottom bracket
  --   of a node along with its root count.

  function Nodes_to_String ( nd : Node ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the node
  --   and all its siblings.

  function Poset_to_String ( p : Array_of_Nodes ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the poset used
  --   for the Pieri root count.

end Localization_Poset_Strings;
