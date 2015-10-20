with generic_lists;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Intersection_Posets;                use Intersection_Posets;

package Standard_Solution_Posets is

-- DESCRIPTION :
--   An intersection solution poset adds the solutions to the intersection
--   poset to resolve a general Schubert problem.
--   In particular, to any poset node in an intersection poset there is
--   a corresponding solution node.  At each level in the intersection
--   poset there is a corresponding list of solution nodes.
--   The solutions at the nodes are in standard double precision.

  type Solution_Node;
  type Link_to_Solution_Node is access Solution_Node;

  type Solution_Node is record
    sols : Solution_List;
    lpnd : Link_to_Poset_Node;
  end record;

  package Lists_of_Solution_Nodes is new Generic_Lists(Link_to_Solution_Node);
  type Solnode_List is new Lists_of_Solution_Nodes.List;

  type Array_of_Solnode_Lists is
    array ( integer32 range <> ) of Solnode_List;

  type Solution_Poset ( m : integer32 ) is record
    level : integer32;                      -- level of completion
    nodes : Array_of_Solnode_Lists(1..m);   -- defined for 1..level
    last : Array_of_Solnode_Lists(1..m);    -- pointers to last nodes
  end record;

-- CONSTRUCTORS :

  function Create ( pnd : Link_to_Poset_Node ) return Solution_Node;

  -- DESCRIPTION :
  --   Returns the solution node which contains an empty solution list
  --   and the pointer to the poset node as given by pnd.

  function Create ( ips : Intersection_Poset ) return Solution_Poset;

  -- DESCRIPTION :
  --   Returns the solution poset that corresponds to the given
  --   intersection poset.  As no solutions are computed, the level
  --   of completion of the solution poset on return equals zero.
  --   Furthermore, the solution poset is created up to ips.level,
  --   the level of completion of the given intersection poset.

-- SELECTOR :

  function Retrieve ( snl : Solnode_List;
                      rows,cols : Standard_Natural_Vectors.Vector )
                    return Link_to_Solution_Node;

  -- DESCRIPTION :
  --   Returns the pointer to the solution node in the list
  --   with corresponding poset that has as root rows and columns
  --   the given vectors rows and cols.  If there is no such node
  --   in the list, then the null pointer is returned.

-- DESTRUCTORS :

  procedure Clear ( snd : in out Solution_Node );

  -- DESCRIPTION :
  --   Deallocates the memory allocated to the solution list at the node.

  procedure Clear ( snd : in out Link_to_Solution_Node );

  -- DESCRIPTION :
  --   Clears the content of the solution node and frees the pointer
  --   to the solution node.

  procedure Clear ( nodes : in out Solnode_List );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the list of solution nodes.

  procedure Clear ( nodes : in out Array_of_Solnode_Lists );

  -- DESCRIPTION :
  --   Deallocates the memory from all nodes in the array of node lists.

  procedure Clear ( sps : in out Solution_Poset );

  -- DESCRIPTION :
  --   Clears all solution lists in the solution poset.

end Standard_Solution_Posets;
