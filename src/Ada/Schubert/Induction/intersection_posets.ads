with generic_lists;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Checker_Posets;                     use Checker_Posets;

package Intersection_Posets is

-- DESCRIPTION :
--   An intersection poset stores several checker posets organized
--   successively to count the solutions to many intersection conditions.

-- DATA STRUCTURES :

  type Poset_Node;
  type Link_to_Poset_Node is access Poset_Node;

  type Poset_Node is record
    ps : Poset;                           -- records one checker game
    first_parent : Link_to_Poset_Node;    -- pointer to first parent
    first_child : Link_to_Poset_Node;     -- pointer to first child
  end record;

  package Lists_of_Poset_Nodes is new Generic_Lists(Link_to_Poset_Node);
  type Poset_List is new Lists_of_Poset_Nodes.List;

  type Array_of_Poset_Lists is array ( integer32 range <> ) of Poset_List;

  type Intersection_Poset ( m : integer32 ) is record
    level : integer32;                    -- level of completion
    nodes : Array_of_Poset_Lists(1..m);   -- defined for 1..level
    last : Array_of_Poset_Lists(1..m);    -- pointers to last nodes
  end record;

-- CONSTRUCTORS :

  function Create ( ps : Poset ) return Poset_Node;

  -- DESCRIPTION :
  --   Returns one node in the poset initialized with the given poset ps
  --   and with all links set to zero.

  function Create ( m : integer32; ps : Poset ) return Intersection_Poset;

  -- DESCRIPTION :
  --   Returns an intersection poset of size m,
  --   initialized at the root with the given poset.

  procedure Intersect ( ips : in out Intersection_Poset;
                        pnd : in Link_to_Poset_Node; w : in Vector;
                        silent : in boolean ); 

  -- DESCRIPTION :
  --   Starting at the poset node pnd at ips.level, new poset nodes
  --   are added to the next level in the intersection poset to deal
  --   with the intersection condition imposed by w.
  --   When silent is true, no information is printed,
  --   otherwise one can track the progress of the resolution.

  -- REQUIRED : 1 <= ips.level < m and pnd belongs to ips.nodes(ips.level).

  procedure Intersect ( ips : in out Intersection_Poset; w : in Vector;
                        silent : in boolean );

  -- DESCRIPTION :
  --   Starting at the leaves of all poset nodes at ips.level,
  --   a new level in the poset is created to deal with the intersection
  --   condition imposed by the k-vector w.
  --   If silent then no extra information is printed,
  --   otherwise the outcome of each bracket intersection is shown.

  -- REQUIRED : 1 <= ips.level < m.

-- SELECTORS :

  function Degree_of_Freedom 
             ( ips : Intersection_Poset; k : integer32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the degree of freedom of the poset nodes at level k.

  function Degree_of_Freedom ( ips : Intersection_Poset ) return natural32;

  -- DESCRIPTION :
  --   Returns the degree of freedom of the nodes in the poset at
  --   the current level, i.e.: at ips.level.

  function Final_Sum ( ips : Intersection_Poset ) return Natural_Number;

  -- DESCRIPTION :
  --   Returns the final sum as the root count of the intersection poset.

  function Retrieve ( pl : Poset_List; k : integer32 )
                    return Link_to_Poset_Node;

  -- DESCRIPTION :
  --   Returns the k-th poset node in the list pl,
  --   or null if there is no k-th list.

  procedure Retrieve ( pl : in Poset_List; rows,cols : in Vector;
                       isin : out boolean; pnd : out Link_to_Poset_Node );

  -- DESCRIPTION :
  --   Seeks for the poset in the list whose root is defined by the
  --   given rows and columns of the white checkers.

  -- ON ENTRY :
  --   pl       list of poset nodes, possibly empty;
  --   rows     rows of white checkers;
  --   cols     columns of white checkers.

  -- ON RETURN :
  --   isin     true if there is a poset node in pl whose root
  --            has white checkers at the positions at rows and cols,
  --            false is no such poset node is contained in pl;
  --   pnd      is undefined if isin is false, otherwise, if isin is true
  --            on return, then pnd points to the node whose root of the
  --            poset has white checkers at the given position.

  function Is_Parent ( parent,child : Poset ) return boolean;

  -- DESCRIPTION :
  --   Returns true if one of the leaves of parent has columns of its
  --   white checker equals to the rows of the root white checker of child.
  --   Returns false otherwise.

  function Is_Parent ( pnd,cnd : Poset_Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if pnd is a parent node of the node cnd, false otherwise.
  --   A node P is a parent of node C if one of the leafs of P has columns of
  --   its white checker equal to the rows at the root white checker of C.

  function Is_Child ( child,parent : Poset ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the rows at the root white checkers of child
  --   correspond to the columns of one of the leaves of parent.
  --   Returns false otherwise.

  function Is_Child ( cnd,pnd : Poset_Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if cnd is a child node of the node pnd, false otherwise.
  --   A node C is a child of node P if the rows at the root white checkers
  --   of C correspond to the columns of one of the leafs of the node P.

  generic
    with procedure Process_Parent ( pnd : in Link_to_Poset_Node );
  procedure Enumerate_Parents ( pl : in Poset_List; nd : in Poset_Node );

  -- DESCRIPTION :
  --   Enumerates all parents of the node nd, from the list pl.
  --   Each time a parent of nd is found, Process_Parent is called.

  function Number_of_Parents
             ( pl : Poset_List; nd : Poset_Node ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of parents from the list pl of the given node nd.

  function Retrieve_Parent
             ( pl : Poset_List; nd : Poset_Node; k : integer32 )
             return Link_to_Poset_Node;

  -- DESCRIPTION :
  --   Returns the k-th parent of nd from the list pl.

-- DESTRUCTORS : 

  procedure Clear ( ps : in out Poset_Node );
  procedure Clear ( ps : in out Link_to_Poset_Node );

  -- DESCRIPTION :
  --   Deallocation of memory occupied by ps.

  procedure Clear ( pl : in out Poset_List );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the list.

  procedure Clear ( apl : in out Array_of_Poset_Lists );

  -- DESCRIPTION :
  --   Deallocation of all lists in the array apl.

  procedure Clear ( ps : in out Intersection_Poset );

  -- DESCRIPTION :
  --   Deallocation of all levels in the poset ps.

end Intersection_Posets;
