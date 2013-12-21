with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers; 
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_VecVecs;           use Standard_Natural_VecVecs;

package Checker_Posets is

-- DESCRIPTION :
--   A checker poset records the data for one "checker game".
--   A checker game defines a deformation of one moving flag,
--   starting at the opposite of a fixed flag (at the root) and
--   ending at the same position of the fixed flag.
--   The coefficients obtained at the leaves of the poset are
--   also called Littlewood-Richardson coefficients.

-- DATA STRUCTURES :

  type Node;
  type Link_to_Node is access Node;

  type Node ( k : integer32 ) is record
    coeff : Natural_Number;             -- Littlewood-Richardson coefficient
    rows : Vector(1..k);                -- row positions on board
    cols : Vector(1..k);                -- column positions
    first_parent : Link_to_Node;        -- first parent of the node
    stay_child : Link_to_Node;          -- child with same location
    swap_child : Link_to_Node;          -- child with swapped rows
    next_sibling : Link_to_Node;        -- next sibling at same level
  end record;

  type Array_of_Nodes is array ( integer32 range <> ) of Link_to_Node;
  type Link_to_Array_of_Nodes is access Array_of_Nodes;

  type Poset is record
    black : Link_to_VecVec;             -- location of black checkers
    white : Link_to_Array_of_Nodes;     -- location of white checkers
  end record;

-- CONSTRUCTORS :

  function Specializing_Moves ( n : integer32 ) return VecVec;

  -- DESCRIPTION :
  --   Returns a vector describing the moves of the black checkers,
  --   of range 1..Checker_Moves.Number_of_Moves(n),
  --   starting with the identity permutation (the shortest word)
  --   and ending at its reverse permutation (the longest word).
  --   These moves specialize a general flag into the fixed identity.

  function Generalizing_Moves ( n : integer32 ) return VecVec;

  -- DESCRIPTION :
  --   Returns a vector describing the moves of the black checkers,
  --   of range 1..Checker_Moves.Number_of_Moves(n),
  --   starting with the reverse permutation (the longest word)
  --   and ending at its identity permutation (the shortest word).
  --   These moves generalize the fixed identity to a general flag.

  function Create ( r,c : Vector ) return Node;

  -- DESCRIPTION :
  --   Creates the node with rows and columns in r and c respectively. 

  function Create ( n,k : integer32; root : Node ) return Poset;

  -- DESCRIPTION :
  --   Creates the poset for k-planes in n-space starting with
  --   white checkers as given in the root and all black checkers
  --   on the diagonal.

  function Create ( n : integer32; r,c : Vector ) return Poset;

  -- DESCRIPTION :
  --   Creates the poset for a board of n checkers,
  --   with the black checkers on the diagonal and white checkers
  --   in rows r and columns indicated by c.

  function Create ( n : integer32; cff : Natural_Number;
                    r,c : Vector ) return Poset;

  -- DESCRIPTION :
  --   Creates the poset for a board of n checkers, where the node
  --   at the root has coefficient equal to cff.  The position of
  --   the white checkers are defined by row and column indices
  --   respectively in r and c.

  procedure Add_Multiplicity ( ps : in out Poset; m : in Natural_Number );

  -- DESCRIPTION :
  --   Adds m to every coefficient in the poset.

-- SELECTORS :

  function Root_Rows ( ps : in Poset ) return Vector;

  -- DESCRIPTION :
  --   Returns the rows indicating the position of the white checkers
  --   at the root of the poset.

  function Root_Columns ( ps : in Poset ) return Vector;

  -- DESCRIPTION :
  --   Returns the columns indicating the position of the white checkers
  --   at the root of the poset.

  function Equal ( nd1,nd2 : Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the rows and columns of both nodes are the same,
  --   return false otherwise.

  function Position ( first_nd : Node; nd : Node ) return integer32;

  -- DESCRIPTION :
  --   Returns 0 if the node nd does not occur in the list of nodes,
  --   starting at the current node first_nd,
  --   otherwise returns the index k of the same sibling as nd.

  function Retrieve ( ps : Poset; i,j : integer32 ) return Link_to_Node;

  -- DESCRIPTION :
  --   Returns the j-th node at level i in the poset ps.
  --   If there is no j-th node at level i in ps, then null is returned.

  procedure Retrieve_Leaf ( ps : in Poset; cols : in Vector;
                            ind : out integer32; lnd : out Link_to_Node );

  -- DESCRIPTION :
  --   Returns the leaf node in the poset ps whose column position of 
  --   the white checkers are the same as the given cols.
  --   If no such leaf exists, lnd is null on return and also ind is 0,
  --   otherwise ind on return equals the position of the leaf in ps.

  function Is_Stay_Child ( parent,child : Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if Equal(parent.stay_child.all,child),
  --   returns false otherwise.

  function Is_Swap_Child ( parent,child : Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if Equal(parent.swap_child.all,child),
  --   returns false otherwise.

  generic
    with procedure Process_Parent ( pnd : in Link_to_Node );
  procedure Enumerate_Parents ( nd : in Node );

  -- DESCRIPTION :
  --   Enumerates all parents of the given node nd,
  --   calling the procedure Process_Parent for each parent node to nd.

  function Number_of_Parents ( nd : Node ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of parents of the given node nd.

  function Multiplicity_of_Parents ( nd : Node ) return Natural_Number;

  -- DESCRIPTION :
  --   Returns the sum of the multiplicities of all parents of the node nd.

  function Retrieve_Parent ( nd : Node; k : integer32 ) return Link_to_Node;

  -- DESCRIPTION :
  --   Returns a pointer to the k-th parent of the given node nd.
  --   The pointer on return will be null if there is no k-th parent.

  function Degree_of_Freedom ( ps : Poset ) return natural32;

  -- DESCRIPTION :
  --   Returns the degree of freedom of a localization pattern
  --   at the root of the poset.

  generic
    with procedure Process_Path 
           ( nds : in Array_of_Nodes; cont : out boolean );

    -- DESCRIPTION :
    --   The procedure is called each time a new path is found.

    -- REQUIRED : nds'range = ps.white'range.
    
    -- ON ENTRY :
    --   nds      the nodes along the path, starting at a leaf
    --            and ending at the root of the poset.

    -- ON RETURN :
    --   cont     if true on return, then the enumeration continues,
    --            otherwise the enumeration stops.

  procedure Enumerate_Paths_in_Poset ( ps : in Poset );

  -- DESCRIPTION :
  --   Enumerates all paths starting at all leaves in the poset.
  --   Process_Path is called each time a new path is found.

-- DESTRUCTORS :

  procedure Clear ( lnd : in out Link_to_Node );
  procedure Clear ( arrlnd : in out Array_of_Nodes );
  procedure Clear ( arrlnd : in out Link_to_Array_of_Nodes );
  procedure Clear ( ps : in out Poset );

  -- DESCRIPTION :
  --   Deallocation of memory.

end Checker_Posets;
