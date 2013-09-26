with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Index_Tree_LP;                       use Index_Tree_LP;

package Zero_Index_Tree is

-- DESCRIPTION :
--   A "zero_index_tree" defines the data structures for the specific
--   case of level 0 in the index tree of LP problems.

  type L0IdxNode;
  type Link_to_L0IdxNode is access L0IdxNode;

  type L0IdxNode is record      -- structure to define linked list
    idx : integer32;            -- index of the node in the list is a number
    R : Link_to_IndexNode;      -- linked list of LPdata
    D : Link_to_L0IdxNode;      -- points to the next item in the list
  end record;

  type L0_IML;
  type Link_to_L0_IML is access L0_IML;

  type L0_IML is record         -- index tree of LP problems at level zero
    L0head : Link_to_L0IdxNode; 
    L0curr : Link_to_L0IdxNode; -- Current IT node ptr, for Add(-)
    L0prev : Link_to_L0IdxNode; -- parent of curr, for Add(-)
    curr : Link_to_IndexNode;   -- Current IT node ptr, for Add(-)
    prev : Link_to_IndexNode;   -- parent of curr, for Add(-)
    LP1 : Link_to_LPPL; -- dummy head of the LP address link (used in level 1)
  end record;

  procedure L0IdxNode_Init ( p : in out Link_to_L0IdxNode; i : in integer32 );

  -- DESCRIPTION :
  --   Initializes the linked list of nodes with LPdata.
  --
  -- ON ENTRY :
  --   p         memory allocated for one structure of type L0IdxNode;
  --   i         value for p.idx.
  --
  -- ON RETURN :
  --   p         p.idx == i and the pointers p.R and p.D are 0.

  procedure L0_IML_Init ( p : in out Link_to_L0_IML );

  -- DESCRIPTION :
  --   Allocates memory for p.L0head and p.LP1, initializing both.
  --   The pointers p0.L0curr and p0.L0prev are set to p.L0head.

  procedure L0_Migrate
               ( p : in out Link_to_L0_IML;
                 inp : in out Link_to_IndexNode;
                 status : out integer32 );

  -- DESCRIPTION :
  --   Migrates the data of the index tree at level zero with the node inp.
  --
  -- ON ENTRY :
  --   p         an index tree with LP problems at level zero;
  --   inp       pointer to an index node.
  --
  -- ON RETURN :
  --   p         p.L0head.D is updated and p.L0prev is cleared;
  --   inp       inp.S now has the value of p.L0prev.R;
  --             p.L0prev.R is assigned to inp.S;
  --   status    0 if p.L0head.D is empty, 1 otherwise.

  procedure L0_FindInR
               ( p : in out Link_to_L0_IML;
                 IDX : in integer32; found : out boolean );

  -- DESCRIPTION :
  --   Does horizontal search in p for a node with label idx equal to IDX.
  --
  -- ON ENTRY :
  --   p         index tree with LP problems at level zero;
  --   IDX       label to an index node.
  --
  -- ON RETURN :
  --   p         p.cur.idx equals IDX if found equals true;
  --   found     false if there is no index node in the list starting
  --               with p.prev.S with the given label IDX, otherwise
  --             true if here is an index node with this index

  procedure L0_FindInD
               ( p : in out Link_to_L0_IML;
                 IDX : in integer32; found : out boolean );

  -- DESCRIPTION :
  --   Does vertical search in p for a node with label idx equal to IDX.
  --
  -- ON ENTRY :
  --   p         index tree with LP problems at level zero;
  --   IDX       label to an index node.
  --
  -- ON RETURN :
  --   p         p.L0curr.idx equals IDX if found equals true;
  --   found     false if there is no index node in the list starting
  --               with p.L0prev.D with the given label IDX, otherwise
  --             true if there is an index node with this index.

  procedure L0_Add1
               ( p : in out Link_to_L0_IML; n : in integer32; 
                 J : in Standard_Integer_Vectors.Link_to_Vector; 
                 d : in integer32; 
                 I : in Standard_Integer_Vectors.Link_to_Vector; 
                 X : in Standard_Floating_Vectors.Link_to_Vector;
                 A : in Standard_Floating_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Updates the index tree at level zero with data for a LP problem.
  --
  -- REQUIRED :
  --   The n entries in J are sorted in ascending order.
  --
  -- ON ENTRY :
  --   p         current index tree at level zero of LP problems.
  --   n         dimension of the vector J;
  --   J         index vector of n entries, sorted in ascending order;
  --   d         dimension of the LP problem;
  --   I         index vector of dimension d with constraints involved;
  --   X         solution vector of dimension d;
  --   A         basis inverse, matrix of dimension d.
  --
  -- ON RETURN :
  --   p         updated index tree at level zero of LP problems.

  procedure L0_Add2 
               ( p : in out Link_to_L0_IML;
                 J : in Standard_Integer_Vectors.Link_to_Vector;
                 d : in integer32;
                 I : in Standard_Integer_Vectors.Link_to_Vector;
                 X : in Standard_Floating_Vectors.Link_to_Vector;
                 A : in Standard_Floating_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Updates the index tree at level zero with data for a LP problem,
  --   for only two points.
  --
  -- REQUIRED :
  --   The two entries in J are sorted in ascending order.
  --
  -- ON ENTRY :
  --   p         current index tree at level zero of LP problems;
  --   J         two number sorted in ascending order;
  --   d         dimension of the LP problem;
  --   I         index vector of dimension d with constraints involved;
  --   X         solution vector of dimension d;
  --   A         basis inverse, matrix of dimension d.
  --
  -- ON RETURN :
  --   p         updated index tree at level zero of LP problems.

  procedure L0_IML_Del ( p : in out Link_to_L0_IML );

  -- DESCRIPTION :
  --   Deletes the index tree at level zero,
  --   calls the routine L0_free below.

  procedure L0_Free ( li : in out Link_to_L0_IML );

  -- DESCRIPTION :
  --   Deallocation of all the LPdata in the index tree,
  --   called by the other routine L0_IML_Del.

  procedure Clear ( p : in out Link_to_L0IdxNode );
  procedure Clear ( p : in out Link_to_L0_IML );

  -- DESCRIPTION :
  --   Deallocation of the memory referred to by the pointer.

end Zero_Index_Tree;
