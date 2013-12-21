with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;

package Index_Tree_LP is

-- DESCRIPTION :
--    This package defines the general data structures for an index tree
--    of LP problems, with unused level 0, treated separtedly by the data
--    structures and operations in the package zero_index_tree.

  type LPdata;
  type Link_to_LPdata is access LPdata;

  type LPdata is record    -- LP data = indices to active constraints,
                           --     feasible solution, and basis inverse
    dim : integer32;       -- size of the vectors xxx, JJJ, and matrix INV
    xxx : Standard_Floating_Vectors.Link_to_Vector;    
                           -- vector of range 0..dim-1, solution vector
    INV : Standard_Floating_Matrices.Link_to_Matrix;
                           -- matrix of dimension dim, basis inverse
    JJJ : Standard_Integer_Vectors.Link_to_Vector;
                           -- index vector for the solution to the LP problem
  end record;

  type LPPL;
  type Link_to_LPPL is access LPPL;
  type Array_of_LPPL is array ( integer32 range <> ) of Link_to_LPPL;
  type Link_to_Array_of_LPPL is access Array_of_LPPL;

  type LPPL is record      -- LP pointer linked list
    addr : Link_to_LPdata; -- information of one node in the list
    next : Link_to_LPPL;   -- pointer to the next item in the list
  end record;

  type IndexNode;
  type Link_to_IndexNode is access IndexNode;
  type Array_of_IndexNodes is
    array ( integer32 range <> ) of Link_to_IndexNode;
  type Link_to_Array_of_IndexNodes is access Array_of_IndexNodes;
  
  type IndexNode is record -- node in linked list of LPdata
    idx : integer32;       -- index of the node is just a number
    info : Link_to_LPdata; -- information for an LP problem
    S : Link_to_IndexNode; -- points to next element in linked list
  end record;

  type IT_LP;
  type Link_to_IT_LP is access IT_LP;

  type IT_LP is record       -- index tree for LP problems
    MaxLevels : integer32;   -- maximal #levels for IT[]
    CurLevel : integer32;    -- index to the current level
                     -- the next 5 items are arrays of dimension MaxLevels
    DIM : Standard_Integer_Vectors.Link_to_Vector;
                           -- dimension of LP in each level
    NP : Standard_Integer_Vectors.Link_to_Vector; 
                           -- NP[L] = #nodes in IT[L], including fixed node
    cell : Standard_Integer_Vectors.Link_to_Vector;
                           -- indices in horizontal branch of tree when full
    InSpt : Standard_Integer_Vectors.Link_to_Vector;
                           -- indices to the supports
    minNP : Standard_Integer_Vectors.Link_to_Vector;
                           -- minimal #nodes in IT[L], excluding fixed node
    LP : Link_to_Array_of_LPPL;
                           -- array of dimension MaxLevels of LP problems
    LPlast : Link_to_LPPL; -- points to the last valid node in LP[CurLevel]
    IT,last : Link_to_Array_of_IndexNodes;
                           -- array of size MaxLevels, IT[0] is not used
                           -- last[] points to the last valid node in IT[]
    curr : Link_to_IndexNode;  -- pointer to the current IT node
    prev : Link_to_IndexNode;  -- pointer to the previous IT node
  end record;

  procedure LPdata_Init
               ( p : in out Link_to_LPdata; n : in integer32;
	         J : in Standard_Integer_Vectors.Link_to_Vector;
	         x : in Standard_Floating_Vectors.Link_to_Vector;
	         A : in Standard_Floating_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Initializes the data in p with the given parameters.
  --
  -- ON ENTRY :
  --   p         memory allocated to hold one structure of type LPdata;
  --   n         dimension of the LP problem;
  --   J         index vector for the solution of the LP problem;
  --   x         vector of range 0..n-1, a solution vector;
  --   A         matrix of dimension n, basis inverse used in LP.
  --
  -- ON RETURN :
  --   p         copies of J, x, and A are stored in p.

  procedure IndexNode_Init
               ( p : in out Link_to_IndexNode; i : in integer32;
                 ptr : in Link_to_LPdata );

  -- DESCRIPTION :
  --   Initializes the index node p with i and the LP data.
  --
  -- ON ENTRY :
  --   p         memory allocated to hold one IndexNode structure;
  --   i         value for p.idx;
  --   ptr       data of the LP problem.
  --
  -- ON RETURN :
  --   p         the pointer fields in p are assigned to i and ptr,
  --             no deep copies of the fields in ptr are made.

  procedure LPPL_Init
               ( p : in out Link_to_LPPL; A : in Link_to_LPdata;
                 N : in Link_to_LPPL );

  -- DESCRIPTION :
  --   Initializes the linked list of LP problems with the LP data in A
  --   and pointer to the next item in the list N.
  --
  -- ON ENTRY :
  --   p         memory allocated for one structure;
  --   A         data for one LP problem;
  --   N         pointer information to the next item in the list.
  --
  -- ON RETURN :
  --   p         structure with its two fields assigned to using A and N,
  --             no deep copies of the fields in p are made.

  procedure IT_LP_Init
               ( p : in out Link_to_IT_LP; nSpt : in integer32;
                 mtype : in Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Initialization of the index tree of linear programs p.
  --
  -- ON ENTRY :
  --   p         memory allocated for one structure;
  --   nSpt      number of different supports;
  --   mtype     type of mixture of the supports,
  --             given by an array of range 0..nSpt-1.
  --
  -- ON RETURN :
  --   p         allocated memory for all its members.

  function IT_IsEmpty ( p : Link_to_IT_LP ) return boolean;

  -- DESCRIPTION :
  --   Returns the outcome of the test p.CurLevel < 1.

  function IT_IsFull ( p : Link_to_IT_LP ) return boolean;

  -- DESCRIPTION :
  --   Returns the outcome of the test p.CurLevel+1 >= p.MaxLevels.

  function IT_Level ( p : Link_to_IT_LP ) return integer32;

  -- DESCRIPTION :
  --   Returns the value of p.CurLevel,
  --   i.e.: the value of the current level in the index tree.

  function IT_CurLPdim ( p : Link_to_IT_LP ) return integer32;

  -- DESCRIPTION :
  --   Returns the value of p.DIM[p.CurLevel],
  --   i.e.: the dimension of the LP problems at the current level.

  function IT_Cell ( p : Link_to_IT_LP )
                   return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns p.cell, the pointer containing the indices of
  --   a horizontal branch.  Note that p.cell(0) is not being used.

  function IT_CurSptIdx ( p : Link_to_IT_LP ) return integer32;

  -- DESCRIPTION :
  --   Returns the value of p.InSpt[p.CurLevel],
  --   i.e.: the index of the support concerning the current node.

  function IT_MinNumPt ( p : Link_to_IT_LP ) return integer32;

  --  DESCRIPTION :
  --    Returns the value of p.minNP[p.CurLevel],
  --    i.e.: the minimal #points required for moving to the next level.

  function IT_NumPt ( p : Link_to_IT_LP ) return integer32;

  -- DESCRIPTION :
  --   Returns the value of p.NP[p.CurLevel],
  --   i.e.: the #points in the current vertical branch of the index tree.

  function IT_FixedIdxNdPtr ( p : Link_to_IT_LP ) return Link_to_IndexNode;

  -- DESCRIPTION :
  --   For p nonempty, returns p.IT[p.CurLevel],
  --   i.e.: the current head of the vertical branch in the tree.

  procedure IT_StepBack ( p : in out Link_to_IT_LP );

  -- DESCRIPTION :
  --   Sets the number of points at the current level to zero
  --   and decreases the value of p.CurLevel by one.

  procedure IT_Find ( p : in out Link_to_IT_LP; IDX : in integer32;
                      found : out boolean );

  -- DESCRIPTION :
  --   Searches for a node in p with index equal to IDX.
  --
  -- ON ENTRY :
  --   p         index tree with LP problems;
  --   IDX       label to an index node.
  --
  -- ON RETURN :
  --   p         index tree with updated pointers;
  --   found     false if there is no index node with index equal to IDX,
  --             true otherwise, in which case: p.curr.idx == IDX.

  procedure IT_ResetCurLevelTo1
               ( p : in out Link_to_IT_LP; frst : out Link_to_IndexNode );

  -- DESCRIPTION :
  --   Resets the value of the current level in the index tree to 1.
  --
  -- ON ENTRY :
  --   p         an index tree for LP problems.
  --
  -- ON RETURN :
  --   p         the values of p.CurLevel and p.NP[1] are set to 1, 
  --             and p.prev is set to p.IT[1], other entries are freed;
  --   frst      p.IT[1].

  procedure IT_RenewNP1 ( p : in out Link_to_IT_LP );

  -- DESCRIPTION :
  --   Increases the count of p.NP[1] by the number of elements in p.IT[1].
  --   The values of p.last[1] and p.cell[1] are changed as well.

  procedure IT_NextLevel ( p : in out Link_to_IT_LP; rcase : out integer32 );

  -- DESCRIPTION :
  --   Prepares the index tree to move to the next level.
  --
  -- ON ENTRY :
  --   p         index tree for LP problems.
  --
  -- ON RETURN :
  --   p         updated index tree, with swapped pointers;
  --   rcase     0 if p.CurLevel is equal to the p.MaxLevels, or
  --             if the number of nodes at the current level is fewer
  --                than the minimal number of nodes at that level,
  --             1 otherwise.

  procedure IT_Add1
               ( p : in out Link_to_IT_LP; n : in integer32;
                 J : in Standard_Integer_Vectors.Link_to_Vector;
                 nn : in integer32;
                 JJ : in Standard_Integer_Vectors.Link_to_Vector;
                 X : in Standard_Floating_Vectors.Link_to_Vector;
                 A : in Standard_Floating_Matrices.Link_to_Matrix );
                 
  -- DESCRIPTION :
  --   Updates the index tree of LP problems with data for one LP problem,
  --   for an array of indices to points.
  --
  -- ON ENTRY :
  --   p         index tree of LP problems;
  --   n         dimension of the array J;
  --   J         labels of points;
  --   nn        dimension of the LP problem;
  --   JJ        array of dimension nn with indices
  --             of the constraints involved;
  --   X         solution vector of dimension nn;
  --   A         basis inverse, matrix of dimension nn.
  --
  -- ON RETURN :
  --   p         updated index tree of LP problems.

  procedure IT_Add2
               ( p : in out Link_to_IT_LP; oneidx,nn : in integer32;
                 JJ : in Standard_Integer_Vectors.Link_to_Vector;
                 X : in Standard_Floating_Vectors.Link_to_Vector;
                 A : in Standard_Floating_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Updates the index tree of LP problems with data for one LP problem,
  --   for just one point.
  --
  -- ON ENTRY :
  --   p         index tree of LP problems;
  --   oneidx    index to the point in the support;
  --   nn        dimension of the LP problem;
  --   JJ        array of dimension nn with indices
  --             of the constraints involved;
  --   X         solution vector of dimension nn;
  --   A         basis inverse, matrix of dimension nn.
  --
  -- ON RETURN :
  --   p         updated index tree of LP problems.

  procedure IT_LP_DEL ( p : in out Link_to_IT_LP );
 
  -- DESCRIPTION :
  --  This is the main function to deallocate all memory occupied by p,
  --  it calls the other two IT_Free routines below.

  procedure IT_FreeIT ( p : in out Link_to_IT_LP );

  -- DESCRIPTION :
  --   Deallocation of the index tree structures in p.IT.

  procedure IT_FreeLP ( p : in out Link_to_IT_LP );
 
  -- DESCRIPTION :
  --    Deallocation of all the LP data in p.LP.

  procedure Clear ( p : in out Link_to_IndexNode );
  procedure Clear ( p : in out Link_to_LPdata );
  procedure Clear ( p : in out Link_to_LPPL );

  -- DESCRIPTION :
  --   Deallocation of the memory referred to by the pointer.

  procedure Clear ( p : in out Array_of_LPPL );
  procedure Clear ( p : in out Link_to_Array_of_LPPL );
  procedure Clear ( p : in out Array_of_IndexNodes );
  procedure Clear ( p : in out Link_to_Array_of_IndexNodes );

  -- DESCRIPTION :
  --   Deallocation of the memory occupied by the arrays.

end Index_Tree_LP;
