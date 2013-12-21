with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;

package Cell_Stack is

-- DESCRIPTION :
--   This package defines a stack of mixed cells
--   and defines the operations on this stack.

  type Cell;
  type Link_to_Cell is access Cell;

  type Cell is record
    idx : Standard_Integer_Vectors.Link_to_Vector;
                           -- labels of the points which span the mixed cell
    next : Link_to_Cell;   -- pointer to the next cell
  end record;

  type CellStack is record
    size : integer32;      -- total #points which span a mixed cell
    count : integer32;     -- number of elements in the stack
    top : Link_to_Cell;    -- top of the stack
    cur : Link_to_Cell;    -- pointer to the current cell in the stack
  end record;

  procedure Cell_Init1 ( c : in out Link_to_Cell );

  -- DESCRIPTION :
  --   Sets c.next to zero.
         
  procedure Cell_Init2 
               ( c : in out Link_to_Cell; n : in integer32;
                 J : in Standard_Integer_Vectors.Link_to_Vector;
                 ptr : in Link_to_Cell );

  -- DESCRIPTION :
  --   Defines the content of the cell with the indices in J.
  --
  -- ON ENTRY :
  --   c         memory allocated for one cell;
  --   n         number of indices of J;
  --   J         labels to the points spanning the mixed cell;
  --   ptr       pointer for the the next cell.
  --
  -- ON RETURN :
  --   c         c.idx contains the values of J and
  --             c.next has the same value as ptr.

  procedure Cs_Init ( cs : in out CellStack; n : in integer32 );

  -- DESCRIPTION :
  --   Initializes the stack cs to contain mixed cells spanned by n points.
     
  procedure Cs_Push ( cs : in out CellStack;
                      J : in Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --  Pushes a mixed cell defined by the indices in J to the stack.

  procedure Cs_Pop ( cs : in out CellStack );

  -- DESCRIPTION :
  --   Pops the top element from the cell stack.
      
  procedure Cs_Next ( cs : in out CellStack; okay : out integer32 );

  -- DESCRIPTION :

  -- ON ENTRY :
  --   cs        a stack of mixed cells;

  -- ON RETURN :
  --   cs        if the next cell to the current cell is not empty,
  --             then the pointer to the current cell is set to the
  --             address of the next cell;
  --   okay      0 if the next cell to the current cell is empty,
  --             otherwise 1.

  function Cs_Cur ( cs : CellStack )
                  return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the content of the current cell.

  procedure Cs_Top ( cs : in out CellStack );

  -- DESCRIPTION :
  --   Assigns the current cell pointer to the top of the stack.
      
  function Cs_IsEmpty ( cs : CellStack ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the cell stack is empty, false otherwise.

  function Cs_Count ( cs : CellStack ) return integer32;

  -- DESCRIPTION :
  --   Returns the value of cs.count, the #cells in the stack.
      
  procedure Csi ( cs : in out CellStack; i : in integer32;
                  J : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Returns the content of the i-th cell in the stack.

  -- ON ENTRY :
  --   cs        a stack of mixed cells;
  --   i         number of a cell, starting to count at zero.

  -- ON RETURN :
  --   cs        cs.cur has been used to traverse the stack,
  --             and cs.cur.idx is J on return;
  --   J         null if there is no i-th cell, otherwise J
  --             contains the labels of the i-th mixed cell.

  procedure Cs_Del ( cs : in out CellStack );

  -- DESCRIPTION :
  --   Pops all the cells of the stack, deallocating all memory.

  procedure Clear ( c : in out Link_to_Cell );

  -- DESCRIPTION :
  --   Deallocation of the pointer to a cell.

end Cell_Stack;
