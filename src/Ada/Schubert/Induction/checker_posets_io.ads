with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Checker_Posets;                     use Checker_Posets;

package Checker_Posets_io is

  procedure Write ( rows,cols : in Vector );
  procedure Write ( file : in file_type; rows,cols : in Vector );

  -- DESCRIPTION :
  --   Writes the rows and the column indices of a white checker
  --   using a bracketed notation on file or standard output.

  procedure Write_Node ( nd : in Node );
  procedure Write_Node ( file : in file_type; nd : in Node );

  -- DESCRIPTION :
  --   Writes the coefficient, row and column indices of the white
  --   checker as stored in the node.

  procedure Write_Parents ( p : in Vector; nd : in Node );

  -- DESCRIPTION :
  --   Writes all parents to the node nd, preceded by the position p
  --   of the black checkers.

  procedure Write_Nodes_in_Poset ( ps : in Poset; i : in integer32 );
  procedure Write_Nodes_in_Poset 
              ( file : in file_type; ps : in Poset; i : in integer32 );

  -- DECRIPTION :
  --   Writes the nodes at level i in the poset ps,
  --   showing the level, the black checkers and then the elaboration
  --   with the red checkers and corresponding root counts.

  procedure Write ( ps : in Poset );

  -- DESCRIPTION :
  --   Writes the poset of black and white checkers.

  procedure Write_Patterns ( ps : in Poset );

  -- DESCRIPTION :
  --   Writes the poset of black and white checkers,
  --   as well as the localization patterns.

  procedure Write_Formal_Sum ( nd : in Link_to_Node );
  procedure Write_Formal_Sum
              ( file : in file_type; nd : in Link_to_Node );

  -- DESCRIPTION :
  --   Writes the formal sum of Littlewood-Richardson coefficients
  --   corresponding to the list of white checkers headed by nd.

  procedure Write_Formal_Sums ( ps : in Poset );

  -- DESRIPTION :
  --   Writes the formal sums of LR coefficients and white brackets
  --   at every level of the poset ps.

  procedure Write_Final_Sum ( ps : in Poset );
  procedure Write_Final_Sum
              ( file : in file_type; ps : in Poset );

  -- DESCRIPTION :
  --   Writes the final form sum at the leaves of the poset.

  procedure Write_Formal_Product ( ps : in Poset );
  procedure Write_Formal_Product
              ( file : in file_type; ps : in Poset );

  -- DESCRIPTION :
  --   Writes the formal product elaborated by the poset,
  --   either to standard output or to file,
  --   which appears at the lefthand side of the formal equation.

  procedure Write_Formal_Equation ( ps : in Poset );
  procedure Write_Formal_Equation
              ( file : in file_type; ps : in Poset );

  -- DESCRIPTION :
  --   Writes the formal equation produced by the poset.

  procedure Write_Node_in_Path
               ( n,k : in integer32; ps : in Poset;
                 path : in Array_of_Nodes; index : in integer32 );
  procedure Write_Node_in_Path
               ( file : in file_type; n,k : in integer32; ps : in Poset;
                 path : in Array_of_Nodes; index : in integer32 );

  -- DESCRIPTION :
  --   Shows all the game information at the current node path(index).

  -- REQUIRED : index >= path'first + 1.

  -- ON ENTRY :
  --   file     by default standard_output;
  --   n        dimension of the ambient space;
  --   k        dimension of the k-plane;
  --   ps       checker poset needed for black checkers;
  --   path     sequence of nodes in a path of white checkers;
  --   index    current position in the path.

  procedure Write_Path_in_Poset
              ( n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32 );

  -- DESCRIPTION :
  --   Shows checker game information for one path through the poset.

  -- ON ENTRY :
  --   n        dimension of the ambient space;
  --   k        dimension of the k-plane;
  --   ps       checker poset needed for black checkers;
  --   path     sequence of nodes in a path of white checkers;
  --   count    counter for path number.

end Checker_Posets_io;
