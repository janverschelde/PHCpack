with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Boolean_Vectors;
with Boolean_Matrices;

package Pivot_Selection is

-- DESCRIPTION :
--   Given the pattern of a sparse matrix, represented by zeros and ones,
--   decide if the determinant of the sparse matrix is zero or not.
--   Replacing the ones in the pattern by random numbers and then compute
--   the determinant of that matrix is one way to solve this problem,
--   The problem is equivalent to the maximum matching problem in
--   a bipartite graph, matching each row to a unique column.
--   For this problem, this package offers exact and efficient algorithms.

-- PREPROCESSING :
--   If the matrix is not so sparse, then a greedy search, which runs
--   linearly in the dimension of the matrix, often finds the maximum
--   matching and solves the pivot selection problem.

  function Is_In ( vec : Standard_Natural_Vectors.Vector;
                   vec_end : integer32; nbr : in natural32 )
                 return boolean;

  -- DESCRIPTION :
  --   Returns true if nbr equals vec(k) for k in range vec'first..vec_end,
  --   returns false otherwise.

  procedure Greedy_Search
              ( mat : in Boolean_Matrices.Matrix;
                prm : out Standard_Natural_Vectors.Vector;
                size : out integer32 );

  -- DESCRIPTION :
  --   Applies a greedy search to find an initial matching
  --   of the graph defined by the matrix mat.

  -- ON ENTRY :
  --   mat      adjacency matrix of a bipartite graph of rows and columns.

  -- ON RETURN :
  --   prm      permutation with the selection of the columns,
  --            zero entries indicate the column was not selected,
  --            in case of zero entries, the matching is not maximal;
  --   size     the size of the permutation equals the number of nonzero
  --            entries in prm, if size = mat'last(1), then the matching
  --            represented by prm is maximal.

-- FORD-FULKERSON :
--   The Ford-Fulkerson algorithm finds a maximum matching with a cost
--   proportional to the number of vertices times the number of edges.

  procedure dfs4bpm
              ( graph : in Boolean_Matrices.Matrix;
                vtx : in integer32;
                seen : in out Boolean_Vectors.Vector;
                match : in out Standard_Integer_Vectors.Vector;
                found : out boolean );

  -- DESCRIPTION :
  --   A depth first search based recursive procedure that returns
  --   found as true if a matching for vertex vtx is possible.

  -- ON INPUT :
  --   graph    adjacency matrix of a bipartite graph;
  --   vtx      a vertex in the graph;
  --   seen     marks the visited vertices;
  --   match    current matching in the graph.

  -- ON RETURN :
  --   seen     updated list of visited vertices;
  --   match    updated matching;
  --   found    true if a matching for vertex vtx is possible,
  --            false otherwise.

  procedure maxBPM
              ( graph : in Boolean_Matrices.Matrix;
                match : out Standard_Integer_Vectors.Vector;
                size : out natural32 );

  -- DESCRIPTION :
  --   Computes the maximum matching in a bipartite graph.

  -- ON ENTRY :
  --   graph    the adjacency matrix of a bipartite graph.

  -- ON RETURN :
  --   match    match(i) is the row assigned to column i,
  --            if -1 then the column is unassigned;
  --   size     number of assigned columns.

end Pivot_Selection;
