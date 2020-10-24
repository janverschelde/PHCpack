with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Natural_VecVecs;
with Standard_Natural_Matrices;
with Boolean_Vectors;
with Boolean_Matrices;

package Test_Pivot_Selection is

-- DESCRIPTION :
--   Interactive development of code for the pivot selection problem.
--   This problem takes on input a Boolean matrix,
--   which represents the sparsity in the coefficient matrix
--   of a linear system: 0 for a zero, 1 for a nonzero coefficient
--   and asks for selection of nonzero pivots,
--   exactly one of every row and one of every column.  
--   If the maximum selection of pivots equals the dimension
--   of the matrix, then the linear system has a unique solution.
--   This pivot selection problem corresponds to the maximum
--   matching problem in a bipartite graph, which can be solved
--   efficiently by the Hopcroft-Karp algorithm.
--   Code for the Ford-Fulkerson algorithm is in the package
--   Pivot_Selection.

  procedure Initialize_PredFlag
              ( flag : in out Boolean_Vectors.Vector;
                prm : in Standard_Natural_Vectors.Vector;
                prm_size : in integer32 );

  -- DESCRIPTION :
  --   The pred-flag data structure signals when for a vertex u,
  --   pred(u) is in the first layer.
  --   The pred data structure gives for every vertex the neighbor 
  --   in the previous layer for each vertex not in the matching
  --   permutation prm, in the left vertex set of the graph.

  -- ON ENTRY :
  --   prm        permutation represents the column selection;
  --   prm_size   size of the permutation.

  -- ON RETURN :
  --   flag       initialized pred-flag data structure,
  --              for matched vertices vtx, flag(vtx) is false,
  --              otherwise flag(vtx) is true.

  procedure Recurse ( vtx : in integer32;
                      preds : in out Standard_Natural_VecVecs.VecVec;
                      pred : in out Standard_Natural_Vectors.Vector;
                      flag : in out Boolean_Vectors.Vector;
                      prm : in out Standard_Natural_Vectors.Vector;
                      found : out boolean );
 
  -- DESCRIPTION :
  --   Recursively search backward through layers to find alternating paths.
  --   The recursion returns true if found path, returns false otherwise.

  -- ON ENTRY :
  --   vtx        a vertex;
  --   preds      gives for preds[v] the neighbors in the previous layer 
  --              for v in the right vertex set of the bipartite graph;
  --   pred       gives for pred[u] the neighbor in the previous 
  --              layer for u in the left vertex set of the graph;
  --   prm        permutation to represent the matching.

  -- ON RETURN :
  --   preds      updated preds data structure;
  --   pred       updated pred data structure;
  --   prm        updated permuation of the matched vertices;
  --   found      true if a path was found, false otherwise.

  function Empty ( b : Boolean_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if b is empty, that is if all elements in b
  --   are false.  Returns false otherwise.

  function Empty ( v : Standard_Natural_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if v is empty, that is if all elements in v
  --   are zero.  Returns false otherwise.

  procedure Bipartite_Match
              ( graph : in Boolean_Matrices.Matrix;
                prm : in out Standard_Natural_Vectors.Vector;
                prm_size : in integer32;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Computes the maximum cardinality matching of a bipartite graph.

  -- ON ENTRY :
  --   graph    the rows of the matrix are adjacency lists to columns,
  --            in a bipartite graph between rows and columns;
  --   prm      permutation initialized by a greedy search;
  --   prm_size is the size of the permutation obtained greedily
  --            and this size is expected to be smaller than graph'last(1);
  --   verbose  flag for extra output during the computation.

  -- ON RETURN :
  --   prm      updated permutation;
  --   prm_size is the updated size of the permutation.

-- CODE for testing :

  function Apply_Permutation
             ( mat : Boolean_Matrices.Matrix;
               prm : Standard_Natural_Vectors.Vector;
               check : boolean := true )
             return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   The selected columns are marked in the matrix on return by 2,
  --   the rest of the 0 and 1 entries are copied from mat.
  --   If check is on, then an error will be raised if a 2 will
  --   be placed at a spot where there is no one.

  procedure Test ( mat : in Boolean_Matrices.Matrix );

  -- DESCRIPTION :
  --   Tests the algorithms for the pivot selection
  --   on the matrix mat.

  procedure Random_Test ( nbrows,nbcols : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Boolean matrix and runs the algorithms to
  --   solve the pivot selection problem.

  procedure Test_Input ( nbrows,nbcols : in integer32 );

  -- DESCRIPTION :
  --   Given in nbrows the number of rows and in nbcols the number of
  --   columns, prompts the user for a 0/1 matrix of the given dimensions
  --   and runs the algorithms for the pivot selection problem.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a dimension and for the type of matrix,
  --   random or user given.

end Test_Pivot_Selection;
