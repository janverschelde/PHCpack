with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Boolean_Vectors;
with Boolean_Matrices;
with Boolean_Matrices_io;                use Boolean_Matrices_io;
with Standard_Random_Matrices;

procedure ts_pivsel is

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

  function Is_In ( vec : Standard_Natural_Vectors.Vector;
                   vec_end : integer32; nbr : in natural32 )
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if nbr equals vec(k) for k in range vec'first..vec_end,
  --   returns false otherwise.

  begin
    for k in vec'first..vec_end loop
      if vec(k) = nbr
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  procedure Greedy_Search
              ( mat : in Boolean_Matrices.Matrix;
                prm : out Standard_Natural_Vectors.Vector;
                size : out integer32 ) is

  -- DESCRIPTION :
  --   Applies a greedy search to find an initial matching
  --   of the graph defined by the matrix mat.

  -- ON ENTRY :
  --   mat      matrix of a bipartite graph of rows and columns.

  -- ON RETURN :
  --   prm      permutation with the selection of the column;
  --   size     size of the permutation.

  begin
    prm := (prm'range => 0);
    size := 0;
    for row in mat'range(1) loop
      for col in mat'range(2) loop
        if mat(row,col) then
          if not Is_In(prm,size,natural32(col)) then
            size := size + 1;
            prm(row) := natural32(col); exit;
          end if;
        end if;
      end loop;
    end loop;
  end Greedy_Search;

  function Apply_Permutation
             ( mat : Boolean_Matrices.Matrix;
               prm : Standard_Natural_Vectors.Vector )
             return Standard_Natural_Matrices.Matrix is

  -- DESCRIPTION :
  --   The selected columns are marked in the matrix on return by 2,
  --   the rest of the 0 and 1 entries are copied from mat.

    res : Standard_Natural_Matrices.Matrix(mat'range(1),mat'range(2));

  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        if mat(i,j)
         then res(i,j) := 1;
         else res(i,j) := 0;
        end if;
      end loop;
    end loop;
    for i in prm'range loop
      res(i,integer32(prm(i))) := 2;
    end loop;
    return res;
  end Apply_Permutation;

  procedure Initialize_PredFlag
              ( flag : in out Boolean_Vectors.Vector;
                prm : in Standard_Natural_Vectors.Vector;
                prm_size : in integer32 ) is

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

  begin
    for i in flag'range loop
      flag(i) := true;
    end loop;
    for i in prm'first..prm_size loop
      flag(integer32(prm(i))) := false;
    end loop;
  end Initialize_PredFlag;

  procedure Recurse
                ( vtx : in integer32;
                  preds : in out Standard_Natural_VecVecs.VecVec;
                  pred : in out Standard_Natural_Vectors.Vector;
                  flag : in out Boolean_Vectors.Vector;
                  prm : in out Standard_Natural_Vectors.Vector;
                  found : out boolean ) is
 
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

    use Standard_Natural_Vectors;
    u : integer32;

  begin
    if preds(vtx) /= null then
      declare
        lst : constant Standard_Natural_Vectors.Link_to_Vector := preds(vtx);
      begin
        preds(vtx) := null;
        for k in lst'range loop
          u := integer32(lst(k));
          if u > 0 then
            if flag(u) then
              declare
                pu : constant integer32 := integer32(pred(u));
              begin
                flag(u) := false;
                if flag(pu) then
                  prm(vtx) := natural32(u);
                  found := true; return;
                else
                  Recurse(pu,preds,pred,flag,prm,found);
                  if found
                   then prm(vtx) := natural32(u); return;
                  end if;
                end if;
              end;
            end if;
          end if;
        end loop;
      end;
    end if;
    found := false;
  end Recurse;

  function Empty ( b : Boolean_Vectors.Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if b is empty, that is if all elements in b
  --   are false.  Returns false otherwise.

  begin
    for k in b'range loop
      if b(k)
       then return false;
      end if;
    end loop;
    return true;
  end Empty;

  function Empty ( v : Standard_Natural_Vectors.Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if v is empty, that is if all elements in v
  --   are zero.  Returns false otherwise.

  begin
    for k in v'range loop
      if v(k) > 0
       then return false;
      end if;
    end loop;
    return true;
  end Empty;

  procedure Bipartite_Match
              ( graph : in Boolean_Matrices.Matrix;
                prm : in out Standard_Natural_Vectors.Vector;
                prm_size : in out integer32;
                verbose : in boolean := true ) is

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

    dim : constant integer32 := graph'last(1);
    flag : Boolean_Vectors.Vector(1..dim);
    pred : Standard_Natural_Vectors.Vector(1..dim);
    preds : Standard_Natural_VecVecs.VecVec(1..dim);
    unmatched : boolean_Vectors.Vector(1..dim) := (1..dim => false);
    found : boolean;

  begin
    Initialize_PredFlag(flag,prm,prm_size);
    while true loop
      for k in 1..dim loop
        pred(k) := natural32(k);
      end loop;
      for k in 1..prm_size loop
        pred(integer32(prm(k))) := 0;
      end loop;
      if verbose then
        put("The data pred initialized in while :");
        put(pred); new_line;
      end if;
      -- repeatedly extend layering structure by another pair of layers
      while Empty(unmatched) loop
        null;
      end loop;
      -- did we finish layering without finding any alternating paths?
      if Empty(unmatched) then
        null;
      end if;
      for k in unmatched'range loop
        if unmatched(k)
         then recurse(k,preds,pred,flag,prm,found);
        end if;
      end loop;
    end loop;
  end Bipartite_Match;

  procedure Random_Test ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Boolean matrix and runs the algorithms to
  --   solve the pivot selection problem.

    mat : Boolean_Matrices.Matrix(1..dim,1..dim)
        := Standard_Random_Matrices.Random_Matrix(natural32(dim));
    prm : Standard_Natural_Vectors.Vector(1..dim);
    prm_size : integer32;

  begin
    put_line("A random Boolean matrix :"); put(mat);
    Greedy_Search(mat,prm,prm_size);
    put("The permutation : "); put(prm);
    put(" with size : "); put(prm_size,1); new_line;
    if prm_size = mat'last(1) then
      put_line("The greedy search solved the pivot selection problem :");
      put(Apply_Permutation(mat,prm));
    else
      put_line("Running the Hopcroft-Karp algorithm ...");
      Bipartite_Match(mat,prm,prm_size);
    end if;
  end Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension
  --   and then generates a random Boolean matrix.

    dim : integer32 := 0;

  begin
    put("Give the dimension of the problem : "); get(dim);
    Random_Test(dim);
  end Main;

begin
  Main;
end ts_pivsel;
