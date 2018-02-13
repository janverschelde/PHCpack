with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
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
            prm(size) := natural32(col);
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
    put("The permutation : ");
    put(prm(prm'first..prm_size)); new_line;
    if prm_size = mat'last(1) then
      put_line("The greedy search solved the pivot selection problem :");
      put(Apply_Permutation(mat,prm));
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
