with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Boolean_Matrices_io;                use Boolean_Matrices_io;
with Standard_Random_Matrices;
with Pivot_Selection;

package body Test_Pivot_Selection is

  procedure Initialize_PredFlag
              ( flag : in out Boolean_Vectors.Vector;
                prm : in Standard_Natural_Vectors.Vector;
                prm_size : in integer32 ) is
  begin
    for i in flag'range loop
      flag(i) := true;
    end loop;
    for i in prm'first..prm_size loop
      flag(integer32(prm(i))) := false;
    end loop;
  end Initialize_PredFlag;

  procedure Recurse ( vtx : in integer32;
                      preds : in out Standard_Natural_VecVecs.VecVec;
                      pred : in out Standard_Natural_Vectors.Vector;
                      flag : in out Boolean_Vectors.Vector;
                      prm : in out Standard_Natural_Vectors.Vector;
                      found : out boolean ) is
 
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

  begin
    for k in b'range loop
      if b(k)
       then return false;
      end if;
    end loop;
    return true;
  end Empty;

  function Empty ( v : Standard_Natural_Vectors.Vector ) return boolean is
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
                prm_size : in integer32;
                verbose : in boolean := true ) is

    dim : constant integer32 := graph'last(1);
    flag : Boolean_Vectors.Vector(1..dim);
    pred : Standard_Natural_Vectors.Vector(1..dim);
    preds : Standard_Natural_VecVecs.VecVec(1..dim);
    unmatched : constant boolean_Vectors.Vector(1..dim) := (1..dim => false);
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
     -- while Empty(unmatched) loop
     --   null;
     -- end loop;
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

-- CODE for testing :

  function Apply_Permutation
             ( mat : Boolean_Matrices.Matrix;
               prm : Standard_Natural_Vectors.Vector;
               check : boolean := true )
             return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(mat'range(1),mat'range(2));
    idx : integer32;

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
      if not check then
        res(i,integer32(prm(i))) := 2;
      else
        idx := integer32(prm(i));
        if res(i,idx) = 1 then
          res(i,idx) := 2;
        else
          put("Problem at row "); put(i,1);
          put(" and column "); put(idx,1);
          put_line(" : no edge!");
        end if;
      end if;
    end loop;
    return res;
  end Apply_Permutation;

  procedure Test ( mat : in Boolean_Matrices.Matrix ) is

    dim : constant integer32 := mat'last(1);
    prm : Standard_Natural_Vectors.Vector(1..dim);
    prm_size : integer32;
    selected : Standard_Integer_Vectors.Vector(1..dim);
    selected_size : natural32;

  begin
    put_line("An adjacency matrix :"); put(mat);
    Pivot_Selection.Greedy_Search(mat,prm,prm_size);
    put("The permutation :"); put(prm);
    put(" with size : "); put(prm_size,1); new_line;
    if prm_size = mat'last(1) then
      put_line("The greedy search solved the pivot selection problem :");
      put(Apply_Permutation(mat,prm));
    else
      put_line("Running the Ford-Fulkerson algorithm ...");
      Pivot_Selection.maxBPM(mat,selected,selected_size);
      put("The size of the selection : "); put(selected_size,1); new_line;
      put("The selected columns :"); put(selected); new_line;
      if selected_size = natural32(mat'last(1)) then
        for i in prm'range loop
          prm(i) := natural32(selected(i));
        end loop;
        put_line("The solution applied to the matrix :");
        put(Apply_Permutation(mat,prm));
      end if;
     -- put_line("Running the Hopcroft-Karp algorithm ...");
     -- Bipartite_Match(mat,prm,prm_size);
    end if;
  end Test;

  procedure Random_Test ( nbrows,nbcols : in integer32 ) is

    use Standard_Random_Matrices;

    mat : Boolean_Matrices.Matrix(1..nbrows,1..nbcols);
    prb : double_float := 0.0;

  begin
    put("Give the probability for true : ");
    get(prb);
    mat := Random_Matrix(natural32(nbrows),natural32(nbcols),prb);
    Test(mat);
  end Random_Test;

  procedure Test_Input ( nbrows,nbcols : in integer32 ) is

    mat : Boolean_Matrices.Matrix(1..nbrows,1..nbcols);
    bit : natural32 := 0;

  begin
    for i in 1..nbrows loop
      put("Give "); put(nbcols,1);
      put(" bits for row "); put(i,1); put(" : ");
      for j in 1..nbcols loop
        get(bit);
        mat(i,j) := (bit = 1);
      end loop;
    end loop;
    Test(mat);
  end Test_Input;

  procedure Main is

    nbrows,nbcols : integer32 := 0;
    ans : character;

  begin
    put("Give the number of rows : "); get(nbrows);
    put("Give the number of columns : "); get(nbcols);
    put("Generate a random matrix ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Random_Test(nbrows,nbcols);
     else Test_Input(nbrows,nbcols);
    end if;
  end Main;

end Test_Pivot_Selection;
