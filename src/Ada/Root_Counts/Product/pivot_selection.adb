package body Pivot_Selection is

  function Is_In ( vec : Standard_Natural_Vectors.Vector;
                   vec_end : integer32; nbr : in natural32 )
                 return boolean is
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
  begin
    prm := (prm'range => 0);
    size := 0;
    for row in mat'range(1) loop
      for col in mat'range(2) loop
        if mat(row,col) then
          if not Is_In(prm,prm'last,natural32(col)) then
            size := size + 1;
            prm(row) := natural32(col); exit;
          end if;
        end if;
      end loop;
    end loop;
  end Greedy_Search;

  procedure dfs4bpm
              ( graph : in Boolean_Matrices.Matrix;
                vtx : in integer32;
                seen : in out Boolean_Vectors.Vector;
                match : in out Standard_Integer_Vectors.Vector;
                found : out boolean ) is
  begin
    for col in graph'range(2) loop -- try every column one by one
      -- if there is an edge from col to vtx and col is not visited
      if col <= graph'last(1) then -- may have more columns than rows
        if(graph(col,vtx) and not seen(col)) then
          seen(col) := true; -- mark col as visited
          -- If col is not assigned to a row or previously assigned row
          -- for col (which is match(col)) has an alternate col available.
          -- Since col is marked as visited in the above line, match(col)
          -- in the following recursive call will not get col again.
          if(match(col) < 0) then
            match(col) := vtx;
            found := true; return;
          else
            dfs4bpm(graph,match(col),seen,match,found);
            if found
             then match(col) := vtx; return;
            end if;
          end if;
        end if;
      end if;
    end loop;
    found := false; 
  end dfs4bpm;

  procedure maxBPM
              ( graph : in Boolean_Matrices.Matrix;
                match : out Standard_Integer_Vectors.Vector;
                size : out natural32 ) is

    seen : Boolean_Vectors.Vector(graph'range(1));
    found : boolean;

  begin
    match := (match'range => -1); -- all columns are unassigned
    size := 0;
    for col in graph'range(2) loop
      seen := (seen'range => false);
      dfs4bpm(graph,col,seen,match,found); -- can col be assigned?
      if found
       then size := size + 1;
      end if;
    end loop;
  end maxBPM;

end Pivot_Selection;
