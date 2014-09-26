with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Checker_Boards_io;                  use Checker_Boards_io;
with Checker_Moves;                      use Checker_Moves;
with Checker_Posets,Checker_Posets_io;   use Checker_Posets,Checker_Posets_io;
with Intersection_Posets;                use Intersection_Posets;
with Intersection_Posets_io;             use Intersection_Posets_io;

procedure ts_lrhom is

-- DESCRIPTION :
--   Interactive development of the application of Littlewood-Richardson
--   homotopies to solve general Schubert problems.

  function Read_Conditions
             ( n,k,m : in integer32 ) return Intersection_Poset is

    res : Intersection_Poset(m-1);
    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    rows,cols : Vector(1..k);
    ps : Poset;

  begin
    put_line("Reading the first two intersection conditions...");
    Read_Permutation(rows); Read_Permutation(cols);
    if not Happy_Checkers(p,rows,cols) then
      put_line("Your conditions form an unhappy configuration.");
    else
      ps := Create(n,rows,cols);
      res := Create(m-1,ps);
      for k in 3..m loop
        put("Reading intersection condition "); put(k,1); put_line("...");
        Read_Permutation(cols);
        Intersect(res,cols,false);
      end loop;
    end if;
    return res;
  end Read_Conditions;

  procedure Walk_Intersection_Poset ( n,k,m : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for m intersection conditions on k-planes in n-space,
  --   and writes the expansion to resolve the conditions from the top down.

    ips : Intersection_Poset(m-1) := Read_Conditions(n,k,m);
    tmp : Poset_List;
    lpn : Link_to_Poset_Node;

  begin
    put_line("Resolving the intersection conditions :");
    Write_Expansion(ips);
    for i in 1..ips.m loop
      put("The nodes at level "); put(i,1); put_line(" :");
      tmp := ips.nodes(i);
      for j in 1..Length_Of(ips.nodes(i)) loop
        lpn := Head_Of(tmp);
        put("-> poset node "); put(j,1); put_line(", root and leaves :");
       -- Checker_Posets_io.Write(lpn.ps);
        Checker_Posets_io.Write_Nodes_in_Poset(lpn.ps,lpn.ps.black'first);
        Checker_Posets_io.Write_Nodes_in_Poset(lpn.ps,lpn.ps.black'last);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Walk_Intersection_Poset;

  procedure Main is

  -- DESCRIPTION :
  --   Provides the user with a list of test procedures
  --   and prompts to make a choice.

    n,k,m : integer32 := 0;

  begin
    put_line("MENU to solve Schubert problems with LR homotopies : ");
    put_line("  0. walk through intersection poset from leaves to root.");
    new_line;
    put("Give the ambient dimension : "); get(n);
    put("Give the dimension of the planes : "); get(k);
    put("Give the number of conditions : "); get(m);
    Walk_Intersection_Poset(n,k,m);
  end Main;

begin
  Main;
end ts_lrhom;
