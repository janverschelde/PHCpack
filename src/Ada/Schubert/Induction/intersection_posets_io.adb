with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Checker_Boards_io;                 use Checker_Boards_io;
with Checker_Posets;                    use Checker_Posets;
with Checker_Posets_io;                 use Checker_Posets_io;

package body Intersection_Posets_io is

  procedure Write_Parents ( pl : in Poset_List; nd : in Poset_Node ) is

    cnt : integer32 := 0;

    procedure Write_Parent ( pnd : in Link_to_Poset_Node ) is

      ps : constant Poset := pnd.ps;
      leaf : Link_to_Node := ps.white(ps.white'last);

    begin
      cnt := cnt + 1;
      put("parent #"); put(cnt,1); put(" : ");
      Write_Bracket(Root_Rows(ps));
      Write_Bracket(Root_Columns(ps)); 
      put(" -> ");
      while leaf /= null loop
        if Equal(Root_Rows(nd.ps),leaf.cols) then
          Write_Bracket(leaf.rows);
          Write_Bracket(leaf.cols); new_line;
        end if;
        leaf := leaf.next_sibling;
      end loop;
    end Write_Parent;
    procedure Write is
      new Intersection_Posets.Enumerate_Parents(Write_Parent);

  begin
    Write(pl,nd);
  end Write_Parents;

  procedure Write_Formal_Equations
               ( ips : in Intersection_Poset; k : in integer32 ) is

    tmp : Poset_List := ips.nodes(k);
    pnd : Link_to_Poset_Node;

  begin
    while not Is_Null(tmp) loop
      pnd := Head_Of(tmp);
      Write_Formal_Equation(pnd.ps);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Formal_Equations;

  procedure Write_Formal_Equations
               ( file : in file_type;
                 ips : in Intersection_Poset; k : in integer32 ) is

    tmp : Poset_List := ips.nodes(k);
    pnd : Link_to_Poset_Node;

  begin
    while not Is_Null(tmp) loop
      pnd := Head_Of(tmp);
      Write_Formal_Equation(file,pnd.ps);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Formal_Equations;

  procedure Write_Formal_Equations ( ips : in Intersection_Poset ) is
  begin
    for k in 1..ips.level loop
      Write_Formal_Equations(ips,k);
    end loop;
  end Write_Formal_Equations;

  procedure Write_Formal_Equations
              ( file : in file_type;ips : in Intersection_Poset ) is
  begin
    for k in 1..ips.level loop
      Write_Formal_Equations(file,ips,k);
    end loop;
  end Write_Formal_Equations;

  procedure Write_Lefthand_Product
              ( ips : in Intersection_Poset; k : in integer32 ) is

    pnd : Link_to_Poset_Node := Head_Of(ips.nodes(1));
    ps : Poset;
    kk : integer32;

  begin
    if ips.level >= k then
      if k = 1 then
        Write_Formal_Product(pnd.ps);
        kk := k + 1;
      else
        kk := k;
      end if;
      for i in kk..ips.level loop
        exit when Is_Null(ips.nodes(i));
        pnd := Head_Of(ips.nodes(i));
        ps := pnd.ps;
        put("*");
        Write_Bracket(ps.white(ps.white'first).cols);
      end loop;
    end if;
  end Write_Lefthand_Product;

  procedure Write_Lefthand_Product
              ( file : in file_type;
                ips : in Intersection_Poset; k : in integer32 ) is

    pnd : Link_to_Poset_Node := Head_Of(ips.nodes(1));
    ps : Poset;
    kk : integer32;

  begin
    if ips.level >= k then
      if k = 1 then
        Write_Formal_Product(file,pnd.ps);
        kk := k + 1;
      else
        kk := k;
      end if;
      for i in kk..ips.level loop
        exit when Is_Null(ips.nodes(i));
        pnd := Head_Of(ips.nodes(i));
        ps := pnd.ps;
        put(file,"*");
        Write_Bracket(file,ps.white(ps.white'first).cols);
      end loop;
    end if;
  end Write_Lefthand_Product;

  procedure Write_Lefthand_Product ( ips : in Intersection_Poset ) is
  begin
    Write_Lefthand_Product(ips,1);
  end Write_Lefthand_Product;

  procedure Write_Lefthand_Product
              ( file : in file_type; ips : in Intersection_Poset ) is
  begin
    Write_Lefthand_Product(file,ips,1);
  end Write_Lefthand_Product;

  procedure Write_Final_Sum ( pl : in Poset_List ) is

    tmp : Poset_List := pl;
    pnd : Link_to_Poset_Node;

  begin
    while not Is_Null(tmp) loop
      pnd := Head_Of(tmp);
      Write_Final_Sum(pnd.ps);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Final_Sum;

  procedure Write_Final_Sum ( file : in file_type; pl : in Poset_List ) is

    tmp : Poset_List := pl;
    pnd : Link_to_Poset_Node;

  begin
    while not Is_Null(tmp) loop
      pnd := Head_Of(tmp);
      Write_Final_Sum(file,pnd.ps);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Final_Sum;

  procedure Write_Expansion ( ips : in Intersection_Poset ) is
  begin
    Write_Lefthand_Product(ips);
    new_line;
    for k in 2..ips.level loop
      exit when Is_Null(ips.nodes(k-1));
      put(" = (");
      Write_Final_Sum(ips.nodes(k-1));
      put(")");
      Write_Lefthand_Product(ips,k);
      new_line;
    end loop;
    if Is_Null(ips.nodes(ips.level)) then
      put_line(" = 0");
    else
      put(" = ");
      Write_Final_Sum(ips.nodes(ips.level));
      new_line;
    end if;
  end Write_Expansion;

  procedure Write_Expansion
              ( file : in file_type; ips : in Intersection_Poset ) is
  begin
    Write_Lefthand_Product(file,ips);
    new_line(file);
    for k in 2..ips.level loop
      exit when Is_Null(ips.nodes(k-1));
      put(file," = (");
      Write_Final_Sum(file,ips.nodes(k-1));
      put(file,")");
      Write_Lefthand_Product(file,ips,k);
      new_line(file);
    end loop;
    if Is_Null(ips.nodes(ips.level)) then
      put_line(file," = 0");
    else
      put(file," = ");
      Write_Final_Sum(file,ips.nodes(ips.level));
      new_line(file);
    end if;
  end Write_Expansion;

end Intersection_Posets_io;
