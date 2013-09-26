with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Brackets,Brackets_io;              use Brackets,Brackets_io;
with Checker_Boards_io;                 use Checker_Boards_io;
with Checker_Moves;                     use Checker_Moves;

package body Schubert_Posets is

  procedure Initialize_Leaves_Coefficients ( ps : in out Poset ) is

    nd : Link_to_Node := ps.white(ps.white'last);

  begin
    while nd /= null loop
      Clear(nd.coeff);
      nd.coeff := Create(integer(1));
      nd := nd.next_sibling;
    end loop;
  end Initialize_Leaves_Coefficients;

  procedure Initialize_Leaves_Coefficients
               ( ps : in out Poset; pl : in Poset_List ) is

    nd : Link_to_Node := ps.white(ps.white'last);
    tmp : Poset_List;
    pnd : Link_to_Poset_Node;

  begin
    while nd /= null loop
      Clear(nd.coeff);
      nd.coeff := Create(integer(0));
      tmp := pl;
      while not Is_Null(tmp) loop
        pnd := Head_Of(tmp);
        if Equal(nd.cols,Root_Rows(pnd.ps)) -- nd is parent of pnd.ps
         then Add(nd.coeff,pnd.ps.white(pnd.ps.white'first).coeff);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      nd := nd.next_sibling;
    end loop;
  end Initialize_Leaves_Coefficients;

  procedure Specialize_Leaves ( ps : in out Poset; w : in Vector ) is

    nd : Link_to_Node := ps.white(ps.white'last);

  begin
    while nd /= null loop
      Clear(nd.coeff);
      Write_Bracket(nd.cols); put(" and "); Write_Bracket(w);
      if Happy_Checkers(ps.black(ps.black'first).all,nd.cols,w) then
        put_line(" are happy and contribute one solution");
        nd.coeff := Create(integer(1));
      else
        put_line(" are not happy and do not contribute a solution");
        nd.coeff := Create(integer(0));
      end if;
      nd := nd.next_sibling;
    end loop;
  end Specialize_Leaves;

  procedure Count_Roots ( ps : in out Poset ) is

    nd : Link_to_Node;

  begin
    for i in reverse ps.white'first..ps.white'last-1 loop
      nd := ps.white(i);
      while nd /= null loop
        Clear(nd.coeff); nd.coeff := Create(integer(0));
        if nd.stay_child /= null
         then Add(nd.coeff,nd.stay_child.coeff);
        end if;
        if nd.swap_child /= null
         then Add(nd.coeff,nd.swap_child.coeff);
        end if;
        nd := nd.next_sibling;
      end loop;
    end loop;
  end Count_Roots;

  procedure Specialize ( ips : in out Intersection_Poset;
                         pnd : in Link_to_Poset_Node;
                         w : in Vector; fail : out boolean ) is

    ps : constant Poset := pnd.ps;
    tmp : Link_to_Node := ps.white(ps.white'last);
    p : constant Vector := ps.black(ps.black'first).all;
    n : constant integer32 := p'last;
    lp1 : constant integer32 := ips.level + 1;
    isin : boolean;
    pnd_child : Link_to_Poset_node;
    m : Natural_Number;

  begin
    put("At current level "); put(ips.level,1); put_line(" :");
    fail := true;
    while tmp /= null loop
      Write_Bracket(tmp.cols); put(" and "); Write_Bracket(w); 
      if Happy_Checkers(ps.black(ps.black'first).all,tmp.cols,w) then
        fail := false;
        Retrieve(ips.nodes(lp1),tmp.cols,w,isin,pnd_child);
        if isin then
          put_line(" are happy and have already created children.");
          if pnd.first_child = null
           then pnd.first_child := pnd_child;
          end if;
          m := tmp.coeff;
         -- put("multiplicity at the current poset : "); put(m,1); new_line;
         -- put("multiplicity at the child : ");
         -- put(pnd_child.ps.white(pnd_child.ps.white'first).coeff,1);
         -- new_line;
          Add_Multiplicity(pnd_child.ps,m);
        else
          put_line(" are happy and will create children...");
          declare
            child : constant Poset := Create(n,tmp.coeff,tmp.cols,w);
            child_node : constant Link_to_Poset_Node
                       := new Poset_Node'(Create(child));
          begin
           -- put_line("The child poset : "); Write(child);
            child_node.first_parent := pnd;
            if pnd.first_child = null
             then pnd.first_child := child_node;
            end if;
            Append(ips.nodes(lp1),ips.last(lp1),child_node);
          end;
        end if;
      else
        put_line(" are not happy and will not create any children.");
        Clear(tmp.coeff); tmp.coeff := Create(integer(0)); -- reset coefficient
      end if;
      tmp := tmp.next_sibling;
    end loop;
  end Specialize;

  procedure Specialize ( ips : in out Intersection_Poset;
                         w : in Vector; fail : out boolean ) is

    tmp : Poset_List := ips.nodes(ips.level);
    pnd : Link_to_Poset_Node;
    sub_fail : boolean := true;

  begin
    fail := true;
    while not Is_Null(tmp) loop
      pnd := Head_Of(tmp);
      Specialize(ips,pnd,w,sub_fail);
      if not sub_fail
       then fail := false;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    if not fail 
     then ips.level := ips.level + 1;
    end if;
  end Specialize;

  function Specialize ( n,k : integer32; b : Bracket_Monomial )
                      return Intersection_Poset is

    nb : constant integer32 := integer32(Number_of_Brackets(b));
    cd : constant Array_of_Brackets(1..nb) := Create(b);
    p : constant Standard_Natural_Vectors.Vector(1..n)
      := Identity_Permutation(natural32(n));
    r : constant Standard_Natural_Vectors.Vector(1..k)
      := Standard_Natural_Vectors.Vector(cd(1).all);
    c : constant Standard_Natural_Vectors.Vector(1..k)
      := Standard_Natural_Vectors.Vector(cd(2).all);
    ps : Poset;
    res : Intersection_Poset(nb-1);
    w : Standard_Natural_Vectors.Vector(1..k);
    fail : boolean;

  begin
    if not Happy_Checkers(p,r,c) then
      put(cd(1).all); put(" and "); put(cd(2).all);
      put_line(" are not happy, no solutions.");
      res.level := 0;
    else
      ps := Create(n,r,c);
      w := Standard_Natural_Vectors.Vector(cd(3).all);
      res := Create(nb-1,ps);
      fail := true;
      for i in 3..nb loop
        w := Standard_Natural_Vectors.Vector(cd(i).all);
        Specialize(res,w,fail);
        exit when fail;
      end loop;
    end if;
    return res;
  end Specialize;

  function Root_Count_at_Leaves
              ( ips : Intersection_Poset ) return Natural_Number is

    res : Natural_Number := Create(integer(0));
    tmp : Poset_List;
    pnd : Link_to_Poset_Node;
 
  begin
    if ips.level = ips.m then
      tmp := ips.nodes(ips.level);
      while not Is_Null(tmp) loop
        pnd := Head_Of(tmp);
        declare
          ps : constant Poset := pnd.ps;
          nd : constant Link_to_Node := ps.white(ps.white'last);
        begin
          Add(res,nd.coeff);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    return res;
  end Root_Count_at_Leaves;

  procedure Count_Roots ( ips : in Intersection_Poset ) is

    tmp : Poset_List;
    pnd : Link_to_Poset_Node;

  begin
    if ips.level = ips.m then
      tmp := ips.nodes(ips.level);
      while not Is_Null(tmp) loop
        pnd := Head_Of(tmp);
        Initialize_Leaves_Coefficients(pnd.ps);
        Count_Roots(pnd.ps);
        tmp := Tail_Of(tmp);
      end loop;
      for k in reverse ips.nodes'first..ips.level-1 loop
        tmp := ips.nodes(k);
        while not Is_Null(tmp) loop
          pnd := Head_Of(tmp);
          Initialize_Leaves_Coefficients(pnd.ps,ips.nodes(k+1));
          Count_Roots(pnd.ps);
          tmp := Tail_Of(tmp);
        end loop;
      end loop;
    end if;
  end Count_Roots;

end Schubert_Posets;
