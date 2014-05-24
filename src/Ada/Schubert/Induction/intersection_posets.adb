with unchecked_deallocation;
with text_io;                           use text_io;
with Checker_Moves;                     use Checker_Moves;
with Checker_Boards_io;                 use Checker_Boards_io;

package body Intersection_Posets is

  function Create ( ps : Poset ) return Poset_Node is

    res : Poset_Node;

  begin
    res.ps := ps;
    res.first_parent := null;   
    res.first_child := null;
    return res;
  end Create;

  function Create ( m : integer32; ps : Poset ) return Intersection_Poset is

    ips : Intersection_Poset(m);
    root : constant Link_to_Poset_Node := new Poset_Node'(Create(ps));

  begin
    ips.level := 1;
    Construct(root,ips.nodes(ips.level));
    ips.last(ips.level) := ips.nodes(ips.level);
    for k in 2..m loop
      ips.nodes(k) := Poset_List(Lists_of_Poset_Nodes.Null_List);
      ips.last(k) := ips.nodes(k);
    end loop;
    return ips;
  end Create;

  procedure Intersect ( ips : in out Intersection_Poset; 
                        pnd : in Link_to_Poset_Node; w : in Vector;
                        silent : in boolean ) is

    ps : constant Poset := pnd.ps;
    tmp : Link_to_Node := ps.white(ps.white'last);
    p : constant Vector := ps.black(ps.black'first).all;
    n : constant integer32 := p'last;
    lp1 : constant integer32 := ips.level + 1;
    isin : boolean;
    pnd_child : Link_to_Poset_node;
    m : Natural_Number;

  begin
    while tmp /= null loop
      if not silent
       then Write_Bracket(tmp.cols); put(" and "); Write_Bracket(w); 
      end if;
      if Happy_Checkers(ps.black(ps.black'first).all,tmp.cols,w) then
        Retrieve(ips.nodes(lp1),tmp.cols,w,isin,pnd_child);
        if isin then
          if not silent
           then put_line(" are happy and have already created children.");
          end if;
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
          if not silent
           then put_line(" are happy and will create children...");
          end if;
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
        if not silent
         then put_line(" are not happy and will not create any children.");
        end if;
      end if;
      tmp := tmp.next_sibling;
    end loop;
  end Intersect;

  procedure Intersect ( ips : in out Intersection_Poset; w : in Vector;
                        silent : in boolean ) is

    tmp : Poset_List := ips.nodes(ips.level);
    pnd : Link_to_Poset_Node;

  begin
    while not Is_Null(tmp) loop
      pnd := Head_Of(tmp);
      Intersect(ips,pnd,w,silent);
      tmp := Tail_Of(tmp);
    end loop;
    ips.level := ips.level + 1;
  end Intersect;

-- SELECTORS :

  function Degree_of_Freedom
             ( ips : Intersection_Poset; k : integer32 ) return natural32 is

    pnd : constant Link_to_Poset_Node := Head_Of(ips.nodes(k));

  begin
    return Degree_of_Freedom(pnd.ps);
  end Degree_of_Freedom;

  function Degree_of_Freedom ( ips : Intersection_Poset ) return natural32 is
  begin
    return Degree_of_Freedom(ips,ips.level);
  end Degree_of_Freedom;

  function Final_Sum ( ips : Intersection_Poset ) return Natural_Number is

    res : Natural_Number := create(natural32(0));
    pl : Poset_List;
    tmp : Poset_List;
    pnd : Link_to_Poset_Node;

  begin
    if ips.level > 0 then
      pl := ips.nodes(ips.level);
      tmp := pl;
      while not Is_Null(tmp) loop
        pnd := Head_Of(tmp);
        Add(res,pnd.ps.white(pnd.ps.white'last).coeff);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    return res;
  end Final_Sum;

  function Retrieve ( pl : Poset_List; k : integer32 )
                    return Link_to_Poset_Node is

    res : Link_to_Poset_Node := null;
    tmp : Poset_List := pl;

  begin
    for i in 1..(k-1) loop
      exit when Is_Null(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    if not Is_Null(tmp)
     then res := Head_Of(tmp);
    end if;
    return res;
  end Retrieve;

  procedure Retrieve ( pl : in Poset_List; rows,cols : in Vector;
                       isin : out boolean; pnd : out Link_to_Poset_Node ) is

    tmp : Poset_List := pl;
    lnd : Link_to_Node;

  begin
    isin := false;
    while not Is_Null(tmp) and not isin loop
      pnd := Head_Of(tmp);
      lnd := pnd.ps.white(pnd.ps.white'first);
      if Equal(lnd.rows,rows) and Equal(lnd.cols,cols)
       then isin := true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
  end Retrieve;

  function Is_Parent ( parent,child : Poset ) return boolean is

    child_rows : constant Vector := Root_Rows(child);
    ptr : Link_to_Node := parent.white(parent.white'last);

  begin
    while ptr /= null loop
      if Equal(ptr.cols,child_rows)
       then return true;
       else ptr := ptr.next_sibling;
      end if;
    end loop;
    return false;
  end Is_Parent;

  function Is_Parent ( pnd,cnd : Poset_Node ) return boolean is

    parent : constant Poset := pnd.ps;
    child : constant Vector := Root_Rows(cnd.ps);
    ptr : Link_to_Node := parent.white(parent.white'last);

  begin
    while ptr /= null loop
      if Equal(ptr.cols,child)
       then return true;
       else ptr := ptr.next_sibling;
      end if;
    end loop;
    return false;
  end Is_Parent;

  function Is_Child ( child,parent : Poset ) return boolean is
  begin
    return Is_Parent(parent,child);
  end Is_Child;

  function Is_Child ( cnd,pnd : Poset_Node ) return boolean is
  begin
    return Is_Parent(pnd,cnd);
  end Is_Child;

  procedure Enumerate_Parents ( pl : in Poset_List; nd : in Poset_Node ) is

    tmp : Poset_List := pl;
    parent : Link_to_Poset_Node;

  begin
    while not Is_Null(tmp) loop
      parent := Head_Of(tmp);
      if Is_Parent(parent.all,nd)
       then Process_Parent(parent);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Enumerate_Parents;

  function Number_of_Parents
             ( pl : Poset_List; nd : Poset_Node ) return natural32 is

    res : natural32 := 0;

    procedure Count_Parent ( pnd : in Link_to_Poset_Node ) is
    begin
      res := res + 1;
    end Count_Parent;
    procedure Count_Parents is new Enumerate_Parents(Count_Parent);

  begin
    Count_Parents(pl,nd);
    return res;
  end Number_of_Parents;

  function Retrieve_Parent 
             ( pl : Poset_List; nd : Poset_Node; k : integer32 )
             return Link_to_Poset_Node is

    res : Link_to_Poset_Node;
    cnt : integer32 := 0;

    procedure Get_Parent ( pnd : in Link_to_Poset_Node ) is
    begin
      cnt := cnt + 1;
      if cnt = k
       then res := pnd;
      end if;
    end Get_Parent;
    procedure Get_Parents is new Enumerate_Parents(Get_Parent);

  begin
    Get_Parents(pl,nd);
    return res;
  end Retrieve_Parent;

-- DESTRUCTORS :

  procedure Clear ( ps : in out Poset_Node ) is
  begin
    Clear(ps.ps);
  end Clear;

  procedure Clear ( ps : in out Link_to_Poset_Node ) is

    procedure free is
      new unchecked_deallocation(Poset_Node,Link_to_Poset_Node);

  begin
    if ps /= null then
      Clear(ps.all);
      free(ps);
    end if;
  end Clear;

  procedure Clear ( pl : in out Poset_List ) is

    tmp : Poset_List := pl;
    ps : Link_to_Poset_Node;

  begin
    while not Is_Null(tmp) loop
      ps := Head_Of(tmp);
      Clear(ps);
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Poset_Nodes.Clear(Lists_of_Poset_Nodes.List(pl));
  end Clear;

  procedure Clear ( apl : in out Array_of_Poset_Lists ) is
  begin
    for i in apl'range loop
      Clear(apl(i));
    end loop;
  end Clear;

  procedure Clear ( ps : in out Intersection_Poset ) is
  begin
    Clear(ps.nodes);
  end Clear;

end Intersection_Posets;
