--with unchecked_deallocation;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Brackets;                           use Brackets;
with Brackets_io;                        use Brackets_io;

package body Pieri_Root_Counts is

-- procedure free is new unchecked_deallocation(Nodal_Pair,Link_to_Nodal_Pair);

  type Boolean_Array is array ( integer32 range <> ) of boolean;

  function Create ( n,d : integer32; t1,t2 : Pieri_Tree )
                  return List_of_Paired_Nodes is

    res,res_last : List_of_Paired_Nodes;
    h1 : constant natural32 := Height(t1);
    h2 : constant natural32 := Height(t2);
    b1,b2 : Bracket(1..d);
    firstlnd : Link_to_Pieri_Node;

    procedure Check_Pair ( lnd : in Link_to_Pieri_Node;
                           continue : out boolean ) is
    begin
      b2 := lnd.node;
      if Pieri_Condition(natural32(n),b1,b2) then
        declare
          lpnd : Paired_Nodes;
        begin
          lpnd.left := firstlnd;
          lpnd.right := lnd;
          Append(res,res_last,lpnd);
        end;
      end if;
      continue := true;
    end Check_Pair;
    procedure Check_Pairs is new Enumerate_Nodes(Check_Pair);

    procedure Count_First ( lnd : in Link_to_Pieri_Node;
                            continue : out boolean ) is
    begin
      b1 := lnd.node;
      firstlnd := lnd;
      Check_Pairs(t2,h2);
      continue := true;
    end Count_First;
    procedure First_Leaves is new Enumerate_Nodes(Count_First);

  begin
    First_Leaves(t1,h1);
    return res;
  end Create;

  function Create ( pnd : Paired_Nodes ) return Paired_Chain is

    res : Paired_Chain(1..integer32(Height(pnd)));
    ind : integer32 := res'last;

  begin
    res(ind) := pnd;
    while not At_First_Branch_Point(res(ind)) loop         -- fill in
      ind := ind - 1;
      res(ind) := Ancestor(res(ind+1));
    end loop;
    if ind = 1 then
      return res;
    else
      for i in 1..res'last-ind+1 loop                  -- shift down
        res(i) := res(i+ind-1);
      end loop;
      return res(1..res'last-ind+1);
    end if;
  end Create;

  procedure Connect ( ancnp,np : in Link_to_Nodal_Pair ) is

  -- DESCRIPTION :
  --   Connects the ancestor paired nodes with the paired nodes np.

    ancpnd : constant Paired_Nodes := Ancestor(np.pnd);
    j1 : constant integer32 := Jump(ancpnd.left.node,np.pnd.left.node);
    j2 : constant integer32 := Jump(ancpnd.right.node,np.pnd.right.node);

  begin
    ancnp.pnd := ancpnd;
    ancnp.children(j1,j2) := np;
    np.ancestor := ancnp;
  end Connect;

  procedure Initial_Branch ( root : in out Link_to_Nodal_Pair;
                             np : in Link_to_Nodal_Pair ) is

  -- DESCRIPTION :
  --   Constructs the initial branch in the tree of paired nodes.

  begin
    if At_First_Branch_Point(np.pnd) then
      root := np;
    else
      declare
        acc : constant Link_to_Nodal_Pair := new Nodal_Pair(np.d);
      begin
        acc.sols := 1;
        Connect(acc,np);
        Initial_Branch(root,acc);
      end;
    end if;
  end Initial_Branch;

  procedure Merge ( root : in Nodal_Pair;
                    current : in Link_to_Nodal_Pair; k : in integer32;
                    chain : in Paired_Chain ) is

  -- DESCRIPTION :
  --   Merges the chain with the root of the tree, at level k.

    j1,j2 : integer32;

  begin
    j1 := Jump(chain(k).left.node,chain(k+1).left.node);
    j2 := Jump(chain(k).right.node,chain(k+1).right.node); 
    if current.children(j1,j2) = null then
      declare
        newnp : constant Link_to_Nodal_Pair := new Nodal_Pair(current.d);
      begin
        newnp.pnd := chain(k+1);
        if Is_In(root,newnp.pnd)
         then newnp.sols := 0;
         else newnp.sols := 1;
        end if;
        current.children(j1,j2) := newnp;
        newnp.ancestor := current;
      end;
    else
      if current.children(j1,j2).sols > 0 then
        current.children(j1,j2).sols
          := current.children(j1,j2).sols + 1;
      end if;
    end if;
    if k+1 < chain'last
     then Merge(root,current.children(j1,j2),k+1,chain);
    end if;
  end Merge;

  function Create ( d : integer32; lp : List_of_Paired_Nodes )
                  return Nodal_Pair is

    root : Nodal_Pair(d);
    lroot : Link_to_Nodal_Pair := new Nodal_Pair'(root);
    first : constant Link_to_Nodal_Pair := new Nodal_Pair(d);
    tmp : List_of_Paired_Nodes := Tail_Of(lp);

  begin
    first.pnd := Head_Of(lp);
    first.sols := 1;
    lroot.sols := 1;
    Initial_Branch(lroot,first);
    while not Is_Null(tmp) loop
      declare
        pnd : constant Paired_Nodes := Head_Of(tmp);
        chn : constant Paired_Chain := Create(pnd);
      begin
        lroot.sols := lroot.sols + 1;
        Merge(lroot.all,lroot,1,chn);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return lroot.all;
  end Create;

-- SELECTORS :

  function Height ( pnd : Paired_Nodes ) return natural32 is
  begin
    if pnd.left.h >= pnd.right.h
     then return pnd.left.h;
     else return pnd.right.h;
    end if;
  end Height;

  function Equal ( pnd1,pnd2 : Paired_Nodes ) return boolean is
  begin
    return (Is_Equal(pnd1.left.node,pnd2.left.node)
        and Is_Equal(pnd1.right.node,pnd2.right.node));
  end Equal;

  function At_First_Branch_Point ( pnd : Paired_Nodes ) return boolean is
  begin
    if pnd.left.h /= pnd.right.h then
      return false;
    elsif ((pnd.left.c > 1) or (pnd.right.c > 1)) then
      return false;
    else return (((pnd.left.i = 0) and (pnd.left.c = 1))
           or else ((pnd.right.i = 0) and (pnd.right.c = 1)));
    end if;
  end At_First_Branch_Point;

  function At_Leaves ( pnd : Paired_Nodes ) return boolean is
  begin
    return (Is_Leaf(pnd.left.all) and Is_Leaf(pnd.right.all));
  end At_Leaves;

  function Ancestor ( pnd : Paired_Nodes ) return Paired_Nodes is

    res : Paired_Nodes;

  begin
    if pnd.left.h = pnd.right.h then
      res.left := pnd.left.ancestor;
      res.right := pnd.right.ancestor;
    elsif pnd.left.h > pnd.right.h then
      res.left := pnd.left.ancestor;
      res.right := pnd.right;
    else
      res.left := pnd.left;
      res.right := pnd.right.ancestor;
    end if;
    return res;
  end Ancestor;

  function First_Branch_Point ( pnd : Paired_Nodes ) return Paired_Nodes is
  begin
    if At_First_Branch_Point(pnd)
     then return pnd;
     else return First_Branch_Point(Ancestor(pnd));
    end if;
  end First_Branch_Point;

  function Height ( np : Nodal_Pair ) return natural32 is
  begin
    if np.pnd.left.h >= np.pnd.right.h
     then return np.pnd.left.h;
     else return np.pnd.right.h;
    end if;
  end Height;

  function Is_In ( root : Nodal_Pair; pnd : Paired_Nodes ) return boolean is
  begin
    if Equal(root.pnd,pnd) then
      return true;
    else
      for j1 in root.children'range(1) loop
        for j2 in root.children'range(2) loop
          if root.children(j1,j2) /= null then
            if Is_In(root.children(j1,j2).all,pnd)
             then return true;
            end if;
          end if;
        end loop;
      end loop;
    end if;
    return false;
  end Is_In;

  function Number_of_Paths ( root : Nodal_Pair ) return natural32 is

    res : natural32 := root.sols;

  begin
    for j1 in root.children'range(1) loop
      for j2 in root.children'range(2) loop
        if root.children(j1,j2) /= null then
          if not At_Leaves(root.children(j1,j2).pnd)
           then res := res + Number_of_Paths(root.children(j1,j2).all);
          end if;
        end if;
      end loop;
    end loop;
    return res;
  end Number_of_Paths;

-- FORMATTED OUTPUT :

  procedure Write ( file : in file_type; chn : in Paired_Chain ) is
  begin
    for i in chn'first..(chn'last-1) loop
      put(file,"("); put(file,chn(i).left.node);
      put(file,","); put(file,chn(i).right.node); put(file,") < ");
    end loop;
    put(file,"("); put(file,chn(chn'last).left.node);
    put(file,","); put(file,chn(chn'last).right.node); put_line(file,")");
  end Write;

  function Last_Child ( np : Nodal_Pair; i,j : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the (i,j)th child is the last child of the node.

  begin
    for j1 in j+1..np.children'last(2) loop
      if np.children(i,j1) /= null
       then return false;
      end if;
    end loop;
    for i1 in i+1..np.children'last(1) loop
      for j1 in np.children'range(2) loop
        if np.children(i1,j1) /= null
         then return false;
        end if;
      end loop;
    end loop;
    return true;
  end Last_Child;

  procedure Write_Labels ( file : in file_type; np : in Nodal_Pair;
                           j1,j2,h : in integer32; last : in Boolean_Array ) is

  -- DESCRIPTION :
  --   Writes the contents of the nodal pair with the jumps, taking into
  --   account which children appeared last.
  --   The current node is at height h in the nodal pair tree.

   -- first : Paired_Nodes := First_Branch_Point(np.pnd);

  begin
    if h /= 0
     then put(file,"   ");
    end if;
    for i in 1..h-1 loop
      if last(i)
       then put(file,"     ");
       else put(file,"|    ");       
      end if;
    end loop;
    if h /= 0 then
      put(file,"!-+(");
      put(file,j1,1); put(file,","); put(file,j2,1);
      put(file,")");
    end if;
    put(file,"("); put(file,np.pnd.left.node);
    put(file,","); put(file,np.pnd.right.node);
    put(file,") ");
    put(file,np.sols,1);
    new_line(file);
  end Write_Labels;

  procedure Write_Nodes ( file : in file_type; np : in Nodal_Pair;
                          j1,j2,h : in integer32;
                          last : in out Boolean_Array ) is

  -- DESCRIPTION :
  --   Writes the contents of the nodal pair, followed by the children.

  begin
    Write_Labels(file,np,j1,j2,h,last);
    for jj1 in np.children'range(1) loop
      for jj2 in np.children'range(2) loop
        if np.children(jj1,jj2) /= null then
          last(h+1) := Last_Child(np,jj1,jj2);
          Write_Nodes(file,np.children(jj1,jj2).all,jj1,jj2,h+1,last);
        end if;
      end loop;
    end loop;
  end Write_Nodes;

  procedure Write ( file : in file_type; root : in Nodal_Pair ) is

    last : Boolean_Array(1..integer32(Height(root))+1);

  begin
    Write_Nodes(file,root,1,1,0,last);
  end Write;

end Pieri_Root_Counts;
