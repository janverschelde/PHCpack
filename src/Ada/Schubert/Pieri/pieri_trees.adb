with unchecked_deallocation;

package body Pieri_Trees is

-- UTILITIES FOR CREATION OF Pieri Trees :

  function Index_of_Increase ( nd : Pieri_Node ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index of increase between the current node nd and the
  --   ancestor node.  If the current node is the root, then the index
  --   of increase equals zero.

    bnd : Link_to_Pieri_Node;

  begin
    if nd.ancestor = null then
      return 0;
    else
      bnd := nd.ancestor;
      for i in nd.node'range loop
        if bnd.node(i) = nd.node(i)-1
         then return i;
        end if;
      end loop;
      return 0;
    end if;
  end Index_of_increase;

  function Branching_Level ( l : integer32; r : Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the current level l is a level where decreasing
  --   is allowed.

    bl : integer32 := 1;

  begin
    for i in r'first..r'last-1 loop
      bl := bl + integer32(r(i));
      if bl = l then
        return true;
      elsif bl > l then
        return false;
      end if;
    end loop;
    return false;
  end Branching_Level;

  procedure Create_Next ( n,d,l,h : in integer32; r : in Vector;
                          nd : in out Link_to_Pieri_Node ) is

  -- DESCRIPTION :
  --   Creates next level of nodes in the Pieri Tree.

  -- ON ENTRY :
  --   n         maximal entry in a bracket, dimension of whole space;
  --   d         number of entries in bracket;
  --   l         current level, must be strictly lower than h;
  --   h         height of the Pieri tree;
  --   nd        current node.

  -- ON RETURN :
  --   nd        node with updated links.

    indinc : constant integer32 := Index_of_Increase(nd.all);

  begin
    if Branching_Level(l,r) then              -- test if jumping-branching node
      nd.i := 0;
      nd.c := nd.ancestor.c + 1;
    end if;
    if integer32(nd.node(d)) <  n then        -- create right node
      declare
        rnd : Pieri_Node(d);
        lnd : Link_to_Pieri_Node;
      begin
        rnd.node := nd.node;                     -- adjust entries
        rnd.node(d) := rnd.node(d)+1;
        rnd.c := nd.c;
        rnd.i := nd.i + 1;
        rnd.h := nd.h + 1;
        lnd := new Pieri_Node'(rnd);
        lnd.ancestor := nd;                      -- establish connections
        nd.children(d) := lnd;
        if l < h                                 -- go to next level
         then Create_Next(n,d,l+1,h,r,lnd);
        end if;
      end;
    end if;
    for i in nd.node'first..(nd.node'last-1) loop
      if nd.node(i) < nd.node(i+1) - 1 then
        if ((i >= indinc) or else ((nd.i = 0) and (nd.c > 0))) then
          declare                         -- jumping-branching, create node
            rnd : Pieri_Node(d);
            lnd : Link_to_Pieri_Node;
          begin
            rnd.node := nd.node;                     -- adjust entries
            rnd.node(i) := rnd.node(i)+1;
            rnd.c := nd.c;
            rnd.i := nd.i + 1;
            rnd.h := nd.h + 1;
            lnd := new Pieri_Node'(rnd);
            lnd.ancestor := nd;               -- establish connections
            nd.children(i) := lnd;
            if l < h                               -- go to next level
             then Create_Next(n,d,l+1,h,r,lnd);
            end if;
          end;
        end if;
      end if;
    end loop;
  end Create_Next;

-- CREATOR :
  
  function Create ( n,d : natural32; r : Vector ) return Pieri_Tree is

    res : Pieri_Tree(integer32(d),r'last);
    hei : natural32;
    pnd : Pieri_Node(integer32(d));

  begin
    res.branches := r;
    for i in pnd.node'range loop                -- root node = [1 2 .. d]
      pnd.node(i) := natural32(i);
    end loop;
    pnd.c := 0;
    pnd.i := 0;
    pnd.h := 0;
    res.root := new Pieri_Node'(pnd);
    res.root.ancestor := null;
    hei := Height(res);
    if hei > 0                           -- create children
     then Create_Next(integer32(n),integer32(d),1,integer32(hei),r,res.root); 
    end if;
    return res;
  end Create;

-- SELECTORS :

  function Height ( t : Pieri_Tree ) return natural32 is

    res : natural32 := 0;

  begin
    for i in t.branches'range loop
      res := res + t.branches(i);
    end loop;
    return res;
  end Height;

  function Is_Leaf ( nd : Pieri_Node ) return boolean is
  begin
    for i in nd.children'range loop
      if nd.children(i) /= null
       then return false;
      end if;
    end loop;
    return true;
  end Is_Leaf;

  function Jump ( b1,b2 : Bracket ) return integer32 is
  begin
    for i in reverse b1'range loop
      if b1(i) < b2(i)
       then return i;
      end if;
    end loop;
    return 0;
  end Jump;

  function Jump ( nd : Pieri_Node ) return integer32 is
  begin
    if nd.ancestor = null
     then return 0;
     else return Jump(nd.ancestor.node,nd.node);
    end if;
  end Jump;

  function Lower_Jump_Decrease ( nd : Pieri_Node ) return Bracket is
  begin
    if ((nd.i = 0) or else (nd.c = 0)) then
      return nd.node;
    elsif nd.ancestor /= null then
      return Lower_Jump_Decrease(nd.ancestor.all);
    else
      return nd.node;
    end if;
  end Lower_Jump_Decrease;

  function Lowest_Jump_Decrease ( nd : Pieri_Node ) return Bracket is
  begin
    if (nd.c = 0) or ((nd.i = 0) and (nd.c = 1)) then
      return nd.node;
    elsif nd.ancestor /= null then
      return Lowest_Jump_Decrease(nd.ancestor.all);
    else
      return nd.node;
    end if;
  end Lowest_Jump_Decrease;

  function Upper_Jump_Decrease ( nd : Pieri_Node ) return Bracket is
  begin
    if ((nd.i = 0) or else (nd.c = 0)) then
      return nd.node;
    elsif nd.children(nd.node'last) /= null then
      return Upper_Jump_Decrease(nd.children(nd.node'last).all);
    else
      return nd.node;
    end if;
  end Upper_Jump_Decrease;

  procedure Enumerate_Nodes ( t : in Pieri_Tree; level : in natural32 ) is

    continue : boolean := true;

    procedure Visit_Nodes ( nd : in Link_to_Pieri_Node ) is
    begin
      if nd.h = level then
        Visit_Node(nd,continue);
      else
        for i in nd.children'range loop
          if nd.children(i) /= null
            then Visit_Nodes(nd.children(i));
          end if;
          exit when not continue;
        end loop;
      end if;
    end Visit_Nodes;

  begin
    if t.root /= null
     then Visit_Nodes(t.root);
    end if;
  end Enumerate_Nodes;

  procedure Enumerate_Chains ( t : in Pieri_Tree ) is

    b : Bracket_Array(1..integer32(Height(t)));
    continue : boolean := true;

    procedure Visit_Nodes ( nd : in Pieri_Node; ind : in integer32 ) is
    begin
      b(ind) := new Bracket'(nd.node);
      if ind = b'last then
        Visit_Chain(b,continue);
      else
        for i in nd.children'range loop
          if nd.children(i) /= null
           then Visit_Nodes(nd.children(i).all,ind+1);
          end if;
          exit when not continue;
        end loop;
      end if;
    end Visit_Nodes;

  begin
    if t.root /= null
     then Visit_Nodes(t.root.all,1);
    end if;
  end Enumerate_Chains;

  procedure Enumerate_Paired_Chains ( t1,t2 : in Pieri_Tree ) is

    continue : boolean := true;

    procedure Outer_Chain ( ob : in Bracket_Array; cont : out boolean ) is

      procedure Inner_Chain ( ib : in Bracket_Array; cont : out boolean ) is 
      begin
        Visit_Paired_Chain(ob,ib,continue);
        cont := continue;
      end Inner_Chain;
      procedure Inner_Chains is new Enumerate_Chains(Inner_Chain);

    begin
      Inner_Chains(t2);
      cont := continue;
    end Outer_Chain;
    procedure Outer_Chains is new Enumerate_Chains(Outer_Chain);

  begin
    Outer_Chains(t1);
  end Enumerate_Paired_Chains;

  function Pieri_Condition
             ( n : natural32; b1,b2 : Bracket ) return boolean is
  begin
    for i in b2'range loop
      if b2(i) > n+1 - b1(b1'last+1-i)       -- negation of weak inequality
       then return false;
      end if;
    end loop;
    for i in b1'first..b1'last-1 loop
      if n+1-b1(b1'last+1-i) >= b2(i+1)      -- negation of strong inequality
       then return false;
      end if;
    end loop;
    return true;
  end Pieri_Condition;

-- DESTRUCTOR :

  procedure Clear ( nd : in out Link_to_Pieri_Node ) is

    procedure free is new unchecked_deallocation(Pieri_Node,Link_to_Pieri_Node);

  begin
    if nd /= null
     then free(nd);
    end if;
  end Clear;

  procedure Clear_Children ( nd : in out Link_to_Pieri_Node ) is

  -- DESCRIPTION :
  --   Deallocation of the memory of all the children, before the memory
  --   occupied by the current node nd is released.  Applied recursively.

  begin
    for i in nd.children'range loop
      if nd.children(i) /= null
       then Clear_Children(nd.children(i));
      end if;
    end loop;
    Clear(nd);
  end Clear_Children;

  procedure Clear ( t : in out Pieri_Tree ) is
  begin
    Clear_Children(t.root);
  end Clear;

end Pieri_Trees;
