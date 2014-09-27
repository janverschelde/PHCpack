with unchecked_deallocation;
with Checker_Moves;                     use Checker_Moves;
with Checker_Localization_Patterns;     use Checker_Localization_Patterns;

package body Checker_Posets is

-- CREATORS :

  function Specializing_Moves ( n : integer32 ) return VecVec is

    res : VecVec(1..integer32(Number_of_Moves(natural32(n))));
    mvp : Vector(1..n) := Identity_Permutation(natural32(n));
    d,r : integer32;

  begin
    res(1) := new Vector'(mvp);
    for i in 2..res'last loop
      d := Descending_Checker(mvp);
      r := Rising_Checker(mvp,d);
      Move(mvp,d,r);
      res(i) := new Vector'(mvp);
    end loop;
    return res;
  end Specializing_Moves;

  function Generalizing_Moves ( n : integer32 ) return VecVec is

    res : VecVec(1..integer32(Number_of_Moves(natural32(n))));
    mvp : Vector(1..n) := Reverse_Permutation(natural32(n));
    f,a : integer32;

  begin
    res(1) := new Vector'(mvp);
    for i in 2..res'last loop
      f := Falling_Checker(mvp);
      a := Ascending_Checker(mvp,f);
      Move(mvp,a,f);
      res(i) := new Vector'(mvp);
    end loop;
    return res;
  end Generalizing_Moves;

  function Create ( r,c : Vector ) return Node is

    res : Node(r'last-r'first+1);

  begin
    res.coeff := Create(integer(1));
    res.rows := r;
    res.cols := c;
    res.first_parent := null;
    res.stay_child := null;
    res.swap_child := null;
    res.next_sibling := null;
    return res;
  end Create;

  function Copy_Child ( nd : Node ) return Node is

  -- DESCRIPTION :
  --   Returns a copy from the node (except for the links!).

    res : Node(nd.k);

  begin
    Copy(nd.coeff,res.coeff);
    res.rows := nd.rows;
    res.cols := nd.cols;
    res.first_parent := null;
    res.stay_child := null;
    res.swap_child := null;
    res.next_sibling := null;
    return res;
  end Copy_Child;

  function Stay_Child ( nd : Node; np : Vector ) return Node is

  -- DESCRIPTION :
  --   Returns a copy from the node (except for the links!)
  --   after making it happy according to the next location
  --   of the black checkers in np.

    res : Node(nd.k) := Copy_Child(nd);

  begin
    Make_Happy(np,res.rows,res.cols);
    return res;
  end Stay_Child;

  function Swap_Child ( nd : Node; np : Vector; cr,cd : integer32 )
                      return Node is

  -- DESCRIPTION :
  --   Returns the node obtained from swapping rows and after
  --   making it happy according to the next location of the
  --   black checkers in np.

    res : Node(nd.k) := Copy_Child(nd);

  begin
    Swap(res.rows,cr,cd);
    Make_Happy(np,res.rows,res.cols);
    return res;
  end Swap_Child;

  procedure Update_Children ( p : in out Poset; level : in integer32;
                              nd : in Node; lnd : out Link_to_Node ) is

  -- DESCRIPTION :
  --   Updates the nodes at p.white(level) with the node nd,
  --   a newly created child, which may already occur at the level.

  -- ON ENTRY :
  --   p        poset with at p.white a list of white checker positions;
  --   level    current level at the poset to consider;
  --   nd       node of some child generated during specialization.

  -- ON RETURN :
  --   p        poset with updated Littlewood-Richardson coefficients;
  --   lnd      link to the node whose content equals nd if already in p,
  --            otherwise lnd is a pointer to nd.

  begin
    if level <= p.white'last then
      declare
        tmp : Link_to_Node := p.white(level);
        previous : Link_to_Node;
      begin
        while tmp /= null loop
          if Equal(tmp.all,nd) then
            Add(tmp.coeff,nd.coeff);
            lnd := tmp;
            return;
          else
            previous := tmp;
            tmp := tmp.next_sibling;
          end if;
        end loop;
        lnd := new Node'(nd);
        if p.white(level) = null
         then p.white(level) := lnd;
         else previous.next_sibling := lnd;
        end if;
       end;
    else
      lnd := null;
    end if;
  end Update_Children;

  procedure White_Moves ( p : in out Poset; level,n,d,r : in integer32;
                          pb : in Vector; lnd : in Link_to_Node ) is

    cr : constant integer32
       := Critical_Row(integer32(pb(d)),n-d+1,lnd.rows,lnd.cols);
    cd : integer32;
    swap,stay : boolean;

  begin
    if cr = 0 then
      stay := true;
      swap := false;
    else
      cd := Top_White_Checker(integer32(pb(r)),n-r+1,n,lnd.rows,lnd.cols);
      if cd = 0 then
        stay := true;
        swap := false;
      else
       -- case Central_Choice(pb,d,r,lnd.rows,lnd.cols,cr,cd,true) is
        case Central_Choice(pb,d,r,lnd.rows,lnd.cols,cr,cd,false) is
          when 0 => stay := true; swap := false;
          when 1 => stay := false; swap := true;
          when others => stay := true; swap := true;
        end case;
      end if;
    end if;
    if stay then 
      declare
        stay_nd : constant Node(lnd.k)
                := Stay_Child(lnd.all,p.black(level+1).all);
        lk_stay : Link_to_Node;
      begin
        Update_Children(p,level+1,stay_nd,lk_stay); 
        lnd.stay_child := lk_stay;
        if lk_stay.first_parent = null
         then lk_stay.first_parent := lnd;
        end if;
      end;
    end if;
    if swap then
      declare
        swap_nd : constant Node(lnd.k)
                := Swap_Child(lnd.all,p.black(level+1).all,cr,cd);
        lk_swap : Link_to_Node;
      begin
        Update_Children(p,level+1,swap_nd,lk_swap); 
        lnd.swap_child := lk_swap;
        if lk_swap.first_parent = null
         then lk_swap.first_parent := lnd;
        end if;
      end;
    end if;
  end White_Moves;

  procedure Black_Moves ( ps : in out Poset; n,level : in integer32 ) is

    pb : constant Vector := ps.black(level).all;
    d : constant integer32 := Descending_Checker(pb);
    r : constant integer32 := Rising_Checker(pb,d);

  begin
    if r /= 0 then
      declare
        lnd : Link_to_Node := ps.white(level);
      begin
        while lnd /= null loop
          White_Moves(ps,level,n,d,r,pb,lnd);
          lnd := lnd.next_sibling;
        end loop;
      end;
    end if;
  end Black_Moves;

  function Create ( n,k : integer32; root : Node ) return Poset is

    res : Poset;

  begin
    res.black := new VecVec'(Specializing_Moves(n));
    res.white := new Array_of_Nodes(res.black'range);
    res.white(1) := new Node'(root);
    for i in res.white'first+1..res.white'last loop
      res.white(i) := null;
    end loop;
    for i in 1..res.black'last-1 loop
      Black_Moves(res,n,i);
    end loop;
    return res;
  end Create;

  function Create ( n : integer32; r,c : Vector ) return Poset is

    k : integer32 := r'last-r'first+1;
    root : constant Node(k) := Create(r,c);

  begin
    return Create(n,k,root);
  end Create;

  function Create ( n : integer32; cff : Natural_Number;
                    r,c : Vector ) return Poset is

    k : constant integer32 := r'last-r'first+1;
    root : Node(k) := Create(r,c);

  begin
    Copy(cff,root.coeff);
    return Create(n,k,root);
  end Create;

  procedure Add_Multiplicity ( ps : in out Poset; m : in Natural_Number ) is

    lnd : Link_to_Node;

  begin
    lnd := ps.white(ps.white'first);     -- add m to the root
    Add(lnd.coeff,m);
    for i in ps.white'first+1..ps.white'last loop
      lnd := ps.white(i);                -- multiplity of node sums up
      while lnd /= null loop             -- multiplicities of parents
        Clear(lnd.coeff);
        lnd.coeff := Multiplicity_of_Parents(lnd.all);
        lnd := lnd.next_sibling;
      end loop;
    end loop;
  end Add_Multiplicity;

  procedure Set_Coefficients_to_Zero ( ps : in out Poset ) is

    lnd : Link_to_Node;

  begin
    for i in ps.white'range loop
      lnd := ps.white(i);
      while lnd /= null loop 
        Clear(lnd.coeff);
        lnd.coeff := Create(natural32(0));
        lnd := lnd.next_sibling;
      end loop;
    end loop;
  end Set_Coefficients_to_Zero;

  procedure Add_from_Leaves_to_Root ( ps : in out Poset ) is

    lnd : Link_to_Node;

  begin
    for i in reverse ps.white'first..ps.white'last-1 loop
      lnd := ps.white(i);
      while lnd /= null loop
        if lnd.stay_child /= null
         then Add(lnd.coeff,lnd.stay_child.coeff);
        end if;
        if lnd.swap_child /= null
         then Add(lnd.coeff,lnd.swap_child.coeff);
        end if;
        lnd := lnd.next_sibling;
      end loop;
    end loop;
  end Add_from_Leaves_to_Root;

-- SELECTORS :

  function Root_Rows ( ps : in Poset ) return Vector is

    lnd : constant Link_to_Node := ps.white(ps.white'first);

  begin
    return lnd.rows;
  end Root_Rows;

  function Root_Columns ( ps : in Poset ) return Vector is

    lnd : constant Link_to_Node := ps.white(ps.white'first);

  begin
    return lnd.cols;
  end Root_Columns;

  function Equal ( nd1,nd2 : Node ) return boolean is
  begin
    if not Equal(nd1.rows,nd2.rows) then
      return false;
    elsif not Equal(nd1.cols,nd2.cols) then
      return false;
    else
      return true;
    end if;
  end Equal;

  function Position ( first_nd : Node; nd : Node ) return integer32 is

    ind : integer32;

  begin
    if Equal(first_nd,nd) then
      return 1;
    elsif first_nd.next_sibling = null then
      return 0;
    else
      ind := Position(first_nd.next_sibling.all,nd);
      if ind = 0
       then return 0;
       else return ind + 1;
      end if;
    end if;
  end Position;

  function Retrieve ( ps : Poset; i,j : integer32 ) return Link_to_Node is

    res : Link_to_Node := null;
    ptr : Link_to_Node;

  begin
    if (i >= ps.white'first) and (i <= ps.white'last) then
      ptr := ps.white(i);
      for k in 1..(j-1) loop
        exit when (ptr = null);
        ptr := ptr.next_sibling;      
      end loop;
      res := ptr;
    end if;
    return res;
  end Retrieve;

  procedure Retrieve_Leaf ( ps : in Poset; cols : in Vector;
                            ind : out integer32; lnd : out Link_to_Node ) is

    lst : constant integer32 := ps.white'last;
    ptr : Link_to_Node;

  begin
    ptr := ps.white(lst);
    ind := 1;
    loop
      exit when (ptr = null);
      if Equal(ptr.cols,cols)
       then lnd := ptr; return;
      end if;
      ptr := ptr.next_sibling;
      ind := ind + 1;
    end loop;
    lnd := ptr;
    ind := 0;
  end Retrieve_Leaf; 

  function Is_Stay_Child ( parent,child : Node ) return boolean is
  begin
    if parent.stay_child = null
     then return false;
     else return Equal(parent.stay_child.all,child);
    end if;
  end Is_Stay_Child;

  function Is_Swap_Child ( parent,child : Node ) return boolean is
  begin
    if parent.swap_child = null
     then return false;
     else return Equal(parent.swap_child.all,child);
    end if;
  end Is_Swap_Child;

  procedure Enumerate_Parents ( nd : in Node ) is

    parent : constant Link_to_Node := nd.first_parent;
    ptr : Link_to_Node;
    found_parent : boolean;  -- notice: parents are not always siblings

  begin
    if parent /= null then
      Process_Parent(parent);
      ptr := parent.next_sibling;
      while ptr /= null loop
        found_parent := false;
        if ptr.stay_child /= null then
          if Equal(ptr.stay_child.all,nd)
           then found_parent := true;
          end if;
        end if;
        if not found_parent then
          if ptr.swap_child /= null then
            if Equal(ptr.swap_child.all,nd)
             then found_parent := true;
            end if;
          end if;
        end if;
        if found_parent
         then Process_Parent(ptr);
        end if;
        ptr := ptr.next_sibling;
      end loop;
    end if;
  end Enumerate_Parents;

  function Number_of_Parents ( nd : Node ) return natural32 is

    res : natural32 := 0;

    procedure Count ( parent : in Link_to_Node ) is
    begin
      res := res + 1;
    end Count;
    procedure Count_Parents is new Enumerate_Parents(Count);

  begin
    Count_Parents(nd);
    return res;
  end Number_of_Parents;

  function Multiplicity_of_Parents ( nd : Node ) return Natural_Number is

    res : Natural_Number := Create(integer(0));

    procedure Sum_Multiplicity ( parent : in Link_to_Node ) is
    begin
      Add(res,parent.coeff);
    end Sum_Multiplicity;
    procedure Sum_Multiplicities is new Enumerate_Parents(Sum_Multiplicity);

  begin
    Sum_Multiplicities(nd);
    return res;
  end Multiplicity_of_Parents;

  function Retrieve_Parent ( nd : Node; k : integer32 ) return Link_to_Node is

    res : Link_to_Node := null;
    cnt : integer32 := 0;

    procedure Visit_Parent ( parent : in Link_to_Node ) is
    begin
      cnt := cnt + 1;
      if k = cnt
       then res := parent;
      end if;
    end Visit_Parent;
    procedure Visit_Parents is new Enumerate_Parents(Visit_Parent);

  begin
    Visit_Parents(nd);
    return res;
  end Retrieve_Parent;

  function Degree_of_Freedom ( ps : Poset ) return natural32 is

    pp : constant Link_to_Vector := ps.black(ps.black'first);
    nd : constant Link_to_Node := ps.white(ps.white'first);

  begin
    return Degree_of_Freedom(pp.all,nd.rows,nd.cols);
  end Degree_of_Freedom;

  procedure Enumerate_Paths_in_Poset ( ps : in Poset ) is

    path : Array_of_Nodes(ps.white'range);
    leaf : Link_to_Node := ps.white(ps.white'last);
    go : boolean := true;

    procedure Walk_to_Parents ( nd : in Node; level : in integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all parents at the current node for any valid level.

    begin
      if level > path'last then
        Process_Path(path,go);
      else
        for k in 1..integer32(Number_of_Parents(nd)) loop
          path(level) := Retrieve_Parent(nd,k);
          Walk_to_Parents(path(level).all,level+1);
          exit when not go;
        end loop;
      end if;
    end Walk_to_Parents;
 
  begin
    while leaf /= null loop
      path(path'first) := leaf;
      Walk_to_Parents(leaf.all,path'first+1);
      exit when not go;
      leaf := leaf.next_sibling;
    end loop;
  end Enumerate_Paths_in_Poset;

-- DESTRUCTORS :

  procedure Clear ( lnd : in out Link_to_Node ) is

    procedure free is new unchecked_deallocation(Node,Link_to_Node);

  begin
    if lnd /= null then
      Clear(lnd.coeff);
      free(lnd);
    end if;
  end Clear;

  procedure Clear ( arrlnd : in out Array_of_Nodes ) is
  begin
    for i in arrlnd'range loop
      Clear(arrlnd(i));
    end loop;
  end Clear;

  procedure Clear ( arrlnd : in out Link_to_Array_of_Nodes ) is

    procedure free is new
      unchecked_deallocation(Array_of_Nodes,Link_to_Array_of_Nodes);

  begin
    if arrlnd /= null then
      for i in arrlnd'range loop
        Clear(arrlnd(i));
      end loop;
      free(arrlnd);
    end if;
  end Clear;

  procedure Clear ( ps : in out Poset ) is
  begin
    Deep_Clear(ps.black);
    Clear(ps.white);
  end Clear;

end Checker_Posets;
