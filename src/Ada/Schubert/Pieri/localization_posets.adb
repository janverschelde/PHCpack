with unchecked_deallocation;

package body Localization_Posets is

-- NOTE :
--   The field nd.roco is set to -1 if all its children have been created.
--   This flag prevents traversing the poset needlessly.

-- CREATOR AUXILIARIES :

  function Last_Sibling ( root : Link_to_Node; level : integer32 )
                        return Link_to_Node is

  -- DESCRIPTION :
  --   Returns the last sibling at the level, or the empty pointer if
  --   there is no node at that level.

    res : Link_to_Node := null;
    sibnd : constant Link_to_Node := Find_Node(root,level);

    procedure Search_Next ( current : in Link_to_Node ) is
    begin
      if current.next_sibling = null
       then res := current;
       else Search_Next(current.next_sibling);
      end if;
    end Search_Next;

  begin
    if sibnd /= null
     then Search_Next(sibnd);
    end if;
    return res;
  end Last_Sibling;

  procedure Search_Sibling ( root : in Link_to_Node; nd : in Node;
                             lnd : out Link_to_Node; found : out boolean ) is

  -- DESCRIPTION :
  --   Searches the poset for the link to a node with contents nd.
  --   If found is true, then lnd is a pointer to that node, otherwise
  --   lnd points to the last sibling, or is empty when there is no
  --   node at level nd.level.         

    sibnd : constant Link_to_Node := Find_Node(root,nd.level);

    procedure Search_Next ( current : in Link_to_Node ) is
    begin
      if Equal(current.all,nd) then
        found := true;
        lnd := current;
      elsif current.next_sibling = null then
        found := false;
        lnd := current; 
      else
        Search_Next(current.next_sibling);
      end if;
    end Search_Next;

  begin
    if sibnd = null
     then lnd := sibnd; found := false;
     else Search_Next(sibnd);
    end if;
  end Search_Sibling;

  function Create_Child ( root : Link_to_Node; child : Node; share : boolean )
                        return Link_to_Node is

  -- DESCRIPTION :
  --   If the flag share is on, then the poset is searched for a node
  --   with the same contents as the child.  If a sibling is found,
  --   then the pointer to this sibling is returned, otherwise the link
  --   on return is a newly created link to node with contents child.
  --   If the flag share is off, then the link on return points to the
  --   last sibling node on that level, which has now contents child.

    res,lnd : Link_to_Node;
    found : boolean;

  begin
    if share then
      Search_Sibling(root,child,lnd,found);
      if found 
       then res := lnd;
      end if;
    else 
      lnd := Last_Sibling(root,child.level);
      found := false;
    end if;
    if not found then
      res := new Node'(child);
      if lnd /= null then
        lnd.next_sibling := res;
        res.prev_sibling := lnd;
      end if;
    end if;
    return res;
  end Create_Child;

  function Find_Index ( indexed_poset : Array_of_Array_of_Nodes;
                        nd : Link_to_Node ) return integer32 is

  -- DESCRIPTION :
  --   Returns 0 if the node does not occur at indexed_poset(nd.level),
  --   otherwise returns the index of the node nd in that array.
  --   Note that the pointers are compared to deal with sharing.

  begin
    if indexed_poset(nd.level) /= null then
      for i in indexed_poset(nd.level)'range loop
        if indexed_poset(nd.level)(i) = nd
         then return i;
        end if;
      end loop;
    end if;
    return 0;
  end Find_Index;

  function Labels_of_Children ( indexed_poset : Array_of_Array_of_Nodes;
                                nd : Node ) return Link_to_Vector is

  -- DESCRIPTION :
  --   Returns the labels of the children of the current node.

  -- REQUIRED : indexed_poset(i) created for i < nd.level.

    res : Link_to_Vector;
    nbc : constant natural32 := Number_of_Children(nd);
    cnt : integer32;

  begin
    if nbc /= 0 then
      res := new Standard_Natural_Vectors.Vector(1..integer32(nbc));
      cnt := 0;
      for i in nd.children'range(1) loop
        for j in nd.children'range(2) loop
          if nd.children(i,j) /= null then
            cnt := cnt+1;
            res(cnt) := natural32(Find_Index(indexed_poset,nd.children(i,j)));
          end if;
        end loop;
      end loop;
    end if;
    return res;
  end Labels_of_Children;

-- SPECIAL TEST FOR GENERAL QUANTUM PIERI RULE :

  function Special_Plane ( piv : Bracket; lag : natural32 ) return Bracket is

  -- DESCRIPTION :
  --   Returns the indices of the basis vectors that span the special
  --   m-dimensional plane, defined by the complementary indices in piv.

    res : Bracket(1..integer32(lag)-piv'last);
    ind : integer32 := 0;
    found : boolean;

  begin
    for i in 1..lag loop
      found := false;
      for j in piv'range loop
        found := (piv(j) = i);
        exit when found or (piv(j) > i);
      end loop;
      if not found then
        ind := ind+1;
        res(ind) := i;
      end if;
    end loop;
    return res;
  end Special_Plane;

  function Intersect_Spaces ( b1,b2 : Bracket ) return Bracket is

  -- DESCRIPTION :
  --   Returns the pivots that are common to both brackets.

    res : Bracket(b1'range);
    cnt : integer32 := 0;
    found : boolean;

  begin
    for i in b1'range loop
      found := false;
      for j in b2'range loop
        found := (b2(j) = b1(i));
        exit when found;
      end loop;
      if found then
        cnt := cnt+1;
        res(cnt) := b1(i);
      end if;
    end loop;
    return res(1..cnt);
  end Intersect_Spaces;

  function Merging_Top_Pivot_Test ( piv,spc : Bracket ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there exists a decreasing sequence of successive
  --   pivots from piv and spc that has length strictly higher than the
  --   value of the last pivot used, starting at the tails of the brackets.

    max : constant integer32 := piv'last + spc'last;
    acc : Bracket(1..max) := (1..max => 0);
    acc_ind : integer32 := max+1;
    piv_ind : integer32 := piv'last; 
    spc_ind : integer32 := spc'last;
    stop : boolean;

    procedure Merge ( fail : out boolean ) is

    -- DESCRIPTION :
    --   A consecutive pivot is added to the accumulator;
    --   failure is reported when such is not possible.

      procedure Add_from_Pivots is
      begin
        if (acc_ind = max+1) or else (piv(piv_ind) >= acc(acc_ind) - 1) then
          acc_ind := acc_ind-1;
          acc(acc_ind) := piv(piv_ind);
          piv_ind := piv_ind-1;
          fail := false;
        end if;
      end Add_from_Pivots;

      procedure Add_from_Space is
      begin
        if (acc_ind = max+1) or else (spc(spc_ind) >= acc(acc_ind) - 1) then
          acc_ind := acc_ind-1;
          acc(acc_ind) := spc(spc_ind);
          spc_ind := spc_ind-1;
          fail := false;
        end if;
      end Add_from_Space;
 
    begin
      fail := true;
      if piv_ind >= piv'first then
        if spc_ind >= spc'first then
          if piv(piv_ind) >= spc(spc_ind)
           then Add_from_Pivots;
           else Add_from_Space;
          end if;
        else 
          Add_from_Pivots;
        end if;
      elsif spc_ind >= spc'first then
        Add_from_Space;
      end if;
    end Merge;

  begin
    loop
      Merge(stop);
      if integer32(acc(acc_ind)) > (acc_ind + (integer32(acc(max)) - max))
       then return true;
      end if;
      exit when stop;
    end loop;
    return false;
  end Merging_Top_Pivot_Test;

  function Merging_Bottom_Pivot_Test ( piv,spc : Bracket ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there exists a increasing sequence of successive
  --   pivots from piv and spc that has length strictly higher than the
  --   value of the last pivot used, starting at the heads of the brackets.

    max : constant integer32 := piv'last + spc'last;
    acc : Bracket(1..max) := (1..max => 0);
    acc_ind : integer32 := 0;
    piv_ind : integer32 := piv'first; 
    spc_ind : integer32 := spc'first;
    stop : boolean;

    procedure Merge ( fail : out boolean ) is

    -- DESCRIPTION :
    --   A consecutive pivot is added to the accumulator;
    --   failure is reported when such is not possible.

      procedure Add_from_Pivots is
      begin
        if (acc_ind = 0) or else (piv(piv_ind) <= acc(acc_ind) + 1) then
          acc_ind := acc_ind+1;
          acc(acc_ind) := piv(piv_ind);
          piv_ind := piv_ind+1;
          fail := false;
        end if;
      end Add_from_Pivots;

      procedure Add_from_Space is
      begin
        if (acc_ind = 0) or else (spc(spc_ind) <= acc(acc_ind) + 1) then
          acc_ind := acc_ind+1;
          acc(acc_ind) := spc(spc_ind);
          spc_ind := spc_ind+1;
          fail := false;
        end if;
      end Add_from_Space;
 
    begin
      fail := true;
      if piv_ind <= piv'last then
        if spc_ind <= spc'last then
          if piv(piv_ind) <= spc(spc_ind)
           then Add_from_Pivots;
           else Add_from_Space;
          end if;
        else
          Add_from_Pivots;
        end if;
      elsif spc_ind <= spc'last then
        Add_from_Space;
      end if;
    end Merge;

  begin
    loop
      Merge(stop);
      if acc(acc_ind) < (natural32(acc_ind) + (acc(1) - 1))
       then return true;
      end if;
      exit when stop;
    end loop;
    return false;
  end Merging_Bottom_Pivot_Test;

-- CREATOR PRIMITIVES I : CHECK IF CREATION IS POSSIBLE AND ALLOWED

  function Top_Creatable ( nd : Node; n,i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the i-th top pivot can be incremented.
  --   The n is the dimension of the working space.    

  begin
    if nd.bottom(i) <= nd.top(i) then
      return false;
    elsif i = nd.p then
      return (nd.top(i) < natural32(n));
    else
      return (nd.top(i)+1 < nd.top(i+1));
    end if;
  end Top_Creatable;

  function Q_Top_Creatable ( nd : Node; n,lag,i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the i-th top pivot can be incremented.
  --   The n is the dimension of the working space.    

  begin
    if not Top_Creatable(nd,n,i) then
      return false;
    elsif i < nd.p then
      return true;
    else
      return (nd.top(nd.p) - nd.top(1) + 1 < natural32(lag));
    end if;
  end Q_Top_Creatable;

  function Q_Top_Creatable
               ( nd : Node; modtop,space : Bracket; n,lag,pi,i : integer32 )
               return boolean is

  -- DESCRIPTION :
  --   This is the quantum analogue to implement the modular bottom-left
  --   rule as needed in the general intersection case.

  -- ON ENTRY :
  --   nd        current node;
  --   modtop    top pivots of nd, modulo the lag;
  --   space     generators of the intersection of special m-planes;
  --   n         dimension of the working space;
  --   lag       equals m+p;
  --   pi        index in nd.top, permuted index i used to sort modtop;
  --   i         modtop(i) will be increased to derive the child.

    child : Bracket(modtop'range) := modtop;

  begin
    if not Q_Top_Creatable(nd,n,lag,pi) then             -- valid pattern ?
      return false;
    else -- valid pattern => valid child, only last entry might be zero
      child(i) := modtop(i)+1;   
      if i = child'last and child(i) = natural32(lag)+1 then
        for j in reverse child'first+1..child'last loop
          child(j) := child(j-1);
        end loop;
        child(child'first) := 1;
      end if;
      return Merging_Top_Pivot_Test(child,space);
    end if;
  end Q_Top_Creatable;

  function Bottom_Creatable ( nd : Node; i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the i-th bottom pivot can be decremented.

  begin
    if nd.bottom(i) <= nd.top(i) then
      return false;
    elsif i = 1 then
      return (nd.bottom(i) > 1);
    else
      return (nd.bottom(i)-1 > nd.bottom(i-1));
    end if;
  end Bottom_Creatable;

  function Q_Bottom_Creatable ( nd : Node; lag,i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the i-th bottom pivot can be decremented and if
  --   the spacing between first and last bottom pivot will remain < lag.

  begin
    if not Bottom_Creatable(nd,i) then
      return false;
    elsif i > 1 then
      return true;
    else
      return (nd.bottom(nd.p) - nd.bottom(1) + 1 < natural32(lag));
    end if;
  end Q_Bottom_Creatable;

  function Q_Bottom_Creatable
               ( nd : Node; modbot,space : Bracket; lag,pi,i : integer32 )
               return boolean is

  -- DESCRIPTION :
  --   This is the quantum analogue to implement the modular bottom-left
  --   rule as needed in the general intersection case.

  -- ON ENTRY :
  --   nd        current node;
  --   modbot    bottom pivots of nd, modulo the lag;
  --   space     generators of the intersection of special m-planes;
  --   lag       equals m+p;
  --   pi        index in nd.bottom, permuted index i used to sort modbot;
  --   i         modbot(i) will be decreased to derive the child.

    child : Bracket(modbot'range) := modbot;

  begin
    if not Q_Bottom_Creatable(nd,lag,pi) then       -- valid pattern ?
      return false;
    else -- valid pattern => valid child, only 1st entry might be zero
      child(i) := modbot(i)-1;   
      if i = 1 and child(i) = 0 then
        for j in child'first..child'last-1 loop
          child(j) := child(j+1);
        end loop;
        child(child'last) := natural32(lag);
      end if;
      return Merging_Bottom_Pivot_Test(child,space);
    end if;
  end Q_Bottom_Creatable;

  function Top_Bottom_Creatable ( nd : Node; n,i,j : integer32 )
                                return boolean is

  -- DESCRIPTION :
  --   Returns true if the i-th top pivot can be incremented and if
  --   the j-th bottom pivot can be decremented.        

  begin
    if not Top_Creatable(nd,n,i) then
      return false;
    elsif not Bottom_Creatable(nd,j) then
      return false;
    elsif i /= j then
      return true;
    else
      return (nd.bottom(i) - nd.top(i) > 1);
    end if;
  end Top_Bottom_Creatable;

  function Q_Top_Bottom_Creatable ( nd : Node; n,lag,i,j : integer32 )
                                  return boolean is

  -- DESCRIPTION :
  --   Returns true if the i-th top pivot can be incremented and if
  --   the j-th bottom pivot can be decremented.        

  begin
    if not Q_Top_Creatable(nd,n,lag,i) then
      return false;
    elsif not Q_Bottom_Creatable(nd,lag,j) then
      return false;
    elsif i /= j then
      return true;
    else
      return (nd.bottom(i) - nd.top(i) > 1);
    end if;
  end Q_Top_Bottom_Creatable;

  function Q_Top_Bottom_Creatable
              ( nd : Node; modtop,topspc,modbot,botspc : Bracket;
                n,lag,pi,i,pj,j : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the i-th top pivot can be incremented and if
  --   the j-th bottom pivot can be decremented in the general quantum
  --   Pieri homotopy algorithm.

  begin
    if not Q_Top_Creatable(nd,modtop,topspc,n,lag,pi,i) then
      return false;
    elsif not Q_Bottom_Creatable(nd,modbot,botspc,lag,pj,j) then
      return false;
    elsif pi /= pj then
      return true;
    else
      return (nd.bottom(pi) - nd.top(pi) > 1);
    end if;
  end Q_Top_Bottom_Creatable;

-- CREATOR PRIMITIVES II : DERIVE CHILD FROM NODE

  procedure Create_Top_Child ( root,nd : in Link_to_Node;
                               i : in integer32; share : in boolean ) is

  -- DESCRIPTION :
  --   Creates a child of the given node by incrementing the i-th top pivot.

    child : Node(nd.p);

  begin
    child.level := nd.level-1;
    child.roco := 0;
    child.bottom := nd.bottom;
    child.top := nd.top;
    child.top(i) := nd.top(i)+1;
    nd.children(i,0) := Create_Child(root,child,share);
  end Create_Top_Child;

  procedure Create_Bottom_Child
              ( root : in Link_to_Node; nd : in Link_to_Node;
                i : in integer32; share : in boolean ) is

  -- DESCRIPTION :
  --   Creates a child of the node nd by decrementing the i-th bottom pivot.

    child : Node(nd.p);

  begin
    child.level := nd.level-1;
    child.roco := 0;
    child.bottom := nd.bottom;
    child.top := nd.top;
    child.bottom(i) := nd.bottom(i)-1;
    nd.children(0,i) := Create_Child(root,child,share);
  end Create_Bottom_Child;

  procedure Create_Top_Bottom_Child
               ( root,nd : in Link_to_Node;
                 i,j : in integer32; share : in boolean ) is
  
  -- DESCRIPTION :
  --   Creates a child of the node nd by incrementing the i-th top pivot
  --   and decrementing the i-th bottom pivot.

    child : Node(nd.p);

  begin
    child.level := nd.level-2;
    child.roco := 0;
    child.top := nd.top;
    child.top(i) := nd.top(i)+1;
    child.bottom := nd.bottom;
    child.bottom(j) := nd.bottom(j)-1;
    nd.children(i,j) := Create_Child(root,child,share);
  end Create_Top_Bottom_Child;

-- CREATOR PRIMITIVES III : TREAT ONE/TWO DEGREE(S) OF FREEDOM

  procedure Top_Create1 ( root,nd : in Link_to_Node; n : in integer32 ) is

  -- DESCRIPTION :
  --   Creates new nodes by incrementing the top pivots, bounded by n.
  --   The levels of the children nodes decrease by one as this is the
  --   hypersurface case.

  begin
    nd.tp := top;
    for i in nd.top'range loop
      if Top_Creatable(nd.all,n,i)
       then Create_Top_Child(root,nd,i,true);
      end if;
    end loop;
  end Top_Create1;

  procedure Q_Top_Create1 ( root,nd : in Link_to_Node;
                            n,lag : in integer32 ) is

  -- DESCRIPTION :
  --   Creates new nodes by incrementing the top pivots, for general q,
  --   where we need the parameters n = dimension of working space
  --   and lag = m+p, to bound the space between first and last entry.

  begin
    nd.tp := top;
    for i in nd.top'range loop
      if Q_Top_Creatable(nd.all,n,lag,i)
       then Create_Top_Child(root,nd,i,true);
      end if;
    end loop;
  end Q_Top_Create1;

  procedure Top_Create1 ( root,nd : in Link_to_Node;
                          k,n,c : in integer32 ) is

  -- DESCRIPTION :
  --   Does k steps of the other Top_Create1 taking pivots larger than c.
  --   This is the general case, for k=1 we have the hypersurface case.

    share : constant boolean := (k = 1);

  begin
    if k > 0 then
      nd.tp := top;
      for i in c..nd.top'last loop
        if Top_Creatable(nd.all,n,i) then
          Create_Top_Child(root,nd,i,share);
          if k > 1
           then Top_Create1(root,nd.children(i,0),k-1,n,i);
          end if;
        end if;
      end loop;
    end if;
  end Top_Create1;

  procedure Q_Top_Create1 ( root,nd : in Link_to_Node;
                            first : in boolean; space : in Bracket;
                            k,n,lag : in integer32 ) is

  -- DESCRIPTION :
  --   Does k steps in a top-right chain on modular brackets.
  --   The top-right rule is enforced by the merging pivot test involving
  --   top pivots and the indices of the vectors that span the space of
  --   intersection of special m-planes.

  -- ON ENTRY :
  --   root       root of the poset where the construction started;
  --   nd         current node;
  --   first      if true, then this is the first step in the sequence,
  --              and the space has yet to be determined;
  --   space      contains generators of the intersection of special m-planes;
  --   k          number of steps still left to do;
  --   n          dimension of the space;
  --   lag        m+p.

    share : constant boolean := (k=1);
    modtop : Bracket(nd.top'range);
    perm : Standard_Natural_Vectors.Vector(modtop'range);
    special : Bracket(1..lag-nd.p);

    procedure Recursive_Top_Create1 ( new_space : in Bracket ) is

    -- DESCRIPTION :
    --   Additional layer needed for the determination of the updated space.

    begin
      for i in modtop'range loop
        if Q_Top_Creatable(nd.all,modtop,new_space,n,lag,
                           integer32(perm(i)),i) then
          Create_Top_Child(root,nd,integer32(perm(i)),share);
          if k > 1 then
            Q_Top_Create1(root,nd.children(integer32(perm(i)),0),
                          false,new_space,k-1,n,lag);
          end if;
        end if;
      end loop;
    end Recursive_Top_Create1;

  begin
    if k > 0 then
      nd.tp := top;
      Modulo(nd.top,natural32(lag),perm,modtop);
      special := Special_Plane(modtop,natural32(lag));
      if first then
        Recursive_Top_Create1(special);
      else
        declare
          int_spc : constant Bracket := Intersect_Spaces(space,special);
        begin
          Recursive_Top_Create1(int_spc);
        end;
      end if;
    end if;
  end Q_Top_Create1;

  procedure Bottom_Create1 ( root,nd : in Link_to_Node ) is

  -- DESCRIPTION :
  --   Creates new nodes by decrementing the bottom pivots.
  --   The levels of the children nodes decrease by one as this is
  --   the hypersurface case.

  begin
    nd.tp := bottom;
    for i in nd.top'range loop
      if Bottom_Creatable(nd.all,i)
       then Create_Bottom_Child(root,nd,i,true);
      end if;
    end loop;
  end Bottom_Create1;

  procedure Q_Bottom_Create1
                ( root,nd : in Link_to_Node; lag : in integer32 ) is

  -- DESCRIPTION :
  --   Creates new nodes by decrementing the bottom pivots for general q,
  --   where the parameter lag > max space between first and last entry.

  begin
    nd.tp := bottom;
    for i in nd.top'range loop
      if Q_Bottom_Creatable(nd.all,lag,i)
       then Create_Bottom_Child(root,nd,i,true);
      end if;
    end loop;
  end Q_Bottom_Create1;

  procedure Bottom_Create1 ( root,nd : in Link_to_Node;
                             k,c : in integer32 ) is

  -- DESCRIPTION :
  --   Does k steps of the other Bottom_Create1 taking pivots smaller than c.
  --   This is the general case, for k=1 we have the hypersurface case.

    share : constant boolean := (k=1);

  begin
    if k > 0 then
      nd.tp := bottom;
      for i in nd.bottom'first..c loop
        if Bottom_Creatable(nd.all,i) then
          Create_Bottom_Child(root,nd,i,share);
          if k > 1
           then Bottom_Create1(root,nd.children(0,i),k-1,i);
          end if;
        end if;
      end loop;
    end if;
  end Bottom_Create1;

  procedure Q_Bottom_Create1 ( root,nd : in Link_to_Node;
                               first : in boolean; space : in Bracket;
                               k,lag : in integer32 ) is

  -- DESCRIPTION :
  --   Does k steps in a bottom-left chain on modular brackets.
  --   The bottom-left rule is enforced by the merging pivot test involving
  --   bottom pivots and the indices of the vectors that span the space of
  --   intersection of special m-planes.

  -- ON ENTRY :
  --   root       root of the poset where the construction started;
  --   nd         current node;
  --   first      if true, then this is the first step in the sequence,
  --              and the space has yet to be determined;
  --   space      contains generators of the intersection of special m-planes;
  --   k          number of steps still left to do;
  --   lag        m+p.

    share : constant boolean := (k=1);
    modbot : Bracket(nd.bottom'range);
    perm : Standard_Natural_Vectors.Vector(modbot'range);
    special : Bracket(1..lag-nd.p);

    procedure Recursive_Bottom_Create1 ( new_space : in Bracket ) is

    -- DESCRIPTION :
    --   Additional layer needed for the determination of the updated space.

    begin
      for i in modbot'range loop
        if Q_Bottom_Creatable(nd.all,modbot,new_space,lag,
                              integer32(perm(i)),i) then
          Create_Bottom_Child(root,nd,integer32(perm(i)),share);
          if k > 1 then
            Q_Bottom_Create1(root,nd.children(0,integer32(perm(i))),
                             false,new_space,k-1,lag);
          end if;
        end if;
      end loop;
    end Recursive_Bottom_Create1;

  begin
    if k > 0 then
      nd.tp := bottom;
      Modulo(nd.bottom,natural32(lag),perm,modbot);
      special := Special_Plane(modbot,natural32(lag));
      if first then
        Recursive_Bottom_Create1(special);
      else
        declare
          int_spc : constant Bracket := Intersect_Spaces(space,special);
        begin
          Recursive_Bottom_Create1(int_spc);
        end;
      end if;
    end if;
  end Q_Bottom_Create1;

  procedure Top_Bottom_Create1 ( root,nd : in Link_to_Node;
                                 n : in integer32 ) is

  -- DESCRIPTION :
  --   Creates new nodes by incrementing top pivots and decrementing bottom
  --   pivots, with n the maximal entry in any pivot.
  --   If no top create is possible, then a bottom create will be done,
  --   and we have only a bottom create when no top create is possible.    

    nocreate : boolean := true;

  begin
    nd.tp := mixed;
    for i in nd.top'range loop                      -- first do top+bottom
      for j in nd.bottom'range loop
        if Top_Bottom_Creatable(nd.all,n,i,j) then
          Create_Top_Bottom_Child(root,nd,i,j,true);
          nocreate := false;
        end if;
      end loop;
    end loop;
    if nocreate then                       -- no top+bottom create possible
      Bottom_Create1(root,nd);
      if Is_Leaf(nd.all)                   -- no bottom create possible
       then Top_Create1(root,nd,n);
      end if;
    end if;
  end Top_Bottom_Create1;

  procedure Q_Top_Bottom_Create1 ( root,nd : in Link_to_Node;
                                   n,lag : in integer32 ) is

  -- DESCRIPTION :
  --   Creates new nodes by incrementing top pivots and decrementing bottom
  --   pivots, with n the maximal entry in any pivot.
  --   If no top create is possible, then a bottom create will be done,
  --   and we have only a bottom create when no top create is possible.    

    nocreate : boolean := true;

  begin
    nd.tp := mixed;
    for i in nd.top'range loop                      -- first do top+bottom
      for j in nd.bottom'range loop
        if Q_Top_Bottom_Creatable(nd.all,n,lag,i,j) then
          Create_Top_Bottom_Child(root,nd,i,j,true);
          nocreate := false;
        end if;
      end loop;
    end loop;
    if nocreate then                       -- no top+bottom create possible
      Q_Bottom_Create1(root,nd,lag);
      if Is_Leaf(nd.all)                   -- no bottom create possible
       then Q_Top_Create1(root,nd,n,lag);
      end if;
    end if;
  end Q_Top_Bottom_Create1;

  procedure Top_Bottom_Create1 ( root : in out Link_to_Node;
                                 nd : in Link_to_Node;
                                 k1,k2,n,c1,c2 : in integer32 ) is

  -- DESCRIPTION :
  --   Applies the hypersurface Top_Bottom_Create max(k1,k2) times,
  --   taking top pivots in c1..p and bottom pivots in 1..c2.
  --   This is the top-bottom create that takes the codimensions in pairs,
  --   which allows more possibilities for sharing.

    share : constant boolean := ((k1=1) and (k2=1));

  begin
    if (k1 > 0) and (k2 > 0) then
      nd.tp := mixed;
      for i in c1..nd.top'last loop                    -- first do top+bottom
        for j in nd.bottom'first..c2 loop
          if Top_Bottom_Creatable(nd.all,n,i,j) then
            Create_Top_Bottom_Child(root,nd,i,j,share);
            if ((k1 > 1) or (k2 > 1))
             then Top_Bottom_Create1(root,nd.children(i,j),k1-1,k2-1,n,i,j);
            end if;
          end if;
        end loop;
      end loop;
    end if;
    if ((k1 = 0) and (k2 > 0)) then
      Bottom_Create1(root,nd,k2,c2);
    elsif ((k1 > 0) and (k2 = 0)) then
      Top_Create1(root,nd,k1,n,c1);
    end if;
  end Top_Bottom_Create1;

  procedure Recursive_Top_Bottom_Create
              ( root,nd : in Link_to_Node;
                codim : in Bracket; ind,k1,k2,n,c1,c2 : in integer32;
                hyper : in boolean ) is

  -- DESCRIPTION :
  --   Applies the hypersurface Top_Bottom_Create max(k1,k2) times,
  --   taking top pivots in c1..p and bottom pivots in 1..c2.
  --   In case k1 and/or k2 are zero, new conditions will be treated.

  -- ON ENTRY :
  --   root     root of the localization poset;
  --   nd       current node;
  --   codim    list of co-dimension conditions;
  --   ind      index of lowest condition being treated;
  --   k1       co-dimension condition satisfied decrementing top pivots;
  --   k2       co-dimension condition satisfied incrementing bottom pivots;
  --   n        dimension of the working space;
  --   c1       needed to enforce the top-right rule;
  --   c2       needed to enforce the bottom-left rule;
  --   hyper    indicates whether or not in the hypersurface case.

    newhyper : boolean;

  begin
    if (k1 > 0) and (k2 > 0) then
      nd.tp := mixed;
      for i in c1..nd.top'last loop                    -- first do top+bottom
        for j in nd.bottom'first..c2 loop
          if Top_Bottom_Creatable(nd.all,n,i,j) then
            Create_Top_Bottom_Child(root,nd,i,j,hyper);
            Recursive_Top_Bottom_Create
              (root,nd.children(i,j),codim,ind,k1-1,k2-1,n,i,j,false);
          end if;
        end loop;
      end loop;
      nd.roco := -1;
    else
      if ((k1 = 0) and (k2 > 0)) then
        if ind > codim'first then
          Recursive_Top_Bottom_Create
            (root,nd,codim,ind-1,integer32(codim(ind-1)),k2,n,1,c2,false);
        else
          Bottom_Create1(root,nd,k2,c2);
        end if;
      elsif ((k1 > 0) and (k2 = 0)) then
        if ind > codim'first then
          Recursive_Top_Bottom_Create
            (root,nd,codim,ind-1,k1,integer32(codim(ind-1)),n,c1,nd.p,false);
        else
          Top_Create1(root,nd,k1,n,c1);
        end if;
      else -- k1 = 0 and k2 = 0
        if ind > codim'first + 1 then
          newhyper := ((codim(ind-2) = 1) and (codim(ind-1) = 1));
          Recursive_Top_Bottom_Create
            (root,nd,codim,ind-2,integer32(codim(ind-2)),
             integer32(codim(ind-1)),n,1,nd.p,newhyper);
        elsif ind > codim'first then
          Bottom_Create1(root,nd,integer32(codim(ind-1)),nd.p);
        end if;
      end if;
    end if;
  end Recursive_Top_Bottom_Create;

  procedure Q_Recursive_Top_Bottom_Create
              ( root,nd : in Link_to_Node; codim : in Bracket;
                fsttop : in boolean; topspc : in Bracket;
                fstbot : in boolean; botspc : in Bracket;
                ind,k1,k2,n,lag : in integer32; hyper : in boolean ) is

  -- DESCRIPTION :
  --   Applies the hypersurface Q_Top_Bottom_Create max(k1,k2) times,
  --   simulating the bottom-left and top-right rules with the modular
  --   brackets and corresponding spaces.

  -- ON ENTRY :
  --   root     root of the localization poset;
  --   nd       current node;
  --   codim    list of co-dimension conditions;
  --   fsttop   if true, then first step taken using top pivots;
  --   topspc   intersection of special m-planes for top pivots;
  --   fstbot   if true, then first step taken using bottom pivots;
  --   botspc   intersection of special m-planes for bottom pivots;
  --   ind      index of lowest condition being treated;
  --   k1       co-dimension condition satisfied decrementing top pivots;
  --   k2       co-dimension condition satisfied incrementing bottom pivots;
  --   n        dimension of the working space;
  --   lag      space in the poset that is of interest;
  --   hyper    indicates whether or not in the hypersurface case.

    newhyper : boolean;
    modtop,modbot : Bracket(1..nd.p);
    topprm,botprm : Standard_Natural_Vectors.Vector(1..nd.p);
    top_special,bot_special : Bracket(1..lag-nd.p);

    procedure Mixed_Create ( new_top_space,new_bot_space : in Bracket ) is
    begin
      for i in modtop'range loop
        for j in modbot'range loop
          if Q_Top_Bottom_Creatable
               (nd.all,modtop,new_top_space,modbot,new_bot_space,
                n,lag,integer32(topprm(i)),i,integer32(botprm(j)),j)
           then Create_Top_Bottom_Child
                  (root,nd,integer32(topprm(i)),integer32(botprm(j)),hyper);
                Q_Recursive_Top_Bottom_Create
                  (root,nd.children(integer32(topprm(i)),integer32(botprm(j))),
                   codim,false,new_top_space,false,new_bot_space,
                   ind,k1-1,k2-1,n,lag,false);
          end if;
        end loop;
      end loop;
      nd.roco := -1;
    end Mixed_Create;

  begin
    if (k1 > 0) and (k2 > 0) then  -- first do top + bottom
      nd.tp := mixed;
      Modulo(nd.top,natural32(lag),topprm,modtop);
      top_special := Special_Plane(modtop,natural32(lag));
      Modulo(nd.bottom,natural32(lag),botprm,modbot);
      bot_special := Special_Plane(modbot,natural32(lag));
      if fsttop then
        if fstbot then
          Mixed_Create(top_special,bot_special);
        else
          declare
            int_spc : constant Bracket := Intersect_Spaces(botspc,bot_special);
          begin
            Mixed_Create(top_special,int_spc);
          end;
        end if;
      elsif fstbot then
        declare
          int_spc : constant Bracket := Intersect_Spaces(topspc,top_special);
        begin
          Mixed_Create(int_spc,bot_special);
        end;
      else
        declare
          int_top : constant Bracket := Intersect_Spaces(topspc,top_special);
          int_bot : constant Bracket := Intersect_Spaces(botspc,bot_special);
        begin
          Mixed_Create(int_top,int_bot);
        end;
      end if;
    elsif ((k1 = 0) and (k2 > 0)) then
      if ind > codim'first
       then Q_Recursive_Top_Bottom_Create
               (root,nd,codim,true,topspc,fstbot,botspc,
                ind-1,integer32(codim(ind-1)),k2,n,lag,false);
       else Q_Bottom_Create1(root,nd,fstbot,botspc,k2,lag);
      end if;
    elsif ((k1 > 0) and (k2 = 0)) then
      if ind > codim'first
       then Q_Recursive_Top_Bottom_Create
              (root,nd,codim,fsttop,topspc,true,botspc,
               ind-1,k1,integer32(codim(ind-1)),n,lag,false);
       else Q_Top_Create1(root,nd,fsttop,topspc,k1,n,lag);
      end if;
    else -- k1 = 0 and k2 = 0
      if ind > codim'first + 1 then
        newhyper := ((codim(ind-2) = 1) and (codim(ind-1) = 1));
        Q_Recursive_Top_Bottom_Create
          (root,nd,codim,true,topspc,true,botspc, ind-2,
           integer32(codim(ind-2)),integer32(codim(ind-1)),n,lag,newhyper);
      elsif ind > codim'first then
        Q_Bottom_Create1(root,nd,true,botspc,integer32(codim(ind-1)),lag);
      end if;
    end if;
  end Q_Recursive_Top_Bottom_Create;

-- TARGET CREATORS :

  function Trivial_Root ( m,p : natural32 ) return Node is

    nd : Node(integer32(p));

  begin
    nd.level := integer32(m*p);
    nd.roco := 0;
    for i in 1..integer32(p) loop
      nd.top(i) := natural32(i);
      nd.bottom(i) := m+natural32(i);
    end loop;
    return nd;
  end Trivial_Root;

  function Trivial_Root ( m,p,q : natural32 ) return Node is

    nd : Node(integer32(p));
    last : natural32;

  begin
    if q = 0 then
      nd := Trivial_Root(m,p);
    else
      nd := Trivial_Root(m,p,q-1);
      nd.level := nd.level + integer32(m+p);
      last := nd.bottom(1)+m+p;
      for i in 1..integer32(p-1) loop
        nd.bottom(i) := nd.bottom(i+1);
      end loop;
      nd.bottom(integer32(p)) := last;
    end if;
    return nd;
  end Trivial_Root;

  procedure Top_Create ( root : in Link_to_Node; n : in natural32 ) is

    procedure Create_Next ( root,nd : in Link_to_Node ) is
    begin
      if ((nd.level > 0) and (nd.roco >= 0)) then
        Top_Create1(root,nd,integer32(n));
        for i in nd.children'range(1) loop
          if nd.children(i,0) /= null
           then Create_Next(root,nd.children(i,0));
          end if;
        end loop;
        nd.roco := -1;
      end if;
    end Create_Next;

  begin
    Create_Next(root,root);
  end Top_Create;

  procedure Q_Top_Create ( root : in out Link_to_Node;
                           n,lag : in natural32 ) is

    procedure Create_Next ( root : in out Link_to_Node;
                            nd : in Link_to_Node ) is
    begin
      if ((nd.level > 0) and (nd.roco >= 0)) then
        Q_Top_Create1(root,nd,integer32(n),integer32(lag));
        for i in nd.children'range(1) loop
          if nd.children(i,0) /= null
           then Create_Next(root,nd.children(i,0));
          end if;
        end loop;
        nd.roco := -1;
      end if;
    end Create_Next;

  begin
    Create_Next(root,root);
  end Q_Top_Create;

  procedure Top_Create ( root : in Link_to_Node;
                         k : in Bracket; n : in natural32 ) is

    procedure Create ( current : in Link_to_Node; ind : in integer32 );

    -- DESCRIPTION :
    --   Creates k(ind) levels above the current node.

    procedure Create_Children ( child : in Link_to_Node;
                                cnt,ind : in integer32 ) is

    -- DESCRIPTION :
    --   Goes to the topmost child to create, counting down with cnt.

    begin
      if cnt = 0 then
        Create(child,ind);
      else
        for i in child.children'range(1) loop
          if child.children(i,0) /= null
           then Create_Children(child.children(i,0),cnt-1,ind);
          end if;
        end loop;
      end if;
    end Create_Children;

    procedure Create ( current : in Link_to_Node; ind : in integer32 ) is
    begin
      if ((current.level > 0) and (ind <= k'last) and (current.roco >= 0)) then
        Top_Create1(root,current,integer32(k(ind)),integer32(n),1);
        if ind > k'first then
          for i in current.children'range(1) loop
            if current.children(i,0) /= null then
              Create_Children(current.children(i,0),integer32(k(ind))-1,ind-1);
            end if;
          end loop;
        end if;
      current.roco := -1;
      end if;
    end Create;

  begin
    Create(root,k'last);
  end Top_Create;

  procedure Q_Top_Create ( root : in Link_to_Node;
                           k : in Bracket; n,lag : in natural32 ) is

    procedure Create ( current : in Link_to_Node; ind : in integer32 );

    -- DESCRIPTION :
    --   Creates k(ind) levels above the current node.

    procedure Create_Children ( child : in Link_to_Node;
                                cnt,ind : in integer32 ) is

    -- DESCRIPTION :
    --   Goes to the topmost child to create, counting down with cnt.

    begin
      if cnt = 0 then
        Create(child,ind);
      else
        for i in child.children'range(1) loop
          if child.children(i,0) /= null
           then Create_Children(child.children(i,0),cnt-1,ind);
          end if;
        end loop;
      end if;
    end Create_Children;

    procedure Create ( current : in Link_to_Node; ind : in integer32 ) is

      space : constant Bracket(1..integer32(lag)-current.p)
            := (1..integer32(lag)-current.p => 0);

    begin
      if ((current.level > 0) and (ind <= k'last) and (current.roco >= 0)) then
        Q_Top_Create1(root,current,true,space,integer32(k(ind)),
                      integer32(n),integer32(lag));
        if ind > k'first then
          for i in current.children'range(1) loop
            if current.children(i,0) /= null then
              Create_Children(current.children(i,0),integer32(k(ind))-1,ind-1);
            end if;
          end loop;
        end if;
        current.roco := -1;
      end if;
    end Create;

  begin
    Create(root,k'last);
  end Q_Top_Create;

  procedure Bottom_Create ( root : in Link_to_Node ) is

    procedure Create_Next ( root,nd : in Link_to_Node ) is
    begin
      if ((nd.level > 0) and (nd.roco >= 0)) then
        Bottom_Create1(root,nd);
        for i in nd.children'range(2) loop
          if nd.children(0,i) /= null
           then Create_Next(root,nd.children(0,i));
          end if;
        end loop;
        nd.roco := -1;
      end if;
    end Create_Next;

  begin
    Create_Next(root,root);
  end Bottom_Create;

  procedure Q_Bottom_Create ( root : in Link_to_Node;
                              lag : in natural32 ) is

    procedure Create_Next ( root,nd : in Link_to_Node ) is
    begin
      if ((nd.level > 0) and (nd.roco >= 0)) then
        Q_Bottom_Create1(root,nd,integer32(lag));
        for i in nd.children'range(2) loop
          if nd.children(0,i) /= null
           then Create_Next(root,nd.children(0,i));
          end if;
        end loop;
        nd.roco := -1;
      end if;
    end Create_Next;

  begin
    Create_Next(root,root);
  end Q_Bottom_Create;

  procedure Bottom_Create ( root : in Link_to_Node; k : in Bracket ) is

    procedure Create ( current : in Link_to_Node; ind : in integer32 );

    -- DESCRIPTION :
    --   Creates k(ind) levels above the current node.

    procedure Create_Children ( child : in Link_to_Node;
                                cnt,ind : in integer32 ) is

    -- DESCRIPTION :
    --   Goes to the topmost child to create, counting down with cnt.

    begin
      if cnt = 0 then
        Create(child,ind);
      else 
        for i in child.children'range(1) loop
          if child.children(0,i) /= null
           then Create_Children(child.children(0,i),cnt-1,ind);
          end if;
        end loop;
      end if;
    end Create_Children;

    procedure Create ( current : in Link_to_Node; ind : in integer32 ) is
    begin
      if ((current.level > 0) and (ind <= k'last) and (current.roco >= 0)) then
        Bottom_Create1(root,current,integer32(k(ind)),current.p);
        if ind > k'first then
          for i in current.children'range(1) loop
            if current.children(0,i) /= null then
              Create_Children(current.children(0,i),integer32(k(ind))-1,ind-1);
            end if;
          end loop;
        end if;
        current.roco := -1;
      end if;
    end Create;

  begin
    Create(root,k'last);
  end Bottom_Create;

  procedure Q_Bottom_Create ( root : in Link_to_Node; k : in Bracket;
                              lag : in natural32 ) is

    procedure Create ( current : in Link_to_Node; ind : in integer32 );

    -- DESCRIPTION :
    --   Creates k(ind) levels above the current node.

    procedure Create_Children ( child : in Link_to_Node;
                                cnt,ind : in integer32 ) is

    -- DESCRIPTION :
    --   Goes to the topmost child to create, counting down with cnt.

    begin
      if cnt = 0 then
        Create(child,ind);
      else
        for i in child.children'range(1) loop
          if child.children(0,i) /= null
           then Create_Children(child.children(0,i),cnt-1,ind);
          end if;
        end loop;
      end if;
    end Create_Children;

    procedure Create ( current : in Link_to_Node; ind : in integer32 ) is

      space : constant Bracket(1..integer32(lag)-current.p)
            := (1..integer32(lag)-current.p => 0);

    begin
      if ((current.level > 0) and (ind <= k'last) and (current.roco >= 0)) then
        Q_Bottom_Create1(root,current,true,space,
                         integer32(k(ind)),integer32(lag));
        if ind > k'first then
          for i in current.children'range(1) loop
            if current.children(0,i) /= null then
              Create_Children(current.children(0,i),integer32(k(ind))-1,ind-1);
            end if;
          end loop;
        end if;
        current.roco := -1;
      end if;
    end Create;

  begin
    Create(root,k'last);
  end Q_Bottom_Create;

  procedure Top_Bottom_Create ( root : in Link_to_Node;
                                n : in natural32 ) is

    procedure Create_Next ( root,nd : in Link_to_Node ) is
    begin
      if ((nd.level > 0) and (nd.roco >= 0)) then
        Top_Bottom_Create1(root,nd,integer32(n));
        for i in nd.children'range(1) loop
          for j in nd.children'range(2) loop
            if nd.children(i,j) /= null
             then Create_Next(root,nd.children(i,j));
            end if;
          end loop;
        end loop;
        nd.roco := -1;
      end if;
    end Create_Next;

  begin
    Create_Next(root,root);
  end Top_Bottom_Create;

  procedure Q_Top_Bottom_Create ( root : in Link_to_Node;
                                  n,lag : in natural32 ) is

    procedure Create_Next ( root,nd : in Link_to_Node ) is
    begin
      if ((nd.level > 0) and (nd.roco >= 0)) then
        Q_Top_Bottom_Create1(root,nd,integer32(n),integer32(lag));
        for i in nd.children'range(1) loop
          for j in nd.children'range(2) loop
            if nd.children(i,j) /= null
             then Create_Next(root,nd.children(i,j));
            end if;
          end loop;
        end loop;
        nd.roco := -1;
      end if;
    end Create_Next;

  begin
    Create_Next(root,root);
  end Q_Top_Bottom_Create;

  procedure Top_Bottom_Create ( root : in Link_to_Node;
                                k : in Bracket; n : in natural32 ) is

    ind : constant integer32 := k'last;
    hyper : boolean;

  begin
    if ind = k'first then
      Bottom_Create1(root,root,integer32(k(k'last)),root.p);
    elsif ind > k'first then
      hyper := ((k(ind-1) = 1) and (k(ind) = 1));
      Recursive_Top_Bottom_Create
        (root,root,k,ind-1,integer32(k(ind-1)),integer32(k(ind)),
         integer32(n),1,root.p,hyper);
    end if;
  end Top_Bottom_Create;

  procedure Q_Top_Bottom_Create ( root : in Link_to_Node;
                                  k : in Bracket; n,lag : in natural32 ) is

    ind : constant integer32 := k'last;
    hyper : boolean;
    space : constant Bracket(1..integer32(lag)-root.p)
          := (1..integer32(lag)-root.p => 0);

  begin
    if ind = k'first then
      Q_Bottom_Create1(root,root,true,space,
                       integer32(k(k'last)),integer32(lag));
    elsif ind > k'first then
      hyper := ((k(ind-1) = 1) and (k(ind) = 1));
      Q_Recursive_Top_Bottom_Create
        (root,root,k,true,space,true,space,
         ind-1,integer32(k(ind-1)),integer32(k(ind)),integer32(n),
         integer32(lag),hyper);
    end if;
  end Q_Top_Bottom_Create;

  function Create_Leveled_Poset ( root : Link_to_Node )
                                return Array_of_Nodes is

    res : Array_of_Nodes(0..root.level);

  begin
    for i in res'range loop
      res(i) := Find_Node(root,i);
    end loop;
    return res;
  end Create_Leveled_Poset;

  function Create_Indexed_Poset ( poset : Array_of_Nodes )
                                return Array_of_Array_of_Nodes is

    res : Array_of_Array_of_Nodes(poset'range);
    ptr : Link_to_Node;

  begin
    for i in poset'range loop
      if poset(i) /= null then
        res(i) := new Array_of_Nodes
                        (1..integer32(Number_of_Siblings(poset(i))));
        ptr := poset(i);
        for j in res(i)'range loop
          res(i)(j) := ptr;
          res(i)(j).label := j;
          res(i)(j).child_labels := Labels_of_Children(res,ptr.all);
          ptr := ptr.next_sibling;
        end loop;
      end if;
    end loop;
    return res;
  end Create_Indexed_Poset;

-- SELECTORS :

  function Equal ( nd1,nd2 : Node ) return boolean is
  begin
    if nd1.level /= nd2.level then
      return false;
    elsif not Equal(nd1.top,nd2.top) then
      return false;
    else
      return Equal(nd1.bottom,nd2.bottom);
    end if;
  end Equal;

  function Is_Leaf ( nd : Node ) return boolean is
  begin
    for i in nd.children'range(1) loop
      for j in nd.children'range(2) loop
        if nd.children(i,j) /= null
         then return false;
        end if;
      end loop;
    end loop;
    return true;
  end Is_Leaf;

  function Find_Node ( root : Link_to_Node; lvl : integer32 )
                     return Link_to_Node is

    res,fst : Link_to_Node := null;

    procedure Search_First ( current : in Link_to_Node ) is

    -- DESCRIPTION :
    --   Scans the list of previous siblings and sets fst to the node
    --   that does not have any previous siblings.

    -- REQUIRED : current /= null.

    begin
      if current.prev_sibling = null
       then fst := current;
       else Search_First(current.prev_sibling);
      end if;
    end Search_First;

  begin
    if root.level = lvl then
      res := root;
    elsif root.level > lvl then
      for i in root.children'range(1) loop
        for j in root.children'range(2) loop
          if root.children(i,j) /= null
           then res := Find_Node(root.children(i,j),lvl);
          end if;
          exit when (res /= null);
        end loop;
        exit when (res /= null);
      end loop;
    end if;
    if res = null
     then fst := res;
     else Search_First(res);
    end if;
    return fst;
  end Find_Node;

  function Number_of_Siblings ( nd : Link_to_Node ) return natural32 is
  begin
    if nd = null
     then return 0;
     else return 1 + Number_of_Siblings(nd.next_sibling);
    end if;
  end Number_of_Siblings;

  function Number_of_Children ( nd : Node ) return natural32 is

    cnt : natural32 := 0;

  begin
    for i in nd.children'range(1) loop
      for j in nd.children'range(2) loop
        if nd.children(i,j) /= null
         then cnt := cnt + 1;
        end if;
      end loop;
    end loop;
    return cnt;
  end Number_of_Children;

-- ITERATORS :

  procedure Enumerate_Siblings ( nd : in Node ) is

    cont : boolean := true;

  begin
    Report(nd,cont);
    if cont and nd.next_sibling /= null
     then Enumerate_Siblings(nd.next_sibling.all);
    end if;
  end Enumerate_Siblings;

  procedure Enumerate_Grand_Children ( nd : in Node; k : in natural32 ) is

    cont : boolean := true;

    procedure Enumerate_Children ( current : in node; l : in natural32 ) is
    begin
      for i in current.children'range(1) loop
        for j in current.children'range(1) loop
          if current.children(i,j) /= null then
            if l = 1
             then Report(current.children(i,j),cont);
             else Enumerate_Children(current.children(i,j).all,l-1);
            end if;
          end if;
          exit when not cont;
        end loop;
        exit when not cont;
      end loop;
    end Enumerate_Children;

  begin
    Enumerate_Children(nd,k);
  end Enumerate_Grand_Children;

  procedure Modify_Siblings ( nd : in out Node ) is

    cont : boolean := true;

  begin
    Modify(nd,cont);
    if cont and nd.next_sibling /= null
     then Modify_Siblings(nd.next_sibling.all);
    end if;
  end Modify_Siblings;

-- COMBINATORIAL ROOT COUNTING :

  procedure Count_Roots ( poset : in out Array_of_Nodes ) is

    procedure Initialize ( nd : in out Node; continue : out boolean ) is
    begin
      nd.roco := 1;
      continue := true;
    end Initialize;
    procedure Initialize_Leaves is new Modify_Siblings(Initialize);

    procedure Add_Children ( nd : in out Node; continue : out boolean ) is
    begin
      nd.roco := 0;
      for i in nd.children'range(1) loop
        for j in nd.children'range(2) loop
          if nd.children(i,j) /= null
           then nd.roco := nd.roco + nd.children(i,j).roco;
          end if;
        end loop;
      end loop;
      continue := true;
    end Add_Children;
    procedure Add_Children_Counts is new Modify_Siblings(Add_Children);

  begin
    if poset(0) /= null
     then Initialize_Leaves(poset(0).all);
    end if;
    for i in 1..poset'last loop
      if poset(i) /= null
       then Add_Children_Counts(poset(i).all);
      end if;
    end loop;
  end Count_Roots;

  function Row_Root_Count_Sum
             ( poset : Array_of_Nodes; i : natural32 ) return natural32 is

    res : natural32 := 0;

    procedure Count ( lnd : in Link_to_Node ) is
    begin
      if lnd /= null then
        res := res + natural32(lnd.roco);
        Count(lnd.next_sibling);
      end if;
    end Count;

  begin
    Count(poset(integer32(i)));
    return res;
  end Row_Root_Count_Sum;

  function Root_Count_Sum ( poset : Array_of_Nodes ) return natural32 is

    res : natural32 := 0;

  begin
    for i in 1..poset'last loop
      res := res + Row_Root_Count_Sum(poset,natural32(i));
    end loop;
    return res;
  end Root_Count_Sum;

-- DESTRUCTORS :

  procedure free is new unchecked_deallocation(Node,Link_to_Node);

  procedure Clear ( nd : in out Node ) is
  begin
    if nd.next_sibling /= null
     then Clear(nd.next_sibling);
    end if;
  end Clear;

  procedure Clear ( lnd : in out Link_to_Node ) is
  begin
    if lnd /= null
     then Clear(lnd.all); free(lnd);
    end if;
  end Clear;

  procedure Clear ( arrnd : in out Array_of_Nodes ) is
  begin
    for i in arrnd'range loop
      Clear(arrnd(i));
    end loop;
  end Clear;

  procedure Clear ( arrnd : in out Link_to_Array_of_Nodes ) is

    procedure free is
      new unchecked_deallocation(Array_of_Nodes,Link_to_Array_of_Nodes);

  begin
    if arrnd /= null
     then Clear(arrnd.all); free(arrnd);
    end if;
  end Clear;

  procedure Clear ( arrnd : in out Array_of_Array_of_Nodes ) is
  begin
    for i in arrnd'range loop
      Clear(arrnd(i));
    end loop;
  end Clear;

  procedure Clear ( matnd : in out Matrix_of_Nodes ) is
  begin
    for i in matnd'range(1) loop
      for j in matnd'range(2) loop
        if matnd(i,j) /= null
         then free(matnd(i,j));
        end if;
      end loop;
    end loop;
  end Clear;

end Localization_Posets;
