with unchecked_deallocation;

-- for testing only:
with text_io; use text_io;
with Standard_Integer_Numbers_io; use Standard_Integer_Numbers_io;

package body Index_Tree_LP is

  procedure LPdata_Init
               ( p : in out Link_to_LPdata; n : in integer32;
	         J : in Standard_Integer_Vectors.Link_to_Vector;
	         x : in Standard_Floating_Vectors.Link_to_Vector;
	         A : in Standard_Floating_Matrices.Link_to_Matrix ) is
  begin
    p.dim := n;
    p.JJJ := new Standard_Integer_Vectors.Vector(0..n-1);
    p.xxx := new Standard_Floating_Vectors.Vector(0..n-1);
    p.INV := new Standard_Floating_Matrices.Matrix(0..n-1,0..n-1);
    for i in 0..n-1 loop
      p.JJJ(i) := J(i);
      p.xxx(i) := x(i);
      for j in 0..n-1 loop
        p.INV(i,j) := A(i,j);
      end loop;
    end loop;
  end LPdata_Init;

  procedure IndexNode_Init
               ( p : in out Link_to_IndexNode; i : in integer32;
                 ptr : in Link_to_LPdata ) is
  begin
    p.idx := i;
    p.info := ptr;
    p.S := null;
  end IndexNode_Init;

  procedure LPPL_Init
               ( p : in out Link_to_LPPL; A : in Link_to_LPdata;
                 N : in Link_to_LPPL ) is
  begin
    p.addr := A;
    p.next := N;
  end LPPL_Init;

  procedure IT_LP_Init
               ( p : in out Link_to_IT_LP; nSpt : in integer32;
                 mtype : in Standard_Integer_Vectors.Link_to_Vector ) is

    itmp,sum : integer32 := 0;

  begin
    for i in 0..(nSpt-1) loop
      itmp := itmp + mtype(i);
    end loop;
    p.MaxLevels := itmp + nSpt + 1;       -- "+1" is for unused 0 level
    p.CurLevel := 1;
    p.DIM := new Standard_Integer_Vectors.Vector'(0..(p.MaxLevels-1) => 0);
    p.NP := new Standard_Integer_Vectors.Vector'(0..(p.MaxLevels-1) => 0);
    p.cell := new Standard_Integer_Vectors.Vector'(0..(p.MaxLevels-1) => 0);
    p.InSpt := new Standard_Integer_Vectors.Vector'(0..(p.MaxLevels-1) => 0);
    p.minNP := new Standard_Integer_Vectors.Vector'(0..(p.MaxLevels-1) => 0);
    p.LP := new Array_of_LPPL(0..(p.MaxLevels-1));
    for i in 0..(p.MaxLevels-1) loop
      p.LP(i) := new LPPL;
      LPPL_Init(p.LP(i),null,null);
    end loop;
    p.IT := new Array_of_IndexNodes'(0..(p.MaxLevels-1) => null);
    p.last := new Array_of_IndexNodes'(0..(p.MaxLevels-1) => null);
    p.DIM(sum) := itmp;
    itmp := itmp + 1;
    for i in 0..(nSpt-1) loop
      p.minNP(sum) := mtype(i) + 1;
      p.InSpt(sum) := i;
      for j in 1..mtype(i) loop
        itmp := itmp - 1;
        p.DIM(sum+j) := itmp;
        p.minNP(sum+j) := mtype(i)+1-j;
      end loop;
      sum := sum + mtype(i)+1;
      if sum < p.MaxLevels
       then p.DIM(sum) := itmp;
      end if;
    end loop;
    p.IT(1) := new IndexNode;  
    IndexNode_Init(p.IT(1),-1,null);
    p.curr := p.IT(1);                            -- points to dummy node
    p.prev := p.IT(1);
    p.last(1) := p.IT(1);
    p.NP(1) := 1;
  end IT_LP_Init;

  function IT_IsEmpty ( p : Link_to_IT_LP ) return boolean is
  begin
    return (p.CurLevel < 1);
  end IT_IsEmpty;

  function IT_IsFull ( p : Link_to_IT_LP ) return boolean is
  begin
    return (p.CurLevel + 1 >= p.MaxLevels);
  end IT_IsFull;

  function IT_Level ( p : Link_to_IT_LP ) return integer32 is
  begin
    return p.CurLevel;
  end IT_Level;

  function IT_CurLPdim ( p : Link_to_IT_LP ) return integer32 is
  begin
    return p.DIM(p.CurLevel);
  end IT_CurLPdim;

  function IT_Cell ( p : Link_to_IT_LP )
                   return Standard_Integer_Vectors.Link_to_Vector is
  begin
    return p.cell;
  end IT_Cell;

  function IT_CurSptIdx ( p : Link_to_IT_LP ) return integer32 is
  begin
    return p.InSpt(p.CurLevel);
  end IT_CurSptIdx;

  function IT_MinNumPt ( p : Link_to_IT_LP ) return integer32 is
  begin
    return p.minNP(p.CurLevel);
  end IT_MinNumPt;

  function IT_NumPt ( p : Link_to_IT_LP ) return integer32 is
  begin
    return p.NP(p.CurLevel);
  end IT_NumPt;

  function IT_FixedIdxNdPtr ( p : Link_to_IT_LP ) return Link_to_IndexNode is
  begin
    if p.CurLevel >= 0
     then return p.IT(p.CurLevel);
     else return null;
    end if;
  end IT_FixedIdxNdPtr;

  procedure IT_StepBack ( p : in out Link_to_IT_LP ) is
  begin
    p.NP(p.CurLevel) := 0;
    p.CurLevel := p.CurLevel - 1;
  end IT_StepBack;

  procedure IT_Find ( p : in out Link_to_IT_LP; IDX : in integer32;
                      found : out boolean ) is
  begin
   -- put_line("entered IT_Find...");
    p.curr := p.prev.S;
   -- if p.last(p.CurLevel) = null
   --  then put_line("p.last(p.CurLevel) is null");
   --  else put_line("p.last(p.CurLevel) is not null");
   --       if p.last(p.CurLevel).S = null
   --        then put_line("  p.last(p.CurLevel).S is null");
   --        else put_line("  p.last(p.CurLevel).S is not null");
   --       end if;
   -- end if;
    while p.curr /= p.last(p.CurLevel).S loop
     -- put("IT_Find, p.curr.idx = "); put(p.curr.idx,1); new_line;
      if IDX <= p.curr.idx then
        if IDX = p.curr.idx then
          found := true;
          return; 
        else
          found := false;
          return; 
        end if;
      end if;
      p.prev := p.prev.S;
      p.curr := p.curr.S;
    end loop;
    found := false; 
  end IT_Find;

  procedure IT_ResetCurLevelTo1
               ( p : in out Link_to_IT_LP; frst : out Link_to_IndexNode ) is
  begin
   -- put("IT_Reset ... ");
    p.prev := p.IT(1);
    p.curr := p.prev.S;
    while p.curr /= null loop
      p.prev.S := p.curr.S;
      Clear(p.curr);
      p.curr := p.prev.S;
     -- put(" + ");
    end loop;
    p.NP(1) := 1;
    p.CurLevel := 1;
    frst := p.IT(1);
   -- put_line(" leaving IT_Reset.");
  end IT_ResetCurLevelTo1;

  procedure IT_RenewNP1 ( p : in out Link_to_IT_LP ) is
  begin
   -- put("IT_RenewNP1: p.NP(1) = "); put(p.NP(1),1);
    p.prev := p.IT(1);
    while p.prev.S /= null loop
      p.NP(1) := p.NP(1) + 1;
      p.prev := p.prev.S;
    end loop;
    p.last(1) := p.prev;
    p.cell(1) := p.IT(1).idx;
   -- put(" -> p.NP(1) = "); put(p.NP(1),1); new_line;
  end IT_RenewNP1;

  procedure IT_NextLevel ( p : in out Link_to_IT_LP; rcase : out integer32 ) is

    tmp : Link_to_IndexNode;

  begin
    if p.CurLevel+1 >= p.MaxLevels then
      rcase := 0;
     -- put_line("returning with rcase = 0 because p.CurLevel too high");
      return;                                                    -- case 0
    else
     -- put("p.NP("); put(p.CurLevel,1); put(") = ");
     -- put(p.NP(p.CurLevel),1); new_line;
     -- put("p.minNP("); put(p.CurLevel,1); put(") = ");
     -- put(p.minNP(p.CurLevel),1); new_line;
      if p.NP(p.CurLevel) <= p.minNP(p.CurLevel)
       then rcase := 0;                                          -- case 1
           -- put_line("returning with rcase = 0, comparing p.NP with p.minNP");
            return;
      end if;               -- now IT[CurLevel] has enough points to go on
      if p.IT(p.CurLevel+1) /= null then            -- next level nonempty
       -- put_line("IT_NextLevel, backtrack to next nonempty level...");
        tmp := p.IT(p.CurLevel+1);                         -- backtracking
               p.IT(p.CurLevel+1) := tmp.S;
        tmp.S := p.last(p.CurLevel).S;
                 p.last(p.CurLevel).S := tmp;
        tmp := p.IT(p.CurLevel).S;
               p.IT(p.CurLevel).S := tmp.S;
        tmp.S := p.IT(p.CurLevel+1);
                 p.IT(p.CurLevel+1) := tmp;
      else
       -- put_line("IT_NextLevel, next level is empty...");
        tmp := p.IT(p.CurLevel).S;
               p.IT(p.CurLevel+1) := tmp;
        p.IT(p.CurLevel).S := tmp.S;
                              tmp.S := null;
      end if;
      if p.NP(p.CurLevel) = 2
       then p.last(p.CurLevel) := p.IT(p.CurLevel);
      end if;
      p.NP(p.CurLevel) := p.NP(p.CurLevel) - 1;
      p.CurLevel := p.CurLevel + 1;
      p.NP(p.CurLevel) := p.NP(p.CurLevel) + 1;
      p.last(p.CurLevel) := p.IT(p.CurLevel);
      p.curr := p.IT(p.CurLevel);
      p.cell(p.CurLevel) := p.curr.idx;
      p.LPlast := p.LP(p.CurLevel);
      rcase := 1;
     -- put_line("returning from IT_NextLevel with rcase = 1");
      return;                                                   -- case 2
    end if;
  end IT_NextLevel;

  procedure IT_Add1
               ( p : in out Link_to_IT_LP; n : in integer32;
                 J : in Standard_Integer_Vectors.Link_to_Vector;
                 nn : in integer32;
                 JJ : in Standard_Integer_Vectors.Link_to_Vector;
                 X : in Standard_Floating_Vectors.Link_to_Vector;
                 A : in Standard_Floating_Matrices.Link_to_Matrix ) is

    AddIdx : boolean := false;
    lpp : Link_to_LPdata;
    ptr : Link_to_IndexNode;
    found : boolean;
    ind : integer32;

  begin
   -- put("entered IT_Add1, p.CurLevel = "); put(p.CurLevel,1); new_line;
    p.prev := p.IT(p.CurLevel);
   -- if p.prev = null
   --  then put_line("p.prev is null");
   --  else put_line("p.prev is not null");
   --       if p.prev.S = null
   --        then put_line("  p.prev.S is null");
   --        else put_line("  p.prev.S is not null");
   --       end if;
   -- end if;
    for i in 0..(n-1) loop    -- search for the 1st idx of J not in level 0
      IT_Find(p,J(i),found);
      if not found then                              -- add J[i] after prev
        if p.LPlast.next /= null then
          lpp := p.LPlast.next.addr;
          for k in 0..(nn-1) loop
            lpp.JJJ(k) := JJ(k);
            lpp.xxx(k) := X(k);
            for j in 0..(nn-1) loop
              lpp.INV(k,j) := A(k,j);
            end loop;
          end loop;
        else
          lpp := new LPdata;
          LPdata_Init(lpp,nn,JJ,X,A);
          p.LPlast.next := new LPPL;
          LPPL_Init(p.LPlast.next,lpp,null);
        end if;
        p.LPlast := p.LPlast.next;
        AddIdx := true;
        ind := i;
        exit;
      end if;
    end loop;
    while AddIdx loop
      if p.last(p.CurLevel).S = null then           -- IT[CurLevel] is full
        p.curr := new IndexNode;
        IndexNode_Init(p.curr,J(ind),lpp);
        p.curr.S := p.prev.S;
        p.prev.S := p.curr;
        p.prev := p.curr;
        if (p.last(p.CurLevel) = p.IT(p.CurLevel))
            or (p.last(p.CurLevel).idx < J(ind)) then
           p.last(p.CurLevel) := p.curr;
        end if;                                 -- spare slots for LP data
      elsif p.last(p.CurLevel) = p.prev then      -- after *last[CurLevel]
        ptr := p.prev.S;
        ptr.idx := J(ind);
        ptr.info := lpp;
        p.last(p.CurLevel) := ptr;
        p.prev := ptr;
      else                                        -- intermediate position
        ptr := p.last(p.CurLevel).S;
        ptr.idx := J(ind);
        ptr.info := lpp;
        p.prev.S := ptr;
        p.last(p.CurLevel).S := ptr.S;
        ptr.S := p.curr;
        p.prev := ptr;
      end if;
      p.NP(p.CurLevel) := p.NP(p.CurLevel) + 1;
      AddIdx := false;
      ind := ind + 1;
      while ind < n loop
        IT_Find(p,J(ind),found);
        if not found
         then AddIdx := true; exit;
        end if;
        ind := ind + 1;
      end loop;
    end loop;
  end IT_Add1;
                 
  procedure IT_Add2
               ( p : in out Link_to_IT_LP; oneidx,nn : in integer32;
                 JJ : in Standard_Integer_Vectors.Link_to_Vector;
                 X : in Standard_Floating_Vectors.Link_to_Vector;
                 A : in Standard_Floating_Matrices.Link_to_Matrix ) is

    lpp : Link_to_LPdata;
    ptr : Link_to_IndexNode;
    found : boolean;

  begin
   -- put_line("entered IT_Add2 ...");
    p.prev := p.IT(p.CurLevel); 
    IT_Find(p,oneidx,found);
    if not found then                              -- add oneidx after prev
      if p.LPlast.next /= null then
        lpp := p.LPlast.next.addr;
        for i in 0..(nn-1) loop
          lpp.JJJ(i) := JJ(i);
          lpp.xxx(i) := X(i);
          for j in 0..(nn-1) loop
            lpp.INV(i,j) := A(i,j);
          end loop;
        end loop;
      else
        lpp := new LPdata;
        LPdata_Init(lpp,nn,JJ,X,A);
        p.LPlast.next := new LPPL; 
        LPPL_Init(p.LPlast.next,lpp,null);
      end if;
      p.LPlast := p.LPlast.next;
      if p.last(p.CurLevel).S = null then            -- IT[CurLevel] is full
        p.curr := new IndexNode;
        IndexNode_Init(p.curr,oneidx,lpp);
        p.curr.S := p.prev.S;
        p.prev.S := p.curr;
        p.prev := p.curr;
        if p.last(p.CurLevel).idx < oneidx then
          p.last(p.CurLevel) := p.curr;
        end if;                             -- have spare slots for LP data
      elsif p.last(p.CurLevel) = p.prev then       -- after *last[CurLevel]
        ptr := p.prev.S;
        ptr.idx := oneidx;
        ptr.info := lpp;
        p.last(p.CurLevel) := ptr;
        p.prev := ptr;
      else       --  intermediate position before last, though last->S != 0
        ptr := p.last(p.CurLevel).S;
        ptr.idx := oneidx;
        ptr.info := lpp;
        p.prev.S := ptr;
        p.last(p.CurLevel).S := ptr.S;
        ptr.S := p.curr;
        p.prev := ptr;
      end if;
      p.NP(p.CurLevel) := p.NP(p.CurLevel) + 1;
    end if;
  end IT_Add2;

  procedure IT_LP_DEL ( p : in out Link_to_IT_LP ) is
  begin
    IT_FreeIT(p);
    IT_FreeLP(p); 
    Standard_Integer_Vectors.Clear(p.DIM);
    Standard_Integer_Vectors.Clear(p.NP);
    Standard_Integer_Vectors.Clear(p.cell);
    Standard_Integer_Vectors.Clear(p.InSpt);
    Standard_Integer_Vectors.Clear(p.minNP);
    Clear(p.last);
    Clear(p.LP);
    Clear(p.IT);
  end IT_LP_DEL;
 
  procedure IT_FreeIT ( p : in out Link_to_IT_LP ) is
  begin
    p.CurLevel := p.MaxLevels-1;
    while p.CurLevel > 1 loop
      p.prev := p.IT(p.CurLevel);
      p.curr := p.prev.S;
      while p.curr /= null loop
         p.prev.S := p.curr.S;
         Clear(p.curr);
         p.curr := p.prev.S;
      end loop;
      p.CurLevel := p.CurLevel - 1;
    end loop;
    p.CurLevel := 0;
    while p.CurLevel < p.MaxLevels loop
      Clear(p.IT(p.CurLevel));
      p.CurLevel := p.CurLevel + 1;
    end loop;
  end IT_FreeIT;

  procedure IT_FreeLP ( p : in out Link_to_IT_LP ) is

    lpp : Link_to_LPdata;

  begin
    p.CurLevel := p.MaxLevels-1;
    while p.CurLevel > 1 loop
      p.LPlast := p.LP(p.CurLevel).next;
      while p.LPlast /= null loop
        p.LP(p.CurLevel).next := p.LPlast.next;
        lpp := p.LPlast.addr;
        if lpp /= null then
          Standard_Integer_Vectors.Clear(lpp.JJJ);
          Standard_Floating_Vectors.Clear(lpp.xxx);
          Standard_Floating_Matrices.Clear(lpp.INV);
        end if;
        Clear(p.LPlast);
        p.LPlast := p.LP(p.CurLevel).next;
      end loop;
      p.CurLevel := p.CurLevel - 1;
    end loop;
    p.CurLevel := 0;
    while p.CurLevel < p.MaxLevels loop
      Clear(p.LP(p.CurLevel));
      p.CurLevel := p.CurLevel + 1;
    end loop;
  end IT_FreeLP;

  procedure Clear ( p : in out Link_to_IndexNode ) is

    procedure free is
      new unchecked_deallocation(IndexNode,Link_to_IndexNode);

  begin
    free(p);
    p := null;
  end Clear;

  procedure Clear ( p : in out Link_to_LPdata ) is

    procedure free is new unchecked_deallocation(LPdata,Link_to_LPdata);

  begin
    free(p);
    p := null;
  end Clear;

  procedure Clear ( p : in out Link_to_LPPL ) is

    procedure free is new unchecked_deallocation(LPPL,Link_to_LPPL);

  begin
    free(p);
    p := null;
  end Clear;

  procedure Clear ( p : in out Array_of_LPPL ) is
  begin
    for i in p'range loop
      Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Array_of_LPPL ) is

    procedure free is
      new unchecked_deallocation(Array_of_LPPL,Link_to_Array_of_LPPL);

  begin
    if p /= null
     then Clear(p.all); free(p); p := null;
    end if;
  end Clear;

  procedure Clear ( p : in out Array_of_IndexNodes ) is
  begin
    for i in p'range loop
      Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Array_of_IndexNodes ) is

    procedure free is
      new unchecked_deallocation(Array_of_IndexNodes,
                                 Link_to_Array_of_IndexNodes);

  begin
    if p /= null
     then Clear(p.all); free(p); p := null;
    end if;
  end Clear;
 
end Index_Tree_LP;
