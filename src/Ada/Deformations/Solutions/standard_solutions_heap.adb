with unchecked_deallocation;
with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Complex_Norms_Equals;

package body Standard_Solutions_Heap is

-- WEIGHING SOLUTIONS :

  function Random_Weight_Vector
             ( nbr : in integer32 )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..nbr);

  begin
    for k in 1..nbr loop
      res(k) := Standard_Random_Numbers.Random;
    end loop;
    return res;
  end Random_Weight_Vector;

  function Weight ( v : Standard_Complex_Vectors.Vector;
                    w : Standard_Floating_Vectors.Vector )
                  return double_float is

    res : double_float := 0.0;
    vdx : integer32 := v'first;
    wdx : integer32 := w'first;
    val : double_float;

  begin
    while wdx <= w'last loop
      val := Standard_Complex_Numbers.REAL_PART(v(vdx));
      res := res + w(wdx)*val;
      wdx := wdx + 1;
      exit when (wdx > w'last);
      val := Standard_Complex_Numbers.IMAG_PART(v(vdx));
      res := res + w(wdx)*val;
      wdx := wdx + 1;
      exit when (wdx > w'last);
      vdx := vdx + 1;
      exit when (vdx > v'last);
    end loop;
    return res;
  end Weight;

-- OPERATIONS ON BUCKET VECTORS :

  function Size ( b : Bucket_Vector ) return integer32 is
  begin
    if Is_Null(b)
     then return 0;
     else return integer32(Length_Of(b))*bucketsize;
    end if;
  end Size;

  function New_Heap_Items return Link_to_Heap_Items is

    res : Link_to_Heap_Items;
    hit : Heap_Items(0..bucketsize-1);

  begin
    for k in hit'range loop
      hit(k).first := 0.0;
      hit(k).second := 0.0;
      hit(k).idx := 0;
      hit(k).ls := null;
    end loop;
    res := new Heap_Items'(hit);
    return res;
  end New_Heap_Items;

  procedure Assign ( b,b_last : in out Bucket_Vector;
                     i : in integer32; v : in Heap_Item ) is

    lhi : Link_to_Heap_Items;
    done : boolean := false;
    idx : integer32 := i;
    tmp : Bucket_Vector;

  begin
    if Is_Null(b) then
      declare
        newlhi : constant Link_to_Heap_Items := New_Heap_Items;
      begin
        if idx <= newlhi'last
         then newlhi(idx) := v; done := true;
        end if;
        Construct(newlhi,b); b_last := b;
      end;
      while not done loop
        idx := idx - bucketsize;
        declare
          newlhi : constant Link_to_Heap_Items := New_Heap_Items;
        begin
          if idx <= newlhi'last
           then newlhi(idx) := v; done := true;
          end if;
          Append(b,b_last,newlhi);
        end;
      end loop;
    else
      tmp := b;
      while not Is_Null(tmp) loop
        lhi := Head_Of(tmp);
        if idx <= lhi'last
         then lhi(idx) := v; done := true; exit;
         else tmp := Tail_Of(tmp); idx := idx - bucketsize;
        end if;
      end loop;
      while not done loop
        declare
          newlhi : constant Link_to_Heap_Items := New_Heap_Items;
        begin
          if idx <= newlhi'last
           then newlhi(idx) := v; done := true;
          end if;
          Append(b,b_last,newlhi);
        end;
        if not done
         then idx := idx - bucketsize;
        end if;
      end loop;
    end if;
  end Assign;

  function Retrieve ( b : Bucket_Vector; i : integer32 ) return Heap_Item is

    res : Heap_Item;
    tmp : Bucket_Vector := b;
    lhi : Link_to_Heap_Items;
    idx : integer32 := i;

  begin
    res.first := 0.0;
    res.second := 0.0;
    res.idx := -1;
    res.ls := null;
    while not Is_Null(tmp) loop
      lhi := Head_Of(tmp);
      if idx <= lhi'last then
        res := lhi(idx);
        return res;
      else
        tmp := Tail_Of(tmp);
        idx := idx - bucketsize;
      end if;
    end loop;
    return res;
  end Retrieve;

  procedure Clear ( h : in out Link_to_Heap_Items ) is

    procedure free is
      new unchecked_deallocation(Heap_Items,Link_to_Heap_Items);

  begin
    if h /= null
     then free(h);
    end if;
  end Clear;

  procedure Clear ( b : in out Bucket_Vector ) is

    tmp : Bucket_Vector := b;
    lhi : Link_to_Heap_Items;

  begin
    while not Is_Null(tmp) loop
      lhi := Head_Of(tmp);
      Clear(lhi);
      tmp := Tail_Of(tmp);
    end loop;
    Buckets.Clear(Buckets.List(b));
  end Clear;

-- HEAP OPERATIONS :

  procedure Swap_from_Bottom ( h : in out Heap; p : in integer32 ) is

    parent : integer32;
    parent_hi,p_hi : Heap_Item;

  begin
    if p > 0 then
      parent := (p-1)/2;
     -- if h.values(parent).wgt < h.values(p).wgt then
      p_hi := Retrieve(h.values,p);
      parent_hi := Retrieve(h.values,parent);
      if parent_hi.first < p_hi.first then
       -- declare
       --   tmp : constant Heap_Item := h.values(p);
       -- begin
       --   h.values(p) := h.values(parent);
       --   h.values(parent) := tmp;
        Assign(h.values,h.values_last,p,parent_hi);
        Assign(h.values,h.values_last,parent,p_hi);
        Swap_from_Bottom(h,parent);
       -- end;
      end if;
    end if;
  end Swap_from_Bottom;

  procedure Push ( h : in out Heap; i : in Heap_Item ) is
  begin
    if h.bottom = -1 then
      h.bottom := 0;
     -- h.values(0) := i;
      Assign(h.values,h.values_last,0,i);
      h.bottom := 1;
    else -- if h.bottom < h.nbr then
     -- h.values(h.bottom) := i;
      Assign(h.values,h.values_last,h.bottom,i);
      Swap_from_Bottom(h,h.bottom);
      h.bottom := h.bottom + 1;
    end if;
  end Push;

  function Max_Child ( h : Heap; p : integer32 ) return integer32 is

    left,right : integer32;
    left_hi,right_hi : Heap_Item;

  begin
    if h.bottom <= p then
      return -1;
    else
      left := 2*p+1;
      if left >= h.bottom then
        return -1;
      else
        right := 2*p+2;
        if right >= h.bottom then
          return left;
        else
         -- if h.values(left).wgt > h.values(right).wgt
          left_hi := Retrieve(h.values,left);
          right_hi := Retrieve(h.values,right);
          if left_hi.first > right_hi.first
           then return left;
           else return right;
          end if;
        end if;
      end if;
    end if;
  end Max_Child;

  procedure Swap_from_Top ( h : in out Heap; p : in integer32 ) is

    child : integer32;
    child_hi,p_hi : Heap_Item;

  begin
    if h.bottom > 0 then
      child := Max_Child(h,p);
      if child /= -1 then
       -- if h.values(child).wgt > h.values(p).wgt then
        child_hi := Retrieve(h.values,child);
        p_hi := Retrieve(h.values,p);
        if child_hi.first > p_hi.first then
         -- declare
         --   tmp : constant Heap_Item := h.values(p);
         -- begin
         --   h.values(p) := h.values(child);
         --   h.values(child) := tmp;
         -- end;
          Assign(h.values,h.values_last,p,child_hi);
          Assign(h.values,h.values_last,child,p_hi);
          Swap_from_Top(h,child);
        end if;
      end if;
    end if;
  end Swap_from_Top;

  procedure Pop ( h : in out Heap ) is

    bottom_hi : Heap_Item;

  begin
    if h.bottom < 1 then -- make empty heap
      h.bottom := -1;
    else
      h.bottom := h.bottom-1;
     -- h.values(0) := h.values(h.bottom);
      bottom_hi := Retrieve(h.values,h.bottom);
      Assign(h.values,h.values_last,0,bottom_hi);
      Swap_from_Top(h,0);
    end if;
  end Pop;

  procedure Push ( h : in out Heap;
                   w1,w2 : in double_float; i : in integer32;
                   ls : in Link_to_Solution ) is

    item : Heap_Item;

  begin
    item.first := w1;
    item.second := w2;
    item.idx := i;
    item.ls := ls;
    Push(h,item);
  end Push;

  procedure Clear ( h : in out Heap ) is
  begin
    Clear(h.values);
  end Clear;

-- CLUSTER REPORT :

  procedure Count_Clusters
              ( h : in out Heap; tol : in double_float;
                cnt : out natural32; verbose : in boolean := true ) is

    prevtop,top : Heap_Item; 
    isclus : boolean;

    use Standard_Complex_Norms_Equals; -- for the Equal function

  begin
    cnt := 0;
    if h.bottom > 0 then
     -- prevtop := h.values(0);
      prevtop := Retrieve(h.values,0);
      while h.bottom > 0 loop
        Pop(h);
        if h.bottom > 0 then
         -- top := h.values(0);
          top := Retrieve(h.values,0);
          if prevtop.first - top.first > tol then
            if verbose then
              put(prevtop.first); put(" > "); put(top.first);
              put_line("  ok");
            end if;
          else
            if verbose then
              put(prevtop.first); put(" ~ "); put(top.first);
              put("  second check...");
            end if;
            if prevtop.second - top.second > tol then
              if verbose then
                put("  ");
                put(prevtop.second); put(" > "); put(top.second);
                put_line("  ok");
              end if;
            else
              isclus := Equal(prevtop.ls.v,top.ls.v,tol);
              if isclus 
               then cnt := cnt + 1;
              end if;
              if verbose then
                put(prevtop.second); put(" ~ "); put(top.second);
                if isclus then
                  put("  "); put(prevtop.idx,1); put(" clusters with ");
                  put(top.idx,1); new_line;
                else
                  put_line("  not clustered");
                end if;
              end if;
            end if;
          end if;
          prevtop := top;
        end if;
      end loop;
    end if;
  end Count_Clusters;

end Standard_Solutions_Heap;
