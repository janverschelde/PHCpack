with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;
with generic_lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with Standard_Solution_Diagnostics;
with Standard_Condition_Tables;
with Standard_Coefficient_Circuits;      use Standard_Coefficient_Circuits;
with Standard_Circuit_Makers;            use Standard_Circuit_Makers;
with Standard_Newton_Circuits;           use Standard_Newton_Circuits;

procedure ts_newcirc is

-- DESCRIPTION :
--   Test the development of Newton's method on coefficient circuits.

  function Random_Weight_Vector
             ( nbr : in integer32 )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of nbr random floating-point numbers,
  --   in the interval [-1, +1].

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

  -- DESCRIPTION :
  --   Returns the value of the vector v for the weights w,
  --   obtained by multiplying the real and imaginary parts
  --   of the components of v with the weights in w.

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

  type Heap_Item is record
    first,second : double_float;
    idx : integer32;
    ls : Link_to_Solution;
  end record;

  type Heap_Items is array ( integer32 range <> ) of Heap_Item;
  type Link_to_Heap_Items is access Heap_Items;

  package Buckets is new Generic_Lists(Link_to_Heap_Items);
  type Bucket_Vector is new Buckets.List;

  bucketsize : constant integer32 := 1024;

  function Size ( b : Bucket_Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns the number of elements in the bucket vector b.

  begin
    if Is_Null(b)
     then return 0;
     else return integer32(Length_Of(b))*bucketsize;
    end if;
  end Size;

  function New_Heap_Items return Link_to_Heap_Items is

  -- DESCRIPTION :
  --   Returns an allocated vector of heap items,
  --   with everything initialized to zero.

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

  -- DESCRIPTION :
  --   Assigns to b(i) the value v.
  --   The b_last points to the last array in b.

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

  -- DESCRIPTION :
  --   Returns the item at position i in the vector b.
  --   If i exceeds the size of b, then the idx of the returned item is -1.

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

 -- type Heap ( nbr : integer32 ) is record
 --   bottom : integer32 := -1; -- next free spot on heap
 --   values : Heap_Items(0..nbr);
 -- end record;

  type Heap is record
    bottom : integer32 := -1; -- next free spot on heap
    values,values_last : Bucket_Vector;
  end record;

  procedure Swap_from_Bottom ( h : in out Heap; p : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps items of the heap h, starting at position p,
  --   as long as the item at position p has a larger value
  --   than the value of its parent.

  -- REQUIRED : p is a valid index for h.values,
  --   initially the position of an item pushed into the heap h.

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

  -- DESCRIPTION :
  --   Pushes the item i into the heap h.

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

  -- DESCRIPTION :
  --   Returns -1 or the index of the largest child of the node
  --   at position p in the heap h.
  --   The value of h.bottom is the position of the empty spot
  --   after the last valid item on the heap.

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

  -- DESCRIPTION :
  --   Swaps items of the heap h, starting at position p,
  --   as long as the item at position p has a smaller value
  --   than the value of any of its children.

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

  -- DESCRIPTION :
  --   Removes the top item at h.values(0).

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

  -- DESCRIPTION :
  --   Makes an item with weight w and index i
  --   and pushes that item into the heap h.

    item : Heap_Item;

  begin
    item.first := w1;
    item.second := w2;
    item.idx := i;
    item.ls := ls;
    Push(h,item);
  end Push;

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

  procedure LU_Newton_Steps
              ( s : in Link_to_System;
                v : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton steps on the system s starting at the vector v.
 
  -- REQUIRED : The system is square.

    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    info : integer32;
    vxr : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    res,rco,err : double_float;
    ans : character;
    condition : boolean;

  begin
    put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
    condition := (ans = 'y');
    put_line("The vector v :"); put_line(v);
    loop
      if condition
       then LU_Newton_Step(s,v,xr,xi,ipvt,res,rco,err,true);
       else LU_Newton_Step(s,v,xr,xi,ipvt,info,res,err,true);
      end if;
      put_line("The vector v :"); put_line(v);
      put("Another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end LU_Newton_Steps;

  procedure Interactive_Run
              ( s : in Link_to_System; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Interactive application of Newton's method on the system of
  --   coefficient circuits s, starting at the solutions in sols.
  --   Pauses frequently, prompting user to continue.

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    ans : character;
    cnt : integer32 := 0;

  begin
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      put("Running Newton's method on solution ");
      put(cnt,1); put_line(" ...");
      LU_Newton_Steps(s,ls.v);
      put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      ptr := Tail_Of(ptr);
    end loop;
  end Interactive_Run;

  procedure Show_Parameters
              ( maxit : in natural32;
                tolres,tolerr,tolsing : in double_float;
                condition : in boolean ) is

  -- DESCRIPTION :
  --   Displays the values of the parameters.

  -- ON ENTRY :
  --   maxit        maximum number of iterations;
  --   tolres       tolerance on the residual;
  --   tolerr       tolerance on the forward error;
  --   tolsing      tolerance on a singularity;
  --   condition    true if condition number is wanted,
  --                false otherwise.

  begin
    put_line("Parameter Settings :");
    put("  1. maximum number of iterations : "); put(maxit,1); new_line;
    put("  2. tolerance on residual        :"); put(tolres,3); new_line;
    put("  3. tolerance on forward error   :"); put(tolerr,3); new_line;
    put("  4. tolerance on singularity     :"); put(tolsing,3); new_line;
    put("  5. condition number wanted      : ");
    if condition
     then put_line("yes");
     else put_line("no");
    end if;
  end Show_Parameters;

  procedure Set_Parameters
              ( maxit : out natural32;
                tolres,tolerr,tolsing : out double_float;
                condition : out boolean ) is

  -- DESCRIPTION :
  --   Sets the parameters to run several steps with Newton's method.

  -- ON RETURN :
  --   maxit        maximum number of iterations;
  --   tolres       tolerance on the residual;
  --   tolerr       tolerance on the forward error;
  --   tolsing      tolerance on a singularity;
  --   condition    true if condition number is wanted,
  --                false otherwise.

    ans : character;

  begin
    maxit := 4; condition := true;
    tolres := 1.0E-8; tolerr := 1.0E-8; tolsing := 1.0E-8;
    loop
      Show_Parameters(maxit,tolres,tolerr,tolsing,condition);
      put("Type 1, 2, 3, 4, or 5 to set parameter, or 0 to exit : ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      case ans is
        when '1' => put("-> maximum number of iterations : "); get(maxit);
        when '2' => put("-> tolerance on residual :"); get(tolres);
        when '3' => put("-> tolerance on forward error :"); get(tolerr);
        when '4' => put("-> tolerance on singularity :"); get(tolsing);
        when '5' => put("-> condition number wanted ? (y/n) ");
                    Ask_Yes_or_No(ans); condition := (ans = 'y');
        when others => null;
      end case;
    end loop;
  end Set_Parameters;

  procedure Monitor_Report
              ( idx : in integer32; fail,isreal : in boolean;
                err,rco,res,wgt,tolsing : in double_float ) is

  -- DESCRIPTION :
  --   Writes one line to screen with a report on a solution,
  --   to monitor the progress on the verification.

  -- ON ENTRY :
  --   idx      index number of the current solution;
  --   fail     true if Newton's method failed, false otherwise;
  --   isreal   true if real, false otherwise (only if not fail);
  --   err      forward error;
  --   rco      estimate for the inverse condition number;
  --   res      residual;
  --   wgt      weight of the coordinates;
  --   tolsing  tolerance to decide whether a solution is singular.

  begin
    put(idx,1); put(" : ");
    if fail then
      put_line("no solution");
    else
      put("err :"); put(err,2);
      put("  rco :"); put(rco,2);
      put("  res :"); put(res,2);
      put("  wgt :"); put(wgt,2);
      if isreal
       then put(" real");
       else put(" complex");
      end if;
      if rco < tolsing
       then put_line(" singular");
       else put_line(" regular");
      end if;
    end if;
  end Monitor_Report;

  procedure Monitored_Run
              ( s : in Link_to_System; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Runs several steps of Newton's method on the system s,
  --   starting at the solutions in sols.
  --   For each solution writes one line to screen.

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    info : integer32;
    vxr : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    startres,res,rco,err,tolres,tolerr,tolsing : double_float;
    numit,maxit : natural32 := 0;
    fail,condition : boolean;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 

  begin
    new_line;
    Set_Parameters(maxit,tolres,tolerr,tolsing,condition);
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1; put(cnt,1); put(" : ");
      if condition then
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        startres,res,rco,err,numit,fail,false);
        put("  err :"); put(err,3); put("  rco :"); put(rco,3);
        put("  res :"); put(res,3); put("  #steps : "); put(numit);
        if fail
         then put_line("  failure");
         else put_line("  success");
        end if;
        Standard_Condition_Tables.Update_Corrector(t_err,err);
        Standard_Condition_Tables.Update_Condition(t_rco,rco);
        Standard_Condition_Tables.Update_Residuals(t_res,res);
      else
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        info,startres,res,err,numit,fail,false);
        put("  err :"); put(err,3);
        put("  res :"); put(res,3); put("  #steps : "); put(numit);
        if fail
         then put_line("  failure");
         else put_line("  success");
        end if;
        Standard_Condition_Tables.Update_Corrector(t_err,err);
        Standard_Condition_Tables.Update_Residuals(t_res,res);
      end if;
      ptr := Tail_Of(ptr);
    end loop;
    Standard_Condition_Tables.Write_Tables(standard_output,t_err,t_res,t_rco);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end Monitored_Run;

  procedure Monitored_Run
              ( file : in file_type;
                s : in Link_to_System; sols : in Solution_List;
               -- len : in integer32;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs several steps of Newton's method on the system s,
  --   starting at the solutions in sols.
  --   Writes the output to file.
  --   If verbose, then one line is written for every solution.

    timer : Timing_Widget;
    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    info : integer32;
    vxr : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    startres,res,rco,err,tolres,tolerr,tolsing : double_float;
    cntfail,cntreal,cntcmplx,cntregu,cntsing,cntclus : natural32 := 0;
    numit,maxit : natural32 := 0;
    fail,condition,isreal : boolean;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 
    nbr : constant integer32 := 2*s.dim;
    wv1 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Random_Weight_Vector(nbr);
    wv2 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Random_Weight_Vector(nbr);
    val1,val2 : double_float;
    weights : Heap; -- Heap(len);

  begin
    new_line;
    Set_Parameters(maxit,tolres,tolerr,tolsing,condition);
    new_line;
    weights.bottom := -1; -- make sure heap is declared as empty
    new_line;
    put_line("See the output file for results ...");
    new_line;
    tstart(timer);
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      put(file,"Solution "); put(file,cnt,1);
      put(file," :    start residual :");
      if condition then
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        startres,res,rco,err,numit,fail,false);
        put(file,startres,3);
        put(file,"  #iterations : "); put(file,numit,1);
        if fail
         then put_line(file,"  failure");
         else put_line(file,"  success");
        end if;
        Standard_Complex_Solutions_io.put_vector(file,ls.v);
        put(file,"== err :"); put(file,err,3);
        put(file," = rco :"); put(file,rco,3);
        put(file," = res :"); put(file,res,3);
        if fail then
          put_line(file," == no solution"); cntfail := cntfail + 1;
        else
          isreal := Standard_Solution_Diagnostics.Is_Real(ls.all,tolsing);
          val1 := Weight(ls.v,wv1);
          val2 := Weight(ls.v,wv2);
          Push(weights,val1,val2,cnt,ls);
          if isreal
           then put(file," == real");    cntreal := cntreal + 1;
           else put(file," == complex"); cntcmplx := cntcmplx + 1;
          end if;
          if rco < tolsing
           then put_line(file," singular"); cntsing := cntsing + 1;
           else put_line(file," regular");  cntregu := cntregu + 1;
          end if;
        end if;
        if verbose
         then Monitor_Report(cnt,fail,isreal,err,rco,res,val1,tolsing);
        end if;
        Standard_Condition_Tables.Update_Corrector(t_err,err);
        Standard_Condition_Tables.Update_Condition(t_rco,rco);
        Standard_Condition_Tables.Update_Residuals(t_res,res);
        ls.rco := rco;
      else
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        info,startres,res,err,numit,fail,false);
        put(file,startres,3);
        put(file,"  #iterations : "); put(file,numit,1);
        if fail
         then put_line(file,"  failure");
         else put_line(file,"  success");
        end if;
        Standard_Complex_Solutions_io.put_vector(file,ls.v);
        put(file,"== err :"); put(file,err,3);
        put(file," = res :"); put(file,res,3);
        if fail then
          put_line(file," == no solution"); cntfail := cntfail + 1;
        else
          isreal := Standard_Solution_Diagnostics.Is_Real(ls.all,tolsing);
          val1 := Weight(ls.v,wv1);
          val2 := Weight(ls.v,wv2);
          Push(weights,val1,val2,cnt,ls);
          if isreal
           then put(file," == real solution");    cntreal := cntreal + 1;
           else put(file," == complex solution"); cntcmplx := cntcmplx + 1;
          end if;
        end if;
        Standard_Condition_Tables.Update_Corrector(t_err,err);
        Standard_Condition_Tables.Update_Residuals(t_res,res);
      end if;
      ls.err := err; ls.res := res;
      ptr := Tail_Of(ptr);
    end loop;
    tstop(timer);
    Standard_Complex_Solutions_io.put_bar(file);
    Count_Clusters(weights,tolsing,cntclus);
    if condition then
      put(file,"number of regular solutions   : ");
      put(file,cntregu,1); new_line(file);
      put(file,"number of singular solutions  : ");
      put(file,cntsing,1); new_line(file);
    end if;
    put(file,"number of real solutions      : ");
    put(file,cntreal,1); new_line(file);
    put(file,"number of complex solutions   : ");
    put(file,cntcmplx,1); new_line(file);
    put(file,"number of clustered solutions : ");
    put(file,cntclus,1); new_line(file);
    put(file,"number of failures            : ");
    put(file,cntfail,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
    new_line(file);
    print_times(file,timer,"Newton with condition table report");
  end Monitored_Run;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions,
  --   makes a coefficient circuit and then runs Newton's method.

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    nbq,len,dim : integer32 := 0;
    s : Link_to_System;
    ans : character;
    file : file_type;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(p,sols,"THE SOLUTIONS");
    nbq := p'last;
    len := integer32(Length_Of(sols));
    if len = 0 then
      put_line("No solutions found on file.");
    else
      dim := Head_Of(sols).n;
      new_line;
      put("Read system of "); put(nbq,1); put(" polynomials and ");
      put(len,1); put(" solutions in dimension "); put(dim,1); put_line(".");
      s := Make_Coefficient_System(p);
      new_line;
      put("Interactive run ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Interactive_Run(s,sols);
      else
        put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'n' then
          Monitored_Run(s,sols);
        else
          new_line;
          put_line("Reading the name of the output file ...");
          Read_Name_and_Create_File(file);
          Monitored_Run(file,s,sols); --,len);
          Close(file);
        end if;
      end if;
    end if;
  end Main;

begin
  Main;
end ts_newcirc;
