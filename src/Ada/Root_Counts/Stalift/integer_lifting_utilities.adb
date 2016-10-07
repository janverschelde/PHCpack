with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Integer_Lifting_Utilities is

  function Adaptive_Lifting ( L : Array_of_Lists ) return Vector is

    res : Vector(L'range);
    fac : constant integer32 := 3;     -- multiplication factor
    max : constant integer32 := 23;    -- upper bound for lifting

  begin
    for i in L'range loop
      res(i) := fac*integer32(Length_Of(L(i)));
      if res(i) > max
       then res(i) := max;
      end if;
    end loop;
    return res;
  end Adaptive_Lifting;

  function Coefficient ( p : Poly; d : Standard_Integer_Vectors.Vector )
                       return Complex_Number is

    res : Complex_Number := Create(0.0);

    procedure Scan ( t : in Term; continue : out boolean ) is

      found : boolean := true;

    begin
      for i in d'range loop
        if t.dg(i) /= d(i)
         then found := false; exit;
        end if;
      end loop;
      if found then
        res := t.cf;
        continue := false;
      else
        continue := true;
      end if;
    end Scan;
    procedure Scan_Terms is new Visiting_Iterator(Scan);

  begin
    Scan_Terms(p);
    return res;
  end Coefficient;

  function Perform_Lifting
             ( n : integer32; L : List; p : Poly ) return Poly is

    res : Poly := Null_Poly;
    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      declare
        d : constant Link_to_Vector := Head_Of(tmp);
        t : Term;
      begin
        t.cf := Coefficient(p,d(1..n));
        t.dg := new Standard_Integer_Vectors.Vector'(d.all);
        Add(res,t);
        Clear(t);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Perform_Lifting;

  function Perform_Lifting
             ( n : integer32; mix : Vector; lifted : Array_of_Lists;
               p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);
    cnt : integer32 := 1;

  begin
    for k in mix'range loop
      for l in 1..mix(k) loop
        res(cnt) := Perform_Lifting(n,lifted(k),p(cnt));
        cnt := cnt+1;
      end loop;
    end loop;
    return res;
  end Perform_Lifting;

  function Copy_Lifting ( lifted : List; pt : Link_to_Vector )
                        return Link_to_Vector is

  -- DESCRIPTION :
  --   Searches the correspoinding point in the list lifted and returns
  --   the lifted point.  If the corresponding point has not been found,
  --   then the original point pt will be returned.

    tmp : List := lifted;
    lpt,res : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if Equal(lpt(pt'range),pt.all)
       then res := new Standard_Integer_Vectors.Vector'(lpt.all);
            return res;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return pt;
  end Copy_Lifting;

  function Copy_Lifting ( lifted,pts : List ) return List is

  -- DESCRIPTION :
  --   Copies the lifting on the points lifted to the points in pts,
  --   i.e., each point in pts will get the same lifting as the corresponding
  --   lifted point in the list lifted.

    res : List;
    tmp : List := pts;

  begin
    while not Is_Null(tmp) loop
      Construct(Copy_Lifting(lifted,Head_Of(tmp)),res);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Copy_Lifting;

  procedure Search_Lifting ( L : in List; pt : in Vector;
                             found : out boolean; lif : out integer32 ) is

    tmp : List := L;
    lpt : Link_to_Vector;

  begin
    found := false;
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if Equal(lpt(pt'range),pt)
       then found := true;
            lif := lpt(lpt'last); exit;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
  end Search_Lifting;

  function Search_and_Lift ( L : List; pt : Vector ) return Vector is

    tmp : List := L;
    lpt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if Equal(lpt(pt'range),pt)
       then return lpt.all;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return pt;
  end Search_and_Lift;

  function Search_and_Lift ( mic : Mixed_Cell; k : integer32; pt : Vector )
                           return Vector is
  begin
    return Search_and_Lift(mic.pts(k),pt);
  end Search_and_Lift;

  function Induced_Lifting ( mixsub : Mixed_Subdivision; k : integer32;
                             pt : Vector ) return Vector is

    tmp : Mixed_Subdivision := mixsub;
    res : Vector(pt'first..pt'last+1);

  begin
    while not Is_Null(tmp) loop
      declare
        mic : constant Mixed_Cell := Head_Of(tmp);
        lpt : constant Vector := Search_and_Lift(mic,k,pt);
      begin
        if lpt'length = pt'length+1
         then return lpt;
         else tmp := Tail_Of(tmp);
        end if;
      end;
    end loop;
    res(pt'range) := pt;
    res(res'last) := 1;
    res(res'last) := Conservative_Lifting(mixsub,k,res);
    return res;
  end Induced_Lifting;

  function Induced_Lifting
             ( n : integer32; mix : Vector; points : Array_of_Lists;
               mixsub : Mixed_Subdivision ) return Array_of_Lists is

    res,res_last : Array_of_Lists(mix'range);
    cnt : integer32 := 1;
    tmp : List;

  begin
    for k in mix'range loop
      res_last(k) := res(k);
      tmp := points(cnt);
      while not Is_Null(tmp) loop
        declare
          pt : constant Link_to_Vector := Head_Of(tmp);
          lpt : constant Vector := Induced_Lifting(mixsub,k,pt.all);
        begin
          Append(res(k),res_last(k),lpt);
        end;
        tmp := Tail_Of(tmp);
      end loop;
      cnt := cnt + mix(k);
    end loop;
    return res;
  end Induced_Lifting;

  procedure Constant_Lifting
                ( L : in List; liftval : in integer32;
                  lifted,lifted_last : in out List ) is

    tmp : List := L;
    pt : Link_to_Vector;
 
  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      declare
        lpt : constant Link_to_Vector := new Vector(pt'first..pt'last+1);
      begin
        lpt(pt'range) := pt.all;
        lpt(lpt'last) := liftval;
        Append(lifted,lifted_last,lpt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Constant_Lifting;

  procedure Constant_Lifting
               ( al : in Array_of_Lists; liftval : in integer32; 
                 lifted,lifted_last : in out Array_of_Lists ) is
  begin
    for i in al'range loop
      Constant_Lifting(al(i),liftval,lifted(i),lifted_last(i));
    end loop;
  end Constant_Lifting;

  function Conservative_Lifting
               ( mic : Mixed_Cell; k : integer32; point : Vector )
               return integer32 is

    sp : constant integer32 := mic.nor*Head_Of(mic.pts(k));
    spp : integer32 := mic.nor.all*point;
    res : integer32;

  begin
    if sp < spp then
      return point(point'last);
    else
      if mic.nor(mic.nor'last) = 0
       then res := point(point'last);
       else spp := spp - point(point'last)*mic.nor(mic.nor'last);
            res := (sp - spp)/mic.nor(mic.nor'last) + 1;
      end if;
      return res;
    end if;
  end Conservative_Lifting;

  function Conservative_Lifting ( mixsub : Mixed_Subdivision; k : integer32;
                                  point : Vector ) return integer32 is

    tmp : Mixed_Subdivision := mixsub;
    pt : Vector(point'range) := point;
    res : integer32;

  begin
    while not Is_Null(tmp) loop
      pt(pt'last) := Conservative_Lifting(Head_Of(tmp),k,pt);
      tmp := Tail_Of(tmp);
    end loop;
    res := pt(pt'last);
    Clear(pt);
    return res;
  end Conservative_Lifting;

  function Lower_Lifting ( mic : Mixed_Cell; k : integer32; point : Vector )
                         return integer32 is
  begin
    if Is_In(mic.pts(k),point) then
      return 0;
    else
      declare
        pt : Vector(point'range) := point;
      begin
        pt(pt'last) := 0;
        return Conservative_Lifting(mic,k,pt);
      end;
    end if;
  end Lower_Lifting;

  function Lower_Lifting ( mixsub : Mixed_Subdivision; k : integer32;
                           point : Vector ) return integer32 is

    lif : integer32 := point(point'last);
    tmp : Mixed_Subdivision := mixsub;
    max : integer32 := 0;

  begin
    while not Is_Null(tmp) loop
      lif := Lower_Lifting(Head_Of(tmp),k,point);
      if lif > max
       then max := lif;
      end if;
      exit when max = point(point'last);
      tmp := Tail_Of(tmp);
    end loop;
    return max;
  end Lower_Lifting;

  function Lifted_Supports
             ( r : integer32; mixsub : Mixed_Subdivision ) 
             return Array_of_Lists is

    res,res_last : Array_of_Lists(1..r);
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    ptr : List;
    ls : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      for i in res'range loop
        ptr := mic.pts(i);
        while not Is_Null(ptr) loop
          ls := Head_Of(ptr);
          Append_Diff(res(i),res_last(i),ls.all);
          ptr := Tail_Of(ptr);
        end loop;
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Lifted_Supports;

end Integer_Lifting_Utilities;
