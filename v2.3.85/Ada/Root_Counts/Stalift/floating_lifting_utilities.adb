package body Floating_Lifting_Utilities is

  function Adaptive_Lifting ( L : Array_of_Lists ) return Vector is

    res : Vector(L'range);
    fac : constant double_float := 3.0;     -- multiplication factor
    max : constant double_float := 23.0;    -- upper bound for lifting

  begin
    for i in L'range loop
      res(i) := fac*double_float(Length_Of(l(i)));
      if res(i) > max
       then res(i) := max;
      end if;
    end loop;
    return res;
  end Adaptive_Lifting;

  procedure Search_Lifting ( L : in List; pt : in Vector;
                             found : out boolean; lif : out double_float ) is

    tmp : List := L;
    lpt : Link_to_Vector;

  begin
    found := false;
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if Equal(lpt(pt'range),pt) then
        found := true;
        lif := lpt(lpt'last);
        exit;
      else
        tmp := Tail_Of(tmp);
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

  function Occurred_Lifting ( mixsub : Mixed_Subdivision; k : integer32;
                              pt : Vector ) return Vector is

    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        lpt : constant Vector := Search_and_Lift(Head_Of(tmp),k,pt);
      begin
        if lpt'last > pt'last
         then return lpt;
         else tmp := Tail_Of(tmp);
        end if;
      end;
    end loop;
    return pt;
  end Occurred_Lifting;

  function Occurred_Lifting
             ( n : integer32; mix : Standard_Integer_Vectors.Vector;
               points : Array_of_Lists; mixsub : Mixed_Subdivision )
             return Array_of_Lists is

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
          lpt : constant Vector := Occurred_Lifting(mixsub,k,pt.all);
        begin
          if lpt'last > pt'last
           then Append(res(k),res_last(k),lpt);
          end if;
        end;
        tmp := Tail_Of(tmp);
      end loop;
      cnt := cnt + mix(k);
    end loop;
    return res;
  end Occurred_Lifting;

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
    res(res'last) := 1.0;
    res(res'last) := Conservative_Lifting(mixsub,k,res);
    return res;
  end Induced_Lifting;

  function Induced_Lifting
              ( n : integer32; mix : Standard_Integer_Vectors.Vector;
                points : Array_of_Lists; mixsub : Mixed_Subdivision )
              return Array_of_Lists is

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

  function Conservative_Lifting
             ( mic : Mixed_Cell; k : integer32; point : Vector )
             return double_float is

    sp : constant double_float := mic.nor*Head_Of(mic.pts(k));
    spp : double_float:= mic.nor.all*point;
    res : double_float;

  begin
    if sp < spp then
      return point(point'last);
    else
      if mic.nor(mic.nor'last) = 0.0 then
        res := point(point'last);
      else
        spp := spp - point(point'last)*mic.nor(mic.nor'last);
        res := (sp - spp)/mic.nor(mic.nor'last) + 1.0;
      end if;
      return res;
    end if;
  end Conservative_Lifting;

  function Conservative_Lifting
             ( mixsub : Mixed_Subdivision; k : integer32;
               point : Vector ) return double_float is

    tmp : Mixed_Subdivision := mixsub;
    pt : Vector(point'range) := point;
    res : double_float;

  begin
    while not Is_Null(tmp) loop
      pt(pt'last) := Conservative_Lifting(Head_Of(tmp),k,pt);
      tmp := Tail_Of(tmp);
    end loop;
    res := pt(pt'last);
    Clear(pt);
    return res;
  end Conservative_Lifting;

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

end Floating_Lifting_Utilities;
