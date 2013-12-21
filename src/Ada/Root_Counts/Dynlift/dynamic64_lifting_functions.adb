package body Dynamic64_Lifting_Functions is

  function Lift_to_Place ( s : Simplex; x : Vector ) return integer64 is 

    nor : constant Vector := Normal(s);
    ips : integer64 := nor*Vertex(s,1);
    ipx : integer64 := nor*x;

  begin
    if ips < ipx
     then return x(x'last);
     else ipx := ipx - x(x'last)*nor(nor'last);
          return ((ips - ipx)/nor(nor'last) + 1);
    end if;
  end Lift_to_Place;

  function Lift_to_Place ( t : Triangulation; x : Vector ) return integer64 is

    tmp : Triangulation := t;
    wrk : Vector(x'range) := x;

  begin
    while not Is_Null(tmp) loop
      wrk(wrk'last) := Lift_to_Place(Head_Of(tmp),wrk);
      tmp := Tail_Of(tmp);
    end loop;
    return wrk(wrk'last);
  end Lift_to_Place;

  function Lift_to_Pull ( s : Simplex; x : Vector ) return integer64 is

    nor : constant Vector := Normal(s);
    ips : integer64 := nor*Vertex(s,1);
    ipx : integer64 := nor*x;

  begin
    if ipx < ips
     then return x(x'last);
     else ipx := ipx - x(x'last)*nor(nor'last);
          return ((ips - ipx)/nor(nor'last) - 1);
    end if;
  end Lift_to_Pull;

  function Lift_to_Pull ( t : Triangulation; x : Vector ) return integer64 is

    tmp : Triangulation := t;
    wrk : Vector(x'range) := x;

  begin
    while not Is_Null(tmp) loop
      wrk(wrk'last) := Lift_to_Pull(Head_Of(tmp),wrk);
      tmp := Tail_Of(tmp);
    end loop;
    return wrk(wrk'last);
  end Lift_to_Pull;

  function Degenerate ( t : Triangulation; x : Vector ) return boolean is

    tmp : Triangulation := t;
    s : Simplex;

  begin
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      declare
        nor : constant Vector := Normal(s);
        apt : constant Vector := Vertex(s,1);
        ipx : constant integer64 := x*nor;
      begin
        if apt*nor = ipx
         then return true;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return false;
  end Degenerate;

  function Lift_to_Pull
             ( t1,t2 : Triangulation; x : Vector ) return integer64 is

    wrk : Vector(x'range) := x;

  begin
    wrk(wrk'last) := Lift_to_Pull(t1,x);
    while Degenerate(t2,wrk) loop   -- pull the lifting further down
      wrk(wrk'last) := wrk(wrk'last) - 1;
    end loop;
    return wrk(wrk'last);
  end Lift_to_Pull;

end Dynamic64_Lifting_Functions;
