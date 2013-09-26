package body Floating_Support_Functions is

  function Maximal_Support ( l : List; v : Vector ) return double_float is

    sp,max : double_float;
    tmp : List;

  begin
    if not Is_Null(l)
     then max := Head_Of(l).all*v;
          tmp := Tail_Of(l);
          while not Is_Null(tmp) loop
            sp := Head_Of(tmp).all*v;
            if sp > max
             then max := sp;
            end if;
            tmp := Tail_Of(tmp);
          end loop;
          return max;
     else return 0.0;
    end if;
  end Maximal_Support;

  function Minimal_Support ( l : List; v : Vector ) return double_float is

    sp,min : double_float;
    tmp : List;

  begin
    if not Is_Null(l)
     then min := Head_Of(l).all*v;
          tmp := Tail_Of(l);
          while not Is_Null(tmp) loop
            sp := Head_Of(tmp).all*v;
            if sp < min
             then min := sp;
            end if;
            tmp := Tail_Of(tmp);
          end loop;
          return min;
     else return 0.0;
    end if;
  end Minimal_Support;

  function Face ( l : List; v : Vector; m,tol : double_float ) return List is

    res,tmp,res_last : List;
    d : Vector(v'range);

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      d := Head_Of(tmp).all;
      if abs(d*v - m) < tol
       then Append(res,res_last,d);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Face;

end Floating_Support_Functions;
