with Graded_Lexicographic_Order;         use Graded_Lexicographic_Order;

package body Integer_Support_Functions is

  function Maximal_Support ( L : List; v : Vector ) return integer32 is

    sp,max : integer32;
    tmp : List;

  begin
    if not Is_Null(L) then
      max := Head_Of(L).all*v;
      tmp := Tail_Of(L);
      while not Is_Null(tmp) loop
        sp := Head_Of(tmp).all*v;
        if sp > max
         then max := sp;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      return max;
    else
      return 0;
    end if;
  end Maximal_Support;

  function Minimal_Support ( L : List; v : Vector ) return integer32 is

    sp,min : integer32;
    tmp : List;

  begin
    if not Is_Null(L) then
      min := Head_Of(L).all*v;
      tmp := Tail_Of(L);
      while not Is_Null(tmp) loop
        sp := Head_Of(tmp).all*v;
        if sp < min
         then min := sp;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      return min;
    else
     return 0;
    end if;
  end Minimal_Support;

  procedure Min_Max ( L : in List; k : in integer32;
                      min,max : in out integer32 ) is

    tmp : List;
    v : Link_to_Vector;

  begin
    if not Is_Null(L) then
      tmp := L;
      v := Head_Of(tmp);
      min := v(k);  max := min;
      tmp := Tail_Of(tmp);
      while not Is_Null(tmp) loop
        v := Head_Of(tmp);
        if v(k) < min then
          min := v(k);
        elsif v(k) > max then
          max := v(k);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Min_Max;

  function Graded_Max ( L : List ) return Link_to_Vector is

    res : constant Link_to_Vector := new Vector'(Head_Of(l).all);
    tmp : List := Tail_Of(L);
    ele : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      ele := Head_Of(tmp);
      if ele > res
       then res.all := ele.all;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Graded_Max;

  function Face ( L : List; v : Vector; m : integer32 ) return List is

    res,tmp,res_last : List;
    d : Vector(v'range);

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      d := Head_Of(tmp).all;
      if d*v = m
       then Append(res,res_last,d);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Face;

  function Inner_Face ( L : List; v : Vector ) return List is
  begin
    return Face(L,v,Minimal_Support(L,v));
  end Inner_Face;

  function Outer_Face ( L : List; v : Vector ) return List is
  begin
    return Face(L,v,Maximal_Support(L,v));
  end Outer_Face;

end Integer_Support_Functions;
