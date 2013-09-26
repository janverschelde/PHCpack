with Standard_Point_Coordinates;

package body Standard_Intrinsic_Solutions is

  function Project ( s : Solution; b : Vector; v : VecVec )
                   return Solution is

    res : Solution(v'length);

  begin
    res.t := s.t;
    res.m := s.m;
    res.v := Standard_Point_Coordinates.Project(s.v(b'range),b,v);
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Project;

  function Project ( s : Solution; p : Matrix ) return Solution is

    res : Solution(p'last(2));

  begin
    res.t := s.t;
    res.m := s.m;
    res.v := Standard_Point_Coordinates.Project(s.v(p'range(1)),p);
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Project;

  function Project ( sols : Solution_List; b : Vector; v : VecVec )
                   return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_list := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Project(ls.all,b,v));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Project;

  function Project ( sols : Solution_List; p : Matrix )
                   return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_list := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Project(ls.all,p));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Project;

  function Expand ( s : Solution; b : Vector; v : VecVec )
                  return Solution is

    res : Solution(b'length);

  begin
    res.t := s.t;
    res.m := s.m;
    res.v := Standard_Point_Coordinates.Affine_Expand(s.v,b,v);
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Expand;

  function Expand ( s : Solution; p : Matrix ) return Solution is

    res : Solution(p'last(1));

  begin
    res.t := s.t;
    res.m := s.m;
    res.v := Standard_Point_Coordinates.Affine_Expand(s.v,p);
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Expand;

  function Expand ( sols : Solution_List; b : Vector; v : VecVec )
                  return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Expand(ls.all,b,v));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Expand;

  function Expand ( sols : Solution_List; p : Matrix )
                  return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Expand(ls.all,p));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Expand;

  function Transform ( s : Solution; p1,p2 : Matrix ) return Solution is
  begin
    return Project(Expand(s,p1),p2);
  end Transform;

  function Transform ( sols : Solution_List; p1,p2 : Matrix )
                     return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_Last,Transform(Head_Of(tmp).all,p1,p2));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Transform;

  procedure Transform ( sols : in out Solution_List; p1,p2 : in Matrix ) is

    res : constant Solution_List := Transform(sols,p1,p2);

  begin
    Clear(sols);
    sols := res;
  end Transform;

end Standard_Intrinsic_Solutions;
