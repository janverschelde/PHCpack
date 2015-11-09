package body Standard_Diagonal_Solutions is

  function Product ( s1,s2 : Solution ) return Solution is

    res : Solution(s1.n+s2.n);

  begin
    res.m := s1.m;
    res.t := s1.t;
    res.err := s1.err;
    res.rco := s1.rco;
    res.res := s1.res;
    res.v(s1.v'range) := s1.v;
    res.v(s1.v'last+1..res.v'last) := s2.v;
    return res;
  end Product;

  function Product ( s1,s2 : Solution_List ) return Solution_List is

    res,res_last,ptr1,ptr2 : Solution_List;

  begin
    ptr1 := s1;
    while not Is_Null(ptr1) loop
      ptr2 := s2;
      while not Is_Null(ptr2) loop
        Append(res,res_last,Product(Head_Of(ptr1).all,Head_Of(ptr2).all));
        ptr2 := Tail_Of(ptr2);
      end loop;
      ptr1 := Tail_Of(ptr1);
    end loop;
    return res;
  end Product;

  function Truncate ( s : Solution; n : integer32 ) return Solution is

    res : Solution(n);

  begin
    res.t := s.t;
    res.m := s.m;
    res.v := s.v(1..n);
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Truncate;

  function Truncate ( s : Solution_List; n : integer32 ) return Solution_List is

    res,res_last,tmp : Solution_List;

  begin
    tmp := s;
    while not Is_Null(tmp) loop
      Append(res,res_last,Truncate(Head_Of(tmp).all,n));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Truncate;

end Standard_Diagonal_Solutions;
