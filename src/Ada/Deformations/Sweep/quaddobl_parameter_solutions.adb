package body QuadDobl_Parameter_Solutions is

  function Select_Variables
              ( s : Solution; nv : integer32;
                v : Standard_Integer_Vectors.Vector ) return Solution is

    res : Solution(nv);

  begin
    res.t := s.t;
    res.m := s.m;
    for i in v'range loop
      res.v(i) := s.v(v(i));
    end loop;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Select_Variables;

  function Select_Variables
              ( s : Solution_List; nv : integer32;
                v : Standard_Integer_Vectors.Vector ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Select_Variables(ls.all,nv,v));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Select_Variables;

  function Join_Variables
              ( s : Solution; n : integer32;
                vars,pars : Standard_Integer_Vectors.Vector;
                val_pars : QuadDobl_Complex_Vectors.Vector )
              return Solution is

    res : Solution(n);

  begin
    res.t := s.t;
    res.m := s.m;
    for i in vars'range loop
      res.v(vars(i)) := s.v(i);
    end loop;
    for i in pars'range loop
      res.v(pars(i)) := val_pars(i);
    end loop;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Join_Variables;

  function Join_Variables
              ( s : Solution_List; n : integer32;
                vars,pars : Standard_Integer_Vectors.Vector;
                val_pars : QuadDobl_Complex_Vectors.Vector )
              return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := s;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
	ns : constant Solution(n)
           := Join_Variables(ls.all,n,vars,pars,val_pars);
      begin
        Append(res,res_last,ns);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Join_Variables;

end QuadDobl_Parameter_Solutions; 
