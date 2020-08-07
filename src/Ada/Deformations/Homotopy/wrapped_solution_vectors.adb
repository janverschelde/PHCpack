with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;

package body Wrapped_Solution_Vectors is

  function Create ( xt : Standard_Complex_Vectors.Vector )
                  return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(xt'last-1);

  begin
    res.t := xt(xt'last);
    res.m := 1;
    res.v := xt(xt'first..xt'last-1);
    res.err := 0.0;
    res.rco := 1.0;
    res.res := 0.0;
    return res;
  end Create;

  function Create ( xt : DoblDobl_Complex_Vectors.Vector )
                  return DoblDobl_Complex_Solutions.Solution is

    res : DoblDobl_Complex_Solutions.Solution(xt'last-1);

  begin
    res.t := xt(xt'last);
    res.m := 1;
    res.v := xt(xt'first..xt'last-1);
    res.err := create(0.0);
    res.rco := create(1.0);
    res.res := create(0.0);
    return res;
  end Create;

  function Create ( xt : QuadDobl_Complex_Vectors.Vector )
                  return QuadDobl_Complex_Solutions.Solution is

    res : QuadDobl_Complex_Solutions.Solution(xt'last-1);

  begin
    res.t := xt(xt'last);
    res.m := 1;
    res.v := xt(xt'first..xt'last-1);
    res.err := create(0.0);
    res.rco := create(1.0);
    res.res := create(0.0);
    return res;
  end Create;

  function Create ( xt : Standard_Complex_Vectors.Vector ) 
                  return Standard_Complex_Solutions.Solution_List is

    sol : constant Standard_Complex_Solutions.Solution(xt'last-1)
        := Create(xt);
    res : Standard_Complex_Solutions.Solution_List;

  begin
    Standard_Complex_Solutions.Add(res,sol);
    return res;
  end Create;

  function Create ( xt : DoblDobl_Complex_Vectors.Vector ) 
                  return DoblDobl_Complex_Solutions.Solution_List is

    sol : constant DoblDobl_Complex_Solutions.Solution(xt'last-1)
        := Create(xt);
    res : DoblDobl_Complex_Solutions.Solution_List;

  begin
    DoblDobl_Complex_Solutions.Add(res,sol);
    return res;
  end Create;

  function Create ( xt : QuadDobl_Complex_Vectors.Vector ) 
                  return QuadDobl_Complex_Solutions.Solution_List is

    sol : constant QuadDobl_Complex_Solutions.Solution(xt'last-1)
        := Create(xt);
    res : QuadDobl_Complex_Solutions.Solution_List;

  begin
    QuadDobl_Complex_Solutions.Add(res,sol);
    return res;
  end Create;

  function Create ( xt : Standard_Complex_Solutions.Link_to_Solution )
                  return Standard_Complex_Solutions.Link_to_Solution is
 
    sol : Standard_Complex_Solutions.Solution(xt.n-1);
    res : Standard_Complex_Solutions.Link_to_Solution;

  begin
    sol.t := xt.v(xt.v'last);
    sol.m := 1;
    sol.v := xt.v(xt.v'first..xt.v'last-1);
    sol.err := xt.err;
    sol.rco := xt.rco;
    sol.res := xt.res;
    res := new Standard_Complex_Solutions.Solution'(sol);
    return res;
  end Create;

  function Create ( xt : DoblDobl_Complex_Solutions.Link_to_Solution )
                  return DoblDobl_Complex_Solutions.Link_to_Solution is
 
    sol : DoblDobl_Complex_Solutions.Solution(xt.n-1);
    res : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    sol.t := xt.v(xt.v'last);
    sol.m := 1;
    sol.v := xt.v(xt.v'first..xt.v'last-1);
    sol.err := xt.err;
    sol.rco := xt.rco;
    sol.res := xt.res;
    res := new DoblDobl_Complex_Solutions.Solution'(sol);
    return res;
  end Create;

  function Create ( xt : QuadDobl_Complex_Solutions.Link_to_Solution )
                  return QuadDobl_Complex_Solutions.Link_to_Solution is
 
    sol : QuadDobl_Complex_Solutions.Solution(xt.n-1);
    res : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    sol.t := xt.v(xt.v'last);
    sol.m := 1;
    sol.v := xt.v(xt.v'first..xt.v'last-1);
    sol.err := xt.err;
    sol.rco := xt.rco;
    sol.res := xt.res;
    res := new QuadDobl_Complex_Solutions.Solution'(sol);
    return res;
  end Create;

  function Create ( xt : Standard_Complex_Solutions.Solution_List )
                  return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := xt;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Create(ls));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( xt : DoblDobl_Complex_Solutions.Solution_List )
                  return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := xt;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Create(ls));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( xt : QuadDobl_Complex_Solutions.Solution_List )
                  return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := xt;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Create(ls));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

-- UPDATING SOLUTION LISTS :

  procedure Update ( sols : in out Standard_Complex_Solutions.Solution_List;
                     xtsols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    xtp : Solution_List := xtsols;
    ls,xtls : Link_to_Solution;

  begin
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
  end Update;

  procedure Update ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                     xtsols : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    xtp : Solution_List := xtsols;
    ls,xtls : Link_to_Solution;

  begin
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
  end Update;

  procedure Update ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                     xtsols : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    xtp : Solution_List := xtsols;
    ls,xtls : Link_to_Solution;

  begin
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
  end Update;

end Wrapped_Solution_Vectors;
