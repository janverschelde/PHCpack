with DoblDobl_Complex_Numbers_cv;
with DoblDobl_Complex_Vectors_cv;
with QuadDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Vectors_cv;
with Multprec_Floating_Numbers;
with Multprec_Complex_Number_Tools;
with Multprec_DoblDobl_Convertors;
with Multprec_QuadDobl_Convertors;
with Multprec_Complex_Vector_Tools;

package body Varbprec_Complex_Solutions is

  function Multprec_to_Standard_Solution
             ( s : Multprec_Complex_Solutions.Solution )
             return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(s.n);

  begin
    res.t := Multprec_Complex_Number_Tools.Round(s.t);
    res.m := s.m;
    res.v := Multprec_Complex_Vector_Tools.Round(s.v);
    res.err := Multprec_Floating_Numbers.Round(s.err);
    res.rco := Multprec_Floating_Numbers.Round(s.rco);
    res.res := Multprec_Floating_Numbers.Round(s.res);
    return res;
  end Multprec_to_Standard_Solution;

  function Multprec_to_Standard_Solutions
             ( s : Multprec_Complex_Solutions.Solution_List )
             return Standard_Complex_Solutions.Solution_List is

    res,res_last : Standard_Complex_Solutions.Solution_List;
    tmp : Multprec_Complex_Solutions.Solution_List := s;
    ls : Multprec_Complex_Solutions.Link_to_Solution;

  begin
    while not Multprec_Complex_Solutions.Is_Null(tmp) loop
      ls := Multprec_Complex_Solutions.Head_Of(tmp);
      declare
        sol : constant Standard_Complex_Solutions.Solution(ls.n)
            := Multprec_to_Standard_Solution(ls.all);
      begin
        Standard_Complex_Solutions.Append(res,res_last,sol);
      end;
      tmp := Multprec_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Multprec_to_Standard_Solutions;

  function Multprec_to_DoblDobl_Solution
             ( s : Multprec_Complex_Solutions.Solution )
             return DoblDobl_Complex_Solutions.Solution is

    res : DoblDobl_Complex_Solutions.Solution(s.n);

  begin
    res.t := DoblDobl_Complex_Numbers_cv.Multprec_to_DoblDobl_Complex(s.t);
    res.m := s.m;
    res.v := DoblDobl_Complex_Vectors_cv.Multprec_to_DoblDobl_Complex(s.v);
    res.err := Multprec_DoblDobl_Convertors.to_double_double(s.err);
    res.rco := Multprec_DoblDobl_Convertors.to_double_double(s.rco);
    res.res := Multprec_DoblDobl_Convertors.to_double_double(s.res);
    return res;
  end Multprec_to_DoblDobl_Solution;

  function Multprec_to_DoblDobl_Solutions
             ( s : Multprec_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Solutions.Solution_List is

    res,res_last : DoblDobl_Complex_Solutions.Solution_List;
    tmp : Multprec_Complex_Solutions.Solution_List := s;
    ls : Multprec_Complex_Solutions.Link_to_Solution;

  begin
    while not Multprec_Complex_Solutions.Is_Null(tmp) loop
      ls := Multprec_Complex_Solutions.Head_Of(tmp);
      declare
        sol : constant DoblDobl_Complex_Solutions.Solution(ls.n)
            := Multprec_to_DoblDobl_Solution(ls.all);
      begin
        DoblDobl_Complex_Solutions.Append(res,res_last,sol);
      end;
      tmp := Multprec_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Multprec_to_DoblDobl_Solutions;

  function Multprec_to_QuadDobl_Solution
             ( s : Multprec_Complex_Solutions.Solution )
             return QuadDobl_Complex_Solutions.Solution is

    res : QuadDobl_Complex_Solutions.Solution(s.n);

  begin
    res.t := QuadDobl_Complex_Numbers_cv.Multprec_to_QuadDobl_Complex(s.t);
    res.m := s.m;
    res.v := QuadDobl_Complex_Vectors_cv.Multprec_to_QuadDobl_Complex(s.v);
    res.err := Multprec_QuadDobl_Convertors.to_quad_double(s.err);
    res.rco := Multprec_QuadDobl_Convertors.to_quad_double(s.rco);
    res.res := Multprec_QuadDobl_Convertors.to_quad_double(s.res);
    return res;
  end Multprec_to_QuadDobl_Solution;

  function Multprec_to_QuadDobl_Solutions
             ( s : Multprec_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Solutions.Solution_List is

    res,res_last : QuadDobl_Complex_Solutions.Solution_List;
    tmp : Multprec_Complex_Solutions.Solution_List := s;
    ls : Multprec_Complex_Solutions.Link_to_Solution;

  begin
    while not Multprec_Complex_Solutions.Is_Null(tmp) loop
      ls := Multprec_Complex_Solutions.Head_Of(tmp);
      declare
        sol : constant QuadDobl_Complex_Solutions.Solution(ls.n)
            := Multprec_to_QuadDobl_Solution(ls.all);
      begin
        QuadDobl_Complex_Solutions.Append(res,res_last,sol);
      end;
      tmp := Multprec_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Multprec_to_QuadDobl_Solutions;

  function Standard_to_Multprec_Solution
             ( s : Standard_Complex_Solutions.Solution; size : natural32 )
             return Multprec_Complex_Solutions.Solution is

    res : Multprec_Complex_Solutions.Solution(s.n);

  begin
    res.t := Multprec_Complex_Number_Tools.Create(s.t);
    res.m := s.m;
    res.v := Multprec_Complex_Vector_Tools.Create(s.v);
    Multprec_Complex_Vector_Tools.Set_Size(res.v,size);
    res.err := Multprec_Floating_Numbers.create(s.err);
    res.rco := Multprec_Floating_Numbers.create(s.rco);
    res.res := Multprec_Floating_Numbers.create(s.res);
    return res;
  end Standard_to_Multprec_Solution;

  function Standard_to_Multprec_Solutions
             ( s : Standard_Complex_Solutions.Solution_List; size : natural32 )
             return Multprec_Complex_Solutions.Solution_List is

    res,res_last : Multprec_Complex_Solutions.Solution_List;
    tmp : Standard_Complex_Solutions.Solution_List := s;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      declare
        sol : constant Multprec_Complex_Solutions.Solution(ls.n)
            := Standard_to_Multprec_Solution(ls.all,size);
      begin
        Multprec_Complex_Solutions.Append(res,res_last,sol);
      end;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Standard_to_Multprec_Solutions;

  function DoblDobl_to_Multprec_Solution
             ( s : DoblDobl_Complex_Solutions.Solution; size : natural32 )
             return Multprec_Complex_Solutions.Solution is

    res : Multprec_Complex_Solutions.Solution(s.n);

  begin
    res.t := DoblDobl_Complex_Numbers_cv.DoblDobl_Complex_to_Multprec(s.t);
    res.m := s.m;
    res.v := DoblDobl_Complex_Vectors_cv.DoblDobl_Complex_to_Multprec(s.v);
    res.err := Multprec_DoblDobl_Convertors.to_floating_number(s.err);
    res.rco := Multprec_DoblDobl_Convertors.to_floating_number(s.rco);
    res.res := Multprec_DoblDobl_Convertors.to_floating_number(s.res);
    return res;
  end DoblDobl_to_Multprec_Solution;

  function DoblDobl_to_Multprec_Solutions
             ( s : DoblDobl_Complex_Solutions.Solution_List; size : natural32 )
             return Multprec_Complex_Solutions.Solution_List is

    res,res_last : Multprec_Complex_Solutions.Solution_List;
    tmp : DoblDobl_Complex_Solutions.Solution_List := s;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sol : constant Multprec_Complex_Solutions.Solution(ls.n)
            := DoblDobl_to_Multprec_Solution(ls.all,size);
      begin
        Multprec_Complex_Solutions.Append(res,res_last,sol);
      end;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end DoblDobl_to_Multprec_Solutions;

  function QuadDobl_to_Multprec_Solution
             ( s : QuadDobl_Complex_Solutions.Solution; size : natural32 )
             return Multprec_Complex_Solutions.Solution is

    res : Multprec_Complex_Solutions.Solution(s.n);

  begin
    res.t := QuadDobl_Complex_Numbers_cv.QuadDobl_Complex_to_Multprec(s.t);
    res.m := s.m;
    res.v := QuadDobl_Complex_Vectors_cv.QuadDobl_Complex_to_Multprec(s.v);
    res.err := Multprec_QuadDobl_Convertors.to_floating_number(s.err);
    res.rco := Multprec_QuadDobl_Convertors.to_floating_number(s.rco);
    res.res := Multprec_QuadDobl_Convertors.to_floating_number(s.res);
    return res;
  end QuadDobl_to_Multprec_Solution;

  function QuadDobl_to_Multprec_Solutions
             ( s : QuadDobl_Complex_Solutions.Solution_List; size : natural32 )
             return Multprec_Complex_Solutions.Solution_List is

    res,res_last : Multprec_Complex_Solutions.Solution_List;
    tmp : QuadDobl_Complex_Solutions.Solution_List := s;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sol : constant Multprec_Complex_Solutions.Solution(ls.n)
            := QuadDobl_to_Multprec_Solution(ls.all,size);
      begin
        Multprec_Complex_Solutions.Append(res,res_last,sol);
      end;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end QuadDobl_to_Multprec_Solutions;


end Varbprec_Complex_Solutions;
