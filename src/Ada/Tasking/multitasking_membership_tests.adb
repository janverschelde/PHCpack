with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Homotopy;
with Standard_Laurent_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Laurent_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Laurent_Homotopy;
with Homotopy_Membership_Target;
with Multitasking_Continuation;

package body Multitasking_Membership_Tests is

  function Is_Member ( s : Standard_Complex_Solutions.Solution_List;
                       z : Standard_Complex_Vectors.Vector;
                       tol : double_float ) return natural32 is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    ls : Link_to_Solution;
    tmp : Solution_List := s;
    found : boolean;
    diffsolz : Complex_Number;
    val : double_float;

  begin
    for i in 1..Length_Of(s) loop
      ls := Head_Of(tmp);
      found := true;
      for k in z'range loop
        diffsolz := ls.v(k) - z(k);
        val := AbsVal(diffsolz);
        if val > tol
         then found := false; exit;
        end if;
      end loop;
      if found
       then return i;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return 0;
  end Is_Member;

  function Is_Member ( s : DoblDobl_Complex_Solutions.Solution_List;
                       z : DoblDobl_Complex_Vectors.Vector;
                       tol : double_float ) return natural32 is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    ls : Link_to_Solution;
    tmp : Solution_List := s;
    found : boolean;
    diffsolz : Complex_Number;
    val : double_double;

  begin
    for i in 1..Length_Of(s) loop
      ls := Head_Of(tmp);
      found := true;
      for k in z'range loop
        diffsolz := ls.v(k) - z(k);
        val := AbsVal(diffsolz);
        if val > tol
         then found := false; exit;
        end if;
      end loop;
      if found
       then return i;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return 0;
  end Is_Member;

  function Is_Member ( s : QuadDobl_Complex_Solutions.Solution_List;
                       z : QuadDobl_Complex_Vectors.Vector;
                       tol : double_float ) return natural32 is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    ls : Link_to_Solution;
    tmp : Solution_List := s;
    found : boolean;
    diffsolz : Complex_Number;
    val : quad_double;

  begin
    for i in 1..Length_Of(s) loop
      ls := Head_Of(tmp);
      found := true;
      for k in z'range loop
        diffsolz := ls.v(k) - z(k);
        val := AbsVal(diffsolz);
        if val > tol
         then found := false; exit;
        end if;
      end loop;
      if found
       then return i;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return 0;
  end Is_Member;

  function Standard_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : Standard_Complex_Poly_Systems.Poly_Sys;
               gpts : Standard_Complex_Solutions.Solution_List;
               tpnt : Standard_Complex_Vectors.Vector )
             return natural32 is 

    res: natural32;
    target : Standard_Complex_Poly_Systems.Poly_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    newgpts : Standard_Complex_Solutions.Solution_List;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

    use Multitasking_Continuation;

  begin
    Standard_Homotopy.Create(target,start,1,gamma);
    Standard_Complex_Solutions.Copy(gpts,newgpts);
    Standard_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Path_Tracker(newgpts,integer32(nt));
    res := Is_Member(newgpts,tpnt,tol);
    Standard_Homotopy.Clear;
    Standard_Complex_Solutions.Clear(newgpts);
    return res;
  end Standard_Membership_Test;

  function Standard_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : Standard_Complex_Laur_Systems.Laur_Sys;
               gpts : Standard_Complex_Solutions.Solution_List;
               tpnt : Standard_Complex_Vectors.Vector )
             return natural32 is 

    res: natural32;
    target : Standard_Complex_Laur_Systems.Laur_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    newgpts : Standard_Complex_Solutions.Solution_List;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

    use Multitasking_Continuation;

  begin
    Standard_Laurent_Homotopy.Create(target,start,1,gamma);
    Standard_Complex_Solutions.Copy(gpts,newgpts);
    Standard_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Laurent_Path_Tracker(newgpts,integer32(nt));
    res := Is_Member(newgpts,tpnt,tol);
    Standard_Laurent_Homotopy.Clear;
    Standard_Complex_Solutions.Clear(newgpts);
    return res;
  end Standard_Membership_Test;

  function DoblDobl_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               gpts : DoblDobl_Complex_Solutions.Solution_List;
               tpnt : DoblDobl_Complex_Vectors.Vector )
             return natural32 is 

    res: natural32;
    target : DoblDobl_Complex_Poly_Systems.Poly_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;
    newgpts : DoblDobl_Complex_Solutions.Solution_List;
    ddzero : constant double_double := create(0.0);
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(ddzero);

    use Multitasking_Continuation;

  begin
    DoblDobl_Homotopy.Create(target,start,1,gamma);
    DoblDobl_Complex_Solutions.Copy(gpts,newgpts);
    DoblDobl_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Path_Tracker(newgpts,integer32(nt));
    res := Is_Member(newgpts,tpnt,tol);
    DoblDobl_Homotopy.Clear;
    DoblDobl_Complex_Solutions.Clear(newgpts);
    return res;
  end DoblDobl_Membership_Test;

  function DoblDobl_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               gpts : DoblDobl_Complex_Solutions.Solution_List;
               tpnt : DoblDobl_Complex_Vectors.Vector )
             return natural32 is 

    res: natural32;
    target : DoblDobl_Complex_Laur_Systems.Laur_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;
    newgpts : DoblDobl_Complex_Solutions.Solution_List;
    ddzero : constant double_double := create(0.0);
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(ddzero);

    use Multitasking_Continuation;

  begin
    DoblDobl_Laurent_Homotopy.Create(target,start,1,gamma);
    DoblDobl_Complex_Solutions.Copy(gpts,newgpts);
    DoblDobl_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Laurent_Path_Tracker(newgpts,integer32(nt));
    res := Is_Member(newgpts,tpnt,tol);
    DoblDobl_Laurent_Homotopy.Clear;
    DoblDobl_Complex_Solutions.Clear(newgpts);
    return res;
  end DoblDobl_Membership_Test;

  function QuadDobl_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               gpts : QuadDobl_Complex_Solutions.Solution_List;
               tpnt : QuadDobl_Complex_Vectors.Vector )
             return natural32 is 

    res: natural32;
    target : QuadDobl_Complex_Poly_Systems.Poly_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;
    newgpts : QuadDobl_Complex_Solutions.Solution_List;
    qdzero : constant quad_double := create(0.0);
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(qdzero);

    use Multitasking_Continuation;

  begin
    QuadDobl_Homotopy.Create(target,start,1,gamma);
    QuadDobl_Complex_Solutions.Copy(gpts,newgpts);
    QuadDobl_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Path_Tracker(newgpts,integer32(nt));
    res := Is_Member(newgpts,tpnt,tol);
    QuadDobl_Homotopy.Clear;
    QuadDobl_Complex_Solutions.Clear(newgpts);
    return res;
  end QuadDobl_Membership_Test;

  function QuadDobl_Membership_Test
             ( nt,dim : natural32; tol : double_float;
               start : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               gpts : QuadDobl_Complex_Solutions.Solution_List;
               tpnt : QuadDobl_Complex_Vectors.Vector )
             return natural32 is 

    res: natural32;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;
    newgpts : QuadDobl_Complex_Solutions.Solution_List;
    qdzero : constant quad_double := create(0.0);
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(qdzero);

    use Multitasking_Continuation;

  begin
    QuadDobl_Laurent_Homotopy.Create(target,start,1,gamma);
    QuadDobl_Complex_Solutions.Copy(gpts,newgpts);
    QuadDobl_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Laurent_Path_Tracker(newgpts,integer32(nt));
    res := Is_Member(newgpts,tpnt,tol);
    QuadDobl_Laurent_Homotopy.Clear;
    QuadDobl_Complex_Solutions.Clear(newgpts);
    return res;
  end QuadDobl_Membership_Test;

  procedure Standard_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in Standard_Complex_Poly_Systems.Poly_Sys;
               gpts : in Standard_Complex_Solutions.Solution_List;
               tpnt : in Standard_Complex_Vectors.Vector;
               idx : out natural32;
               match : out Standard_Complex_Vectors.Vector ) is

    target : Standard_Complex_Poly_Systems.Poly_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    newgpts : Standard_Complex_Solutions.Solution_List;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    tmp : Standard_Complex_Solutions.Solution_List;
    ls : Standard_Complex_Solutions.Link_to_Solution;

    use Multitasking_Continuation;

  begin
    Standard_Homotopy.Create(target,start,1,gamma);
    Standard_Complex_Solutions.Copy(gpts,newgpts);
    Standard_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Path_Tracker(newgpts,integer32(nt));
    idx := Is_Member(newgpts,tpnt,tol);
    if idx /= 0 then
      tmp := newgpts;
      for i in 1..Standard_Complex_Solutions.Length_Of(tmp) loop
        if i = idx then
          ls := Standard_Complex_Solutions.Head_Of(tmp);
          for k in match'range loop
            match(k) := ls.v(k);
          end loop;
        end if;
        exit when (i = idx);
        tmp := Standard_Complex_Solutions.Tail_Of(tmp);
      end loop;
    end if;
    Standard_Homotopy.Clear;
    Standard_Complex_Solutions.Clear(newgpts);
  end Standard_Membership_Test;

  procedure Standard_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in Standard_Complex_Laur_Systems.Laur_Sys;
               gpts : in Standard_Complex_Solutions.Solution_List;
               tpnt : in Standard_Complex_Vectors.Vector;
               idx : out natural32;
               match : out Standard_Complex_Vectors.Vector ) is

    target : Standard_Complex_Laur_Systems.Laur_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    newgpts : Standard_Complex_Solutions.Solution_List;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    tmp : Standard_Complex_Solutions.Solution_List;
    ls : Standard_Complex_Solutions.Link_to_Solution;

    use Multitasking_Continuation;

  begin
    Standard_Laurent_Homotopy.Create(target,start,1,gamma);
    Standard_Complex_Solutions.Copy(gpts,newgpts);
    Standard_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Laurent_Path_Tracker(newgpts,integer32(nt));
    idx := Is_Member(newgpts,tpnt,tol);
    if idx /= 0 then
      tmp := newgpts;
      for i in 1..Standard_Complex_Solutions.Length_Of(tmp) loop
        if i = idx then
          ls := Standard_Complex_Solutions.Head_Of(tmp);
          for k in match'range loop
            match(k) := ls.v(k);
          end loop;
        end if;
        exit when (i = idx);
        tmp := Standard_Complex_Solutions.Tail_Of(tmp);
      end loop;
    end if;
    Standard_Laurent_Homotopy.Clear;
    Standard_Complex_Solutions.Clear(newgpts);
  end Standard_Membership_Test;

  procedure DoblDobl_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               gpts : in DoblDobl_Complex_Solutions.Solution_List;
               tpnt : in DoblDobl_Complex_Vectors.Vector;
               idx : out natural32;
               match : out DoblDobl_Complex_Vectors.Vector ) is

    target : DoblDobl_Complex_Poly_Systems.Poly_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;
    newgpts : DoblDobl_Complex_Solutions.Solution_List;
    ddzero : constant double_double := create(0.0);
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(ddzero);
    tmp : DoblDobl_Complex_Solutions.Solution_List;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

    use Multitasking_Continuation;

  begin
    DoblDobl_Homotopy.Create(target,start,1,gamma);
    DoblDobl_Complex_Solutions.Copy(gpts,newgpts);
    DoblDobl_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Path_Tracker(newgpts,integer32(nt));
    idx := Is_Member(newgpts,tpnt,tol);
    if idx /= 0 then
      tmp := newgpts;
      for i in 1..DoblDobl_Complex_Solutions.Length_Of(tmp) loop
        if i = idx then
          ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
          for k in match'range loop
            match(k) := ls.v(k);
          end loop;
        end if;
        exit when (i = idx);
        tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
      end loop;
    end if;
    DoblDobl_Homotopy.Clear;
    DoblDobl_Complex_Solutions.Clear(newgpts);
  end DoblDobl_Membership_Test;

  procedure DoblDobl_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
               gpts : in DoblDobl_Complex_Solutions.Solution_List;
               tpnt : in DoblDobl_Complex_Vectors.Vector;
               idx : out natural32;
               match : out DoblDobl_Complex_Vectors.Vector ) is

    target : DoblDobl_Complex_Laur_Systems.Laur_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;
    newgpts : DoblDobl_Complex_Solutions.Solution_List;
    ddzero : constant double_double := create(0.0);
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(ddzero);
    tmp : DoblDobl_Complex_Solutions.Solution_List;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

    use Multitasking_Continuation;

  begin
    DoblDobl_Laurent_Homotopy.Create(target,start,1,gamma);
    DoblDobl_Complex_Solutions.Copy(gpts,newgpts);
    DoblDobl_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Laurent_Path_Tracker(newgpts,integer32(nt));
    idx := Is_Member(newgpts,tpnt,tol);
    if idx /= 0 then
      tmp := newgpts;
      for i in 1..DoblDobl_Complex_Solutions.Length_Of(tmp) loop
        if i = idx then
          ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
          for k in match'range loop
            match(k) := ls.v(k);
          end loop;
        end if;
        exit when (i = idx);
        tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
      end loop;
    end if;
    DoblDobl_Laurent_Homotopy.Clear;
    DoblDobl_Complex_Solutions.Clear(newgpts);
  end DoblDobl_Membership_Test;

  procedure QuadDobl_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               gpts : in QuadDobl_Complex_Solutions.Solution_List;
               tpnt : in QuadDobl_Complex_Vectors.Vector;
               idx : out natural32;
               match : out QuadDobl_Complex_Vectors.Vector ) is

    target : QuadDobl_Complex_Poly_Systems.Poly_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;
    newgpts : QuadDobl_Complex_Solutions.Solution_List;
    qdzero : constant quad_double := create(0.0);
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(qdzero);
    tmp : QuadDobl_Complex_Solutions.Solution_List;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

    use Multitasking_Continuation;

  begin
    QuadDobl_Homotopy.Create(target,start,1,gamma);
    QuadDobl_Complex_Solutions.Copy(gpts,newgpts);
    QuadDobl_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Path_Tracker(newgpts,integer32(nt));
    idx := Is_Member(newgpts,tpnt,tol);
    if idx /= 0 then
      tmp := newgpts;
      for i in 1..QuadDobl_Complex_Solutions.Length_Of(tmp) loop
        if i = idx then
          ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
          for k in match'range loop
            match(k) := ls.v(k);
          end loop;
        end if;
        exit when (i = idx);
        tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
      end loop;
    end if;
    QuadDobl_Homotopy.Clear;
    QuadDobl_Complex_Solutions.Clear(newgpts);
  end QuadDobl_Membership_Test;

  procedure QuadDobl_Membership_Test
             ( nt,dim : in natural32; tol : in double_float;
               start : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
               gpts : in QuadDobl_Complex_Solutions.Solution_List;
               tpnt : in QuadDobl_Complex_Vectors.Vector;
               idx : out natural32;
               match : out QuadDobl_Complex_Vectors.Vector ) is

    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(start'range)
           := Homotopy_Membership_Target.Adjusted_Target(start,dim,tpnt);
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;
    newgpts : QuadDobl_Complex_Solutions.Solution_List;
    qdzero : constant quad_double := create(0.0);
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(qdzero);
    tmp : QuadDobl_Complex_Solutions.Solution_List;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

    use Multitasking_Continuation;

  begin
    QuadDobl_Laurent_Homotopy.Create(target,start,1,gamma);
    QuadDobl_Complex_Solutions.Copy(gpts,newgpts);
    QuadDobl_Complex_Solutions.Set_Continuation_Parameter(newgpts,zero);
    Silent_Multitasking_Laurent_Path_Tracker(newgpts,integer32(nt));
    idx := Is_Member(newgpts,tpnt,tol);
    if idx /= 0 then
      tmp := newgpts;
      for i in 1..QuadDobl_Complex_Solutions.Length_Of(tmp) loop
        if i = idx then
          ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
          for k in match'range loop
            match(k) := ls.v(k);
          end loop;
        end if;
        exit when (i = idx);
        tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
      end loop;
    end if;
    QuadDobl_Laurent_Homotopy.Clear;
    QuadDobl_Complex_Solutions.Clear(newgpts);
  end QuadDobl_Membership_Test;

  function Standard_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : Standard_Complex_Poly_Systems.Poly_Sys;
               gpts,sols : Standard_Complex_Solutions.Solution_List )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    idx : natural32;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      idx := Standard_Membership_Test(nt,dim,tol,start,gpts,ls.v);
      if idx = 0
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Standard_Membership_Filter;

  function Standard_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : Standard_Complex_Laur_Systems.Laur_Sys;
               gpts,sols : Standard_Complex_Solutions.Solution_List )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    idx : natural32;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      idx := Standard_Membership_Test(nt,dim,tol,start,gpts,ls.v);
      if idx = 0
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Standard_Membership_Filter;

  function DoblDobl_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               gpts,sols : DoblDobl_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    idx : natural32;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      idx := DoblDobl_Membership_Test(nt,dim,tol,start,gpts,ls.v);
      if idx = 0
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end DoblDobl_Membership_Filter;

  function DoblDobl_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               gpts,sols : DoblDobl_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    idx : natural32;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      idx := DoblDobl_Membership_Test(nt,dim,tol,start,gpts,ls.v);
      if idx = 0
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end DoblDobl_Membership_Filter;

  function QuadDobl_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               gpts,sols : QuadDobl_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    idx : natural32;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      idx := QuadDobl_Membership_Test(nt,dim,tol,start,gpts,ls.v);
      if idx = 0
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end QuadDobl_Membership_Filter;

  function QuadDobl_Membership_Filter
             ( nt,dim : natural32; tol : double_float;
               start : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               gpts,sols : QuadDobl_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    idx : natural32;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      idx := QuadDobl_Membership_Test(nt,dim,tol,start,gpts,ls.v);
      if idx = 0
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end QuadDobl_Membership_Filter;

end Multitasking_Membership_Tests;
