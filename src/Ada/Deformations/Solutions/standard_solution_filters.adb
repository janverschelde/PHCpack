with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Solution_Diagnostics;      use Standard_Solution_Diagnostics;

package body Standard_Solution_Filters is

-- CRITERIA :

  function Vanishing ( sol : Solution; tol : double_float ) return boolean is
  begin
    return ((abs(sol.err) < tol) or (abs(sol.res) < tol));
  end Vanishing;

  function Zero_Component ( sol : Solution; k : natural32;
                            tol : double_float ) return boolean is
  begin
    return (AbsVal(sol.v(integer32(k))) < tol);
  end Zero_Component;

  function Regular ( sol : Solution; tol : double_float ) return boolean is
  begin
    return (abs(sol.rco) > tol);
  end Regular;

  function On_Target ( sol : Solution; target : Complex_Number;
                       tol : double_float ) return boolean is
  begin
    return (AbsVal(sol.t-target) < tol);
  end On_Target;

-- PROTOTYPE FILTERS :

  function List_Filter ( sols : Solution_List ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      if Criterium(i,ls.all)
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end List_Filter;

  procedure Scan_Filter 
              ( infile,outfile : in file_type;
                len,dim : in natural32; cnt : out natural32 ) is

    s : Solution(integer32(dim));
    i : natural32 := 1;
    f : natural32 := 1024;  -- frequency updater

  begin
    cnt := 0;
    Write_First(outfile,len,dim);
    while i <= len loop
      Read_Next(infile,s);
      if Criterium(i,s)
       then Write_Next(outfile,cnt,s);
      end if;
      i := i + 1;
      if i mod f = 0
       then put(i,1); put(" ... "); f := 2*f;
      end if;
    end loop;
    if i >= 1024
     then new_line;
    end if;
  exception
    when others
       => put("Exception raised when reading solution "); put(i,1);
          put_line("."); return;
  end Scan_Filter;

  procedure Scan_Filter1
              ( infile,outfile : in file_type; len,dim : in natural32;
                ls : in Link_to_Solution; cnt : out natural32 ) is

    s : Solution(integer32(dim));
    i : natural32 := 2;
    f : natural32 := 1024;  -- frequency updater

  begin
    cnt := 0;
    Write_First(outfile,len,dim);
    if Criterium(1,ls.all)
     then Write_Next(outfile,cnt,ls.all);
    end if;
    while i <= len loop
      Read_Next(infile,s);
      if Criterium(i,s)
       then Write_Next(outfile,cnt,s);
      end if;
      i := i + 1;
      if i mod f = 0
       then put(i,1); put(" ... "); f := 2*f;
      end if;
    end loop;
    if i >= 1024
     then new_line;
    end if;
  exception
    when others
       => put("Exception raised when reading solution "); put(i,1);
          put_line("."); return;
  end Scan_Filter1;

-- FILTERS :

  function Vanishing_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List is

    function Vanishes ( i : natural32; s : Solution ) return boolean is
    begin
      return Vanishing(s,tol);
    end Vanishes;
    function Select_Vanishing is new List_Filter(Vanishes);

  begin
    return Select_Vanishing(sols);
  end Vanishing_Filter;

  procedure Vanishing_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 ) is

    function Vanishes ( i : natural32; s : Solution ) return boolean is
    begin
      return Vanishing(s,tol);
    end Vanishes;
    procedure Select_Vanishing is new Scan_Filter(Vanishes);

  begin
    Select_Vanishing(infile,outfile,len,dim,cnt);
  end Vanishing_Filter;

  function Spurious_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List is

    function Spurious ( i : natural32; s : Solution ) return boolean is
    begin
      return not Vanishing(s,tol);
    end Spurious;
    function Select_Spurious is new List_Filter(Spurious);

  begin
    return Select_Spurious(sols);
  end Spurious_Filter;

  procedure Spurious_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
	       tol : in double_float; cnt : out natural32 ) is

    function Spurious ( i : natural32; s : Solution ) return boolean is
    begin
      return not Vanishing(s,tol);
    end Spurious;
    procedure Select_Spurious is new Scan_Filter(Spurious);

  begin
    Select_Spurious(infile,outfile,len,dim,cnt);
  end Spurious_Filter;

  function Zero_Component_Filter
             ( sols : Solution_List; k : natural32; tol : double_float )
             return Solution_List is

    function Zero_Component ( i : natural32; s : Solution ) return boolean is
    begin
      return Zero_Component(s,k,tol);
    end Zero_Component;
    function Select_Zero_Component is new List_Filter(Zero_Component);

  begin
    return Select_Zero_Component(sols);
  end Zero_Component_Filter;

  procedure Zero_Component_Filter
             ( infile,outfile : in file_type; len,dim,k : in natural32;
               tol : in double_float; ls : in Link_to_Solution;
               cnt : out natural32 ) is

    function Zero_Component ( i : natural32; s : Solution ) return boolean is
    begin
      return Zero_Component(s,k,tol);
    end Zero_Component;
    procedure Select_Zero_Component is new Scan_Filter1(Zero_Component);

  begin
    Select_Zero_Component(infile,outfile,len,dim,ls,cnt);
  end Zero_Component_Filter;

  function Free_Component_Filter
             ( sols : Solution_List; k : natural32; tol : double_float )
             return Solution_List is

    function Free_Component ( i : natural32; s : Solution ) return boolean is
    begin
      return not Zero_Component(s,k,tol);
    end Free_Component;
    function Select_Zero_Component is new List_Filter(Free_Component);

  begin
    return Select_Zero_Component(sols);
  end Free_Component_Filter;

  procedure Free_Component_Filter
             ( infile,outfile : in file_type; len,dim,k : natural32;
               tol : in double_float; ls : in Link_to_Solution;
               cnt : out natural32 ) is

    function Free_Component ( i : natural32; s : Solution ) return boolean is
    begin
      return not Zero_Component(s,k,tol);
    end Free_Component;
    procedure Select_Zero_Component is new Scan_Filter1(Free_Component);

  begin
    Select_Zero_Component(infile,outfile,len,dim,ls,cnt);
  end Free_Component_Filter;

  function Regular_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List is

    function Regular ( i : natural32; s : Solution ) return boolean is
    begin
      return Regular(s,tol);
    end Regular;
    function Select_Regular is new List_Filter(Regular);

  begin
    return Select_Regular(sols);
  end Regular_Filter;

  procedure Regular_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 ) is

    function Regular ( i : natural32; s : Solution ) return boolean is
    begin
      return Regular(s,tol);
    end Regular;
    procedure Select_Regular is new Scan_Filter(Regular);

  begin
    Select_Regular(infile,outfile,len,dim,cnt);
  end Regular_Filter;

  function Singular_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List is

    function Singular ( i : natural32; s : Solution ) return boolean is
    begin
      return not Regular(s,tol);
    end Singular;
    function Select_Singular is new List_Filter(Singular);

  begin
    return Select_Singular(sols);
  end Singular_Filter;

  procedure Singular_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 ) is

    function Singular ( i : natural32; s : Solution ) return boolean is
    begin
      return not Regular(s,tol);
    end Singular;
    procedure Select_Singular is new Scan_Filter(Singular);

  begin
    Select_Singular(infile,outfile,len,dim,cnt);
  end Singular_Filter;

  function On_Target_Filter
             ( sols : Solution_List; target : Complex_Number;
               tol : double_float ) return Solution_List is

    function On_Target ( i : natural32; s : Solution ) return boolean is
    begin
      return On_Target(s,target,tol);
    end On_Target;
    function Select_On_Target is new List_Filter(On_Target);

  begin
    return Select_On_Target(sols);
  end On_Target_Filter;

  procedure On_Target_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               target : in Complex_Number; tol : in double_float;
               cnt : out natural32 ) is

    function On_Target ( i : natural32; s : Solution ) return boolean is
    begin
      return On_Target(s,target,tol);
    end On_Target;
    procedure Select_On_Target is new Scan_Filter(On_Target);

  begin
    Select_On_Target(infile,outfile,len,dim,cnt);
  end On_Target_Filter;

  function Off_Target_Filter
             ( sols : Solution_List; target : Complex_Number;
               tol : double_float ) return Solution_List is

    function Off_Target ( i : natural32; s : Solution ) return boolean is
    begin
      return not On_Target(s,target,tol);
    end Off_Target;
    function Select_Off_Target is new List_Filter(Off_Target);

  begin
    return Select_Off_Target(sols);
  end Off_Target_Filter;

  procedure Off_Target_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               target : in Complex_Number; tol : in double_float;
               cnt : out natural32 ) is

    function Off_Target ( i : natural32; s : Solution ) return boolean is
    begin
      return not On_Target(s,target,tol);
    end Off_Target;
    procedure Select_Off_Target is new Scan_Filter(Off_Target);

  begin
    Select_Off_Target(infile,outfile,len,dim,cnt);
  end Off_Target_Filter;

  function Real_Filter
             ( sols : Solution_List; tol : double_float )
             return Solution_List is

    function Real ( i : natural32; s : Solution ) return boolean is
    begin
      return Is_Real(s,tol);
    end Real;
    function Select_Real_Solutions is new List_Filter(Real);

  begin
    return Select_Real_Solutions(sols);
  end Real_Filter;

  procedure Real_Filter
             ( infile,outfile : in file_type; len,dim : in natural32;
               tol : in double_float; cnt : out natural32 ) is

    function Real ( i : natural32; s : Solution ) return boolean is
    begin
      return Is_Real(s,tol);
    end Real;
    procedure Select_Real_Solutions is new Scan_Filter(Real);

  begin
    Select_Real_Solutions(infile,outfile,len,dim,cnt);
  end Real_Filter;

  function Select_Solutions
             ( sols : Solution_List;
               nb : Standard_Natural_Vectors.Vector ) return Solution_List is

    res,res_last,tmp : Solution_List;
    ind : integer32 := nb'first;

  begin
    tmp := sols;
    for i in 1..Length_Of(sols) loop
      if i = nb(ind) then
        Append(res,res_last,Head_Of(tmp).all);
        ind := ind + 1;
      end if;
      exit when (ind > nb'last);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Select_Solutions;

  procedure Select_Solutions
             ( infile,outfile : in file_type; len,dim : in natural32;
               nb : in Standard_Natural_Vectors.Vector;
               cnt : out natural32 ) is

    s : Solution(integer32(dim));
    i : natural32 := 1;
    f : natural32 := 1024;
    ind : integer32 := nb'first;
    pos : natural32;

  begin
    cnt := 0;
    Write_First(outfile,nb'length,dim);
    while i <= len loop
      Read_Next(infile,s);
      if i = nb(ind) then
        pos := i-1;
        Write_Next(outfile,pos,s);
        ind := ind + 1;
        cnt := cnt + 1;
      end if;
      exit when (ind > nb'last);
      i := i + 1;
      if i mod f = 0
       then put(i,1); put(" ... "); f := 2*f;
      end if;
    end loop;
    if i >= 1024
     then new_line;
    end if;
  exception
    when others
       => put("Exception raised when reading solution "); put(i,1);
          put_line("."); return;
  end Select_Solutions;

  function Select_Failed_Solutions
              ( psols,qsols : Solution_List; tol : double_float;
                verbose : boolean := false ) return Solution_List is

    res,res_last : Solution_List;
    plist : Solution_List := psols;
    qlist : Solution_List := qsols;
    pls,qls : Link_to_Solution;
    cnt : natural32 := 0;
    target : constant Complex_Number := Create(1.0);

  begin
    while not Is_Null(plist) loop
      pls := Head_Of(plist);
      cnt := cnt + 1;
      if not On_Target(pls.all,target,tol) then
        qls := Head_Of(qlist);
        Append(res,res_last,qls.all);
        if verbose then
          put("Solution path "); put(cnt,1);
          put_line(" failed: did not reach 1.0.");
        end if;
      elsif not Vanishing(pls.all,tol) then
        qls := Head_Of(qlist);
        Append(res,res_last,qls.all);
        if verbose then
          put("Solution path "); put(cnt,1);
          put_line(" failed: is not vanishing.");
        end if;
      end if;
      plist := Tail_Of(plist);
      qlist := Tail_Of(qlist);
    end loop;
    return res;
  end Select_Failed_Solutions;

end Standard_Solution_Filters;
