with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Symbol_Table;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Floating_Vectors;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_to_Real_Poly;
with Standard_Floating_Poly_Systems;
with Standard_Floating_Poly_SysFun;
with Standard_Floating_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_to_Real_Poly;
with Double_Double_Poly_Systems;
with Double_Double_Poly_SysFun;
with Double_Double_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_to_Real_Poly;
with Quad_Double_Poly_Systems;
with Quad_Double_Poly_SysFun;
with Quad_Double_Jaco_Matrices;
with Standard_Complex_Solutions;
with Standard_Solution_Diagnostics;
with DoblDobl_Complex_Solutions;
with DoblDobl_Solution_Diagnostics;
with QuadDobl_Complex_Solutions;
with QuadDobl_Solution_Diagnostics;
with Standard_Parameter_Systems;
with Standard_Parameter_Solutions;
with Standard_Quad_Sweepers;
with DoblDobl_Parameter_Systems;
with DoblDobl_Parameter_Solutions;
with DoblDobl_Quad_Sweepers;
with QuadDobl_Parameter_Systems;
with QuadDobl_Parameter_Solutions;
with QuadDobl_Quad_Sweepers;
with Complex_Convex_Continuation;
with Parameter_Homotopy_State;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
-- output for debugging purposes:
--with Symbol_Table_io;
with Standard_Natural_Numbers_io; use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io; use Standard_Integer_Numbers_io;
--with Standard_Integer_Vectors_io; use Standard_Integer_Vectors_io;
--with Standard_Complex_Solutions_io; use Standard_Complex_Solutions_io;

function use_sweep ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is -- define parameters numerically

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nbequ : constant integer32 := integer32(v_a(v_a'first));
    nbvar : constant integer32 := integer32(v_a(v_a'first+1));
    nbpar : constant integer32 := integer32(v_a(v_a'first+2));
    idx : Standard_Integer_Vectors.Vector(1..nbpar);

  begin
    Parameter_Homotopy_State.Set_Number_of_Equations(nbequ);
    Parameter_Homotopy_State.Set_Number_of_Variables(nbvar);
    Parameter_Homotopy_State.Set_Number_of_Parameters(nbpar);
    Assign(natural32(nbpar),b,idx);
    Parameter_Homotopy_State.Set_Indices(idx);
    return 0;
  end Job0;

  function Job1 return integer32 is -- define parameters symbolically

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    nbequ : constant integer32 := integer32(v_a(v_a'first));
    nbvar : constant integer32 := integer32(v_a(v_a'first+1));
    nbpar : constant integer32 := integer32(v_a(v_a'first+2));
    nbchr : constant natural32 := natural32(v_a(v_a'first+3));
    idx : Standard_Integer_Vectors.Vector(1..nbpar);
    v_b : constant C_Integer_Array 
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(nbchr));
    s : constant string
      := C_Integer_Array_to_String(nbchr,v_b);
    ind : integer;

  begin
    Parameter_Homotopy_State.Set_Number_of_Equations(nbequ);
    Parameter_Homotopy_State.Set_Number_of_Variables(nbvar);
    Parameter_Homotopy_State.Set_Number_of_Parameters(nbpar);
   -- put_line("The string in the interface : " & s);
    ind := s'first-1;
    for k in 1..nbpar loop
      ind := ind + 1;
      declare
        sb : Symbol_Table.Symbol;
        sx : integer := sb'first-1;
        idxsb : natural32;
      begin
        sb := (sb'range => ' ');
        while ind <= s'last loop
          exit when (s(ind) = ' ');
          sx := sx + 1;
          sb(sx) := s(ind);
          ind := ind + 1;
        end loop;
       -- put("Read symbol : "); Symbol_Table_io.put(sb); new_line;
       -- put("Number of symbols in the table : ");
       -- put(Symbol_Table.Number,1); new_line;
        idxsb := Symbol_Table.Get(sb);
       -- put("The index in the symbol table : "); put(idxsb,1); new_line;
        if idxsb > 0
         then idx(k) := integer32(idxsb);
        end if;
       -- put("the vector idx : "); put(idx); new_line;
      end;
    end loop;
    Parameter_Homotopy_State.Set_Indices(idx);
    return 0;
  end Job1;

  function Job2 return integer32 is -- return number of equations

    n : constant integer32
      :=  Parameter_Homotopy_State.Get_Number_of_Equations;
  
  begin
    Assign(n,a);
    return 0;
  end Job2;

  function Job3 return integer32 is -- return number of variables

    n : constant integer32
      :=  Parameter_Homotopy_State.Get_Number_of_Variables;
  
  begin
    Assign(n,a);
    return 0;
  end Job3;

  function Job4 return integer32 is -- return number of parameters

    n : constant integer32
      :=  Parameter_Homotopy_State.Get_Number_of_Parameters;
  
  begin
    Assign(n,a);
    return 0;
  end Job4;

  function Job5 return integer32 is -- return parameters numerically

    idx : constant Standard_Integer_Vectors.Link_to_Vector
        := Parameter_Homotopy_State.Get_Indices;

    use Standard_Integer_Vectors;

  begin
    if idx /= null
     then Assign(idx.all,a);
    end if;
    return 0;
  end Job5;

  function Job6 return integer32 is -- return parameters symbolically

    idx : constant Standard_Integer_Vectors.Link_to_Vector
        := Parameter_Homotopy_State.Get_Indices;

    use Standard_Integer_Vectors;
    use String_Splitters;

    ls : Link_to_String;
    sb : Symbol_Table.Symbol;
    len : natural32;

  begin
    if idx /= null then
      for i in idx'range loop
        sb := Symbol_Table.Get(natural32(idx(i)));
       -- put("Retrieved symbol : "); Symbol_Table_io.put(sb); new_line;
        len := Symbol_Table.Length(sb);
        declare
          strsb : string(1..integer(len+1));
        begin
          for k in 1..len loop
            strsb(integer(k)) := sb(integer(k));
          end loop;
          strsb(strsb'last) := ' ';
          Append(ls,strsb);
        end;
      end loop;
      len := natural32(ls'last) - 1; -- skip the last space
      Assign(integer32(len),a);
      if len > 0 then
        declare
          v : constant Standard_Integer_Vectors.Vector
            := String_to_Integer_Vector(ls(ls'first..integer(len)));
        begin
          Assign(v,b);
        end;
      end if;
    end if;
    return 0;
  end Job6;

  function Job7 return integer32 is -- clear parameter definitions
  begin
    Parameter_Homotopy_State.Clear;
    return 0;
  end Job7;

  function Job8 return integer32 is -- set start or target parameter values

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    precision : constant natural32 := natural32(v_a(v_a'first));
    startortarget : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    dim : constant integer32 := integer32(v_b(v_b'first));

  begin
    if precision = 0 then
      declare
        cff : Standard_Complex_Vectors.Vector(1..dim);
      begin
        Assign(natural32(2*dim),c,cff);
        if startortarget = 0 then
          Parameter_Homotopy_State.Set_Start(cff);
        elsif startortarget = 1 then
          Parameter_Homotopy_State.Set_Target(cff);
        else
          return 8;
        end if;
      end;
    elsif precision = 1 then
      declare
        cff : DoblDobl_Complex_Vectors.Vector(1..dim);
      begin
        Assign(natural32(4*dim),c,cff);
        if startortarget = 0 then
          Parameter_Homotopy_State.Set_Start(cff);
        elsif startortarget = 1 then
          Parameter_Homotopy_State.Set_Target(cff);
        else
          return 8;
        end if;
      end;
    elsif precision = 2 then
      declare
        cff : QuadDobl_Complex_Vectors.Vector(1..dim);
      begin
        Assign(natural32(8*dim),c,cff);
        if startortarget = 0 then
          Parameter_Homotopy_State.Set_Start(cff);
        elsif startortarget = 1 then
          Parameter_Homotopy_State.Set_Target(cff);
        else
          return 8;
        end if;
      end;
    else
      return 8;
    end if;
    return 0;
  end Job8;

  function Job9 return integer32 is -- get start or target parameter values

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    precision : constant natural32 := natural32(v_a(v_a'first));
    startortarget : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    dim : constant integer32 := integer32(v_b(v_b'first));

    use Standard_Complex_Vectors;
    use DoblDobl_Complex_Vectors;
    use QuadDobl_Complex_Vectors;

  begin
    if precision = 0 then
      declare
        lnkcff : Standard_Complex_Vectors.Link_to_Vector;
      begin
        if startortarget = 0 then
          lnkcff := Parameter_Homotopy_State.Get_Start;
        elsif startortarget = 1 then
          lnkcff := Parameter_Homotopy_State.Get_Target;
        else
          return 9;
        end if;
        if lnkcff /= null
         then Assign(lnkcff.all,c); 
        end if;
      end;
    elsif precision = 1 then
      declare
        lnkcff : DoblDobl_Complex_Vectors.Link_to_Vector;
      begin
        if startortarget = 0 then
          lnkcff := Parameter_Homotopy_State.Get_Start;
        elsif startortarget = 1 then
          lnkcff := Parameter_Homotopy_State.Get_Target;
        else
          return 9;
        end if;
        if lnkcff /= null
         then Assign(lnkcff.all,c); 
        end if;
      end;
    elsif precision = 2 then
      declare
        lnkcff : QuadDobl_Complex_Vectors.Link_to_Vector;
      begin
        if startortarget = 0 then
          lnkcff := Parameter_Homotopy_State.Get_Start;
        elsif startortarget = 1 then
          lnkcff := Parameter_Homotopy_State.Get_Target;
        else
          return 9;
        end if;
        if lnkcff /= null
         then Assign(lnkcff.all,c); 
        end if;
      end;
    else
      return 9;
    end if;
    return 0;
  end Job9;

  function Standard_Complex_Sweep ( x,y : double_float ) return integer32 is

  -- DESCRIPTION :
  --   Runs a complex convex parameter continuation with real and
  --   imaginary parts for gamma in x and y respectively.
  --   If both x and y are zero, then gamma will be taken equal to one.
  --   Computations happen in standard double precision.
  --   On return is the failure code, which is zero if all went well.

    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    nb_equ : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Equations;
    nb_var : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Variables;
    nb_par : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Parameters;
    indpar : constant Standard_Integer_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Indices;
    indvar : constant Standard_Integer_Vectors.Vector
           := Standard_Parameter_Systems.Complement(nb_var,indpar.all);
    nv : constant integer32 := indvar'last;
    vrsols : Standard_Complex_Solutions.Solution_List
           := Standard_Parameter_Solutions.Select_Variables(sols,nv,indvar);
    startv : Standard_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Start;
    target : Standard_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Target;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    newsols : Standard_Complex_Solutions.Solution_List;

    use Complex_Convex_Continuation;

    function Eval_Pars ( t : Standard_Complex_Numbers.Complex_Number )
                       return Standard_Complex_Vectors.Vector is
    begin
     -- return Interpolate(startv.all,target.all,t);
      return Circulate(startv.all,target.all,gamma,t);
    end Eval_Pars;

    function Diff_Pars ( t : Standard_Complex_Numbers.Complex_Number )
                       return Standard_Complex_Vectors.Vector is
    begin
      return Differentiate(startv.all,target.all);
    end Diff_Pars;

    procedure Par_Con is
      new Standard_Silent_Parameter_Continuation(Eval_Pars,Diff_Pars);

    use Standard_Complex_Solutions;

  begin
  -- put("The number of equations : "); put(nb_equ,1); new_line;
  -- put("The number of variables : "); put(nb_var,1); new_line;
  -- put("The number of parameters : "); put(nb_par,1); new_line;
  -- put("The solution list on input :");
  -- put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  -- put("After selection of the variables in the solution list on input :");
  -- put(standard_output,Length_Of(vrsols),natural32(Head_Of(vrsols).n),vrsols);
  -- Par_Con(nb_var,lp(1..nb_var-nb_par),indpar.all,indvar,vrsols);
    Par_Con(nb_var,lp.all,indpar.all,indvar,vrsols);
    newsols := Standard_Parameter_Solutions.Join_Variables
                 (vrsols,nb_var,indvar,indpar.all,target.all);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(newsols);
    return 0;
  exception
    when others => return 10;
  end Standard_Complex_Sweep;

  function DoblDobl_Complex_Sweep ( x,y : double_float ) return integer32 is

  -- DESCRIPTION :
  --   Runs a complex convex parameter continuation with real and
  --   imaginary parts for gamma in x and y respectively.
  --   If both x and y are zero, then gamma will be taken equal to one.
  --   Computations happen in double double precision.
  --   On return is the failure code, which is zero if all went well.

    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;
    nb_equ : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Equations;
    nb_var : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Variables;
    nb_par : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Parameters;
    indpar : constant Standard_Integer_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Indices;
    indvar : constant Standard_Integer_Vectors.Vector
           := Standard_Parameter_Systems.Complement(nb_var,indpar.all);
    nv : constant integer32 := indvar'last;
    vrsols : DoblDobl_Complex_Solutions.Solution_List
           := DoblDobl_Parameter_Solutions.Select_Variables(sols,nv,indvar);
    startv : DoblDobl_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Start;
    target : DoblDobl_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Target;
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;
    newsols : DoblDobl_Complex_Solutions.Solution_List;

    use Complex_Convex_Continuation;

    function Eval_Pars ( t : DoblDobl_Complex_Numbers.Complex_Number )
                       return DoblDobl_Complex_Vectors.Vector is
    begin
     -- return Interpolate(startv.all,target.all,t);
      return Circulate(startv.all,target.all,gamma,t);
    end Eval_Pars;

    function Diff_Pars ( t : DoblDobl_Complex_Numbers.Complex_Number )
                       return DoblDobl_Complex_Vectors.Vector is
    begin
      return Differentiate(startv.all,target.all);
    end Diff_Pars;

    procedure Par_Con is
      new DoblDobl_Silent_Parameter_Continuation(Eval_Pars,Diff_Pars);

  begin
    -- Par_Con(nb_var,lp(1..nb_var-nb_par),indpar.all,indvar,vrsols);
    Par_Con(nb_var,lp.all,indpar.all,indvar,vrsols);
    newsols := DoblDobl_Parameter_Solutions.Join_Variables
                 (vrsols,nb_var,indvar,indpar.all,target.all);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(newsols);
    return 0;
  exception
    when others => return 10;
  end DoblDobl_Complex_Sweep;

  function QuadDobl_Complex_Sweep ( x,y : double_float ) return integer32 is

  -- DESCRIPTION :
  --   Runs a complex convex parameter continuation with real and
  --   imaginary parts for gamma in x and y respectively.
  --   If both x and y are zero, then gamma will be taken equal to one.
  --   Computations happend in quad double precision.
  --   On return is the failure code, which is zero if all went well.

    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    sols : QuadDobl_Complex_Solutions.Solution_List
         := QuadDobl_Solutions_Container.Retrieve;
    nb_equ : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Equations;
    nb_var : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Variables;
    nb_par : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Parameters;
    indpar : constant Standard_Integer_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Indices;
    indvar : constant Standard_Integer_Vectors.Vector
           := Standard_Parameter_Systems.Complement(nb_var,indpar.all);
    nv : constant integer32 := indvar'last;
    vrsols : QuadDobl_Complex_Solutions.Solution_List
           := QuadDobl_Parameter_Solutions.Select_Variables(sols,nv,indvar);
    startv : QuadDobl_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Start;
    target : QuadDobl_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Target;
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;
    newsols : QuadDobl_Complex_Solutions.Solution_List;

    use Complex_Convex_Continuation;

    function Eval_Pars ( t : QuadDobl_Complex_Numbers.Complex_Number )
                       return QuadDobl_Complex_Vectors.Vector is
    begin
     -- return Interpolate(startv.all,target.all,t);
      return Circulate(startv.all,target.all,gamma,t);
    end Eval_Pars;

    function Diff_Pars ( t : QuadDobl_Complex_Numbers.Complex_Number )
                       return QuadDobl_Complex_Vectors.Vector is
    begin
      return Differentiate(startv.all,target.all);
    end Diff_Pars;

    procedure Par_Con is
      new QuadDobl_Silent_Parameter_Continuation(Eval_Pars,Diff_Pars);

  begin
    -- Par_Con(nb_var,lp(1..nb_var-nb_par),indpar.all,indvar,vrsols);
    Par_Con(nb_var,lp.all,indpar.all,indvar,vrsols);
    newsols := QuadDobl_Parameter_Solutions.Join_Variables
                 (vrsols,nb_var,indvar,indpar.all,target.all);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(newsols);
    return 0;
  exception
    when others => return 10;
  end QuadDobl_Complex_Sweep;

  function Job10 return integer32 is -- complex convex-parameter sweep

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    precision : constant natural32 := natural32(v_a(v_a'first));
    randgamma : constant natural32 := natural32(v_a(v_a'first+1));

  begin
   -- put("randgamma = "); put(randgamma,1); new_line;
   -- put("precision = "); put(precision,1); new_line;
    if randgamma < 2 then -- no user given gamma
      case precision is
        when 0 => return Standard_Complex_Sweep(0.0,0.0);
        when 1 => return DoblDobl_Complex_Sweep(0.0,0.0);
        when 2 => return QuadDobl_Complex_Sweep(0.0,0.0);
        when others => null;
      end case;
    else -- extract the user given gamma from c
      declare
        v_c : constant C_Double_Array
            := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2));
        regamma : constant double_float
                := double_float(v_c(v_c'first));
        imgamma : constant double_float
                := double_float(v_c(v_c'first+1));
      begin
        case precision is
          when 0 => return Standard_Complex_Sweep(regamma,imgamma);
          when 1 => return DoblDobl_Complex_Sweep(regamma,imgamma);
          when 2 => return QuadDobl_Complex_Sweep(regamma,imgamma);
          when others => null;
        end case;
      end;
    end if;
    return 0;
  exception
    when others => return 10;
  end Job10;

  function Standard_Real_Sweep return integer32 is

  -- DESCRIPTION :
  --   Runs through the container for the solution in standard double
  --   precision and runs a natural-parameter sweep defined by the
  --   family of polynomial systems in the standard systems container.
  --   The sweep runs from the start to the target values of the
  --   parameters, or until a singularity is encountered.

    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    nb_equ : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Equations;
    nb_var : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Variables;
    nb_par : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Parameters;
    indpar : constant Standard_Integer_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Indices;
    startv : Standard_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Start;
    target : Standard_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Target;
    dx : Standard_Complex_Vectors.Vector(1..nb_var);
    rx,rdx : Standard_Floating_Vectors.Vector(1..nb_var);

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Standard_Quad_Sweepers;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    isreal : constant boolean := Standard_Complex_to_Real_Poly.Is_Real(lp.all);
    realvals : boolean; -- true if everything is real

    spf : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(lp'range);
    sjm : Standard_Complex_Jaco_Matrices.Jaco_Mat(lp'range,1..nb_var);
    sjf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat(lp'range,1..nb_var);
    rp : Standard_Floating_Poly_Systems.Poly_Sys(lp'range);
    rpf : Standard_Floating_Poly_SysFun.Eval_Poly_Sys(lp'range);
    rjm : Standard_Floating_Jaco_Matrices.Jaco_Mat(lp'range,1..nb_var);
    rjf : Standard_Floating_Jaco_Matrices.Eval_Jaco_Mat(lp'range,1..nb_var);

  begin
    spf := Standard_Complex_Poly_SysFun.Create(lp.all);
    sjm := Standard_Complex_Jaco_Matrices.Create(lp.all);
    sjf := Standard_Complex_Jaco_Matrices.Create(sjm);
    if isreal then
      rp := Standard_Complex_to_Real_Poly.Convert_Complex_to_Real(lp.all);
      rpf := Standard_Floating_Poly_SysFun.Create(rp);
      rjm := Standard_Floating_Jaco_Matrices.Create(rp);
      rjf := Standard_Floating_Jaco_Matrices.Create(rjm);
    end if;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ls.t := Create(0.0);
      if not isreal then
        realvals := false;
      else
        if(IMAG_Part(startv(startv'first)) /= 0.0) then
          realvals := false;
        elsif(IMAG_Part(target(target'first)) /= 0.0) then
          realvals := false;
        else
          realvals := Standard_Solution_Diagnostics.Is_Real(ls.all,1.0E-14);
        end if;
      end if;
      if realvals then
        for i in ls.v'range loop
          rx(i) := REAL_PART(ls.v(i));
          rdx(i) := 0.0;
        end loop;
        Silent_Real_Sweep
          (true,natural32(nb_equ),natural32(nb_var),1.0,rpf,rjf,rx,rdx);
        for i in rx'range loop
          ls.v(i) := Create(rx(i));
        end loop;
      else
        for i in ls.v'range loop
          dx(i) := Create(0.0);
        end loop;
        Start_Complex_Sweep
          (natural32(nb_equ),natural32(nb_var),1.0,spf,sjf,ls.v,dx);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Complex_Poly_SysFun.Clear(spf);
    Standard_Complex_Jaco_Matrices.Clear(sjm);
    Standard_Complex_Jaco_Matrices.Clear(sjf);
    if isreal then
      Standard_Floating_Poly_Systems.Clear(rp);
      Standard_Floating_Poly_SysFun.Clear(rpf);
      Standard_Floating_Jaco_Matrices.Clear(rjm);
      Standard_Floating_Jaco_Matrices.Clear(rjf);
    end if;
    return 0;
  end Standard_Real_Sweep;

  function DoblDobl_Real_Sweep return integer32 is

  -- DESCRIPTION :
  --   Runs through the container for the solution in double double
  --   precision and runs a natural-parameter sweep defined by the
  --   family of polynomial systems in the dobldobl systems container.
  --   The sweep runs from the start to the target values of the
  --   parameters, or until a singularity is encountered.

    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;
    nb_equ : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Equations;
    nb_var : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Variables;
    nb_par : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Parameters;
    indpar : constant Standard_Integer_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Indices;
    startv : DoblDobl_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Start;
    target : DoblDobl_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Target;
    dx : DoblDobl_Complex_Vectors.Vector(1..nb_var);
    rx,rdx : Double_Double_Vectors.Vector(1..nb_var);

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Quad_Sweepers;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    isreal : constant boolean := DoblDobl_Complex_to_Real_Poly.Is_Real(lp.all);
    realvals : boolean; -- true if everything is real

    spf : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(lp'range);
    sjm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(lp'range,1..nb_var);
    sjf : DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat(lp'range,1..nb_var);
    rp : Double_Double_Poly_Systems.Poly_Sys(lp'range);
    rpf : Double_Double_Poly_SysFun.Eval_Poly_Sys(lp'range);
    rjm : Double_Double_Jaco_Matrices.Jaco_Mat(lp'range,1..nb_var);
    rjf : Double_Double_Jaco_Matrices.Eval_Jaco_Mat(lp'range,1..nb_var);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);

  begin
    spf := DoblDobl_Complex_Poly_SysFun.Create(lp.all);
    sjm := DoblDobl_Complex_Jaco_Matrices.Create(lp.all);
    sjf := DoblDobl_Complex_Jaco_Matrices.Create(sjm);
    if isreal then
      rp := DoblDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(lp.all);
      rpf := Double_Double_Poly_SysFun.Create(rp);
      rjm := Double_Double_Jaco_Matrices.Create(rp);
      rjf := Double_Double_Jaco_Matrices.Create(rjm);
    end if;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ls.t := Create(integer(0));
      if not isreal then
        realvals := false;
      else
        if(IMAG_Part(startv(startv'first)) /= zero) then
          realvals := false;
        elsif(IMAG_Part(target(target'first)) /= zero) then
          realvals := false;
        else
          realvals := DoblDobl_Solution_Diagnostics.Is_Real(ls.all,1.0E-14);
        end if;
      end if;
      if realvals then
        for i in ls.v'range loop
          rx(i) := REAL_PART(ls.v(i));
          rdx(i) := create(0.0);
        end loop;
        Silent_Real_Sweep
          (true,natural32(nb_equ),natural32(nb_var),one,rpf,rjf,rx,rdx);
        for i in rx'range loop
          ls.v(i) := Create(rx(i));
        end loop;
      else
        for i in ls.v'range loop
          dx(i) := Create(integer(0));
        end loop;
        Start_Complex_Sweep
          (natural32(nb_equ),natural32(nb_var),one,spf,sjf,ls.v,dx);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    DoblDobl_Complex_Poly_SysFun.Clear(spf);
    DoblDobl_Complex_Jaco_Matrices.Clear(sjm);
    DoblDobl_Complex_Jaco_Matrices.Clear(sjf);
    if isreal then
      Double_Double_Poly_Systems.Clear(rp);
      Double_Double_Poly_SysFun.Clear(rpf);
      Double_Double_Jaco_Matrices.Clear(rjm);
      Double_Double_Jaco_Matrices.Clear(rjf);
    end if;
    return 0;
  end DoblDobl_Real_Sweep;

  function QuadDobl_Real_Sweep return integer32 is

  -- DESCRIPTION :
  --   Runs through the container for the solution in quad double
  --   precision and runs a natural-parameter sweep defined by the
  --   family of polynomial systems in the quaddobl systems container.
  --   The sweep runs from the start to the target values of the
  --   parameters, or until a singularity is encountered.

    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    sols : QuadDobl_Complex_Solutions.Solution_List
         := QuadDobl_Solutions_Container.Retrieve;
    nb_equ : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Equations;
    nb_var : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Variables;
    nb_par : constant integer32
           := Parameter_Homotopy_State.Get_Number_of_Parameters;
    indpar : constant Standard_Integer_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Indices;
    startv : QuadDobl_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Start;
    target : QuadDobl_Complex_Vectors.Link_to_Vector
           := Parameter_Homotopy_State.Get_Target;
    dx : QuadDobl_Complex_Vectors.Vector(1..nb_var);
    rx,rdx : Quad_Double_Vectors.Vector(1..nb_var);

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Quad_Sweepers;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    isreal : constant boolean := QuadDobl_Complex_to_Real_Poly.Is_Real(lp.all);
    realvals : boolean; -- true if everything is real

    spf : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(lp'range);
    sjm : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(lp'range,1..nb_var);
    sjf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat(lp'range,1..nb_var);
    rp : Quad_Double_Poly_Systems.Poly_Sys(lp'range);
    rpf : Quad_Double_Poly_SysFun.Eval_Poly_Sys(lp'range);
    rjm : Quad_Double_Jaco_Matrices.Jaco_Mat(lp'range,1..nb_var);
    rjf : Quad_Double_Jaco_Matrices.Eval_Jaco_Mat(lp'range,1..nb_var);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    spf := QuadDobl_Complex_Poly_SysFun.Create(lp.all);
    sjm := QuadDobl_Complex_Jaco_Matrices.Create(lp.all);
    sjf := QuadDobl_Complex_Jaco_Matrices.Create(sjm);
    if isreal then
      rp := QuadDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(lp.all);
      rpf := Quad_Double_Poly_SysFun.Create(rp);
      rjm := Quad_Double_Jaco_Matrices.Create(rp);
      rjf := Quad_Double_Jaco_Matrices.Create(rjm);
    end if;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ls.t := Create(integer(0));
      if not isreal then
        realvals := false;
      else
        if(IMAG_Part(startv(startv'first)) /= zero) then
          realvals := false;
        elsif(IMAG_Part(target(target'first)) /= zero) then
          realvals := false;
        else
          realvals := QuadDobl_Solution_Diagnostics.Is_Real(ls.all,1.0E-14);
        end if;
      end if;
      if realvals then
        for i in ls.v'range loop
          rx(i) := REAL_PART(ls.v(i));
          rdx(i) := create(0.0);
        end loop;
        Silent_Real_Sweep
          (true,natural32(nb_equ),natural32(nb_var),one,rpf,rjf,rx,rdx);
        for i in rx'range loop
          ls.v(i) := Create(rx(i));
        end loop;
      else
        for i in ls.v'range loop
          dx(i) := Create(integer(0));
        end loop;
        Start_Complex_Sweep
          (natural32(nb_equ),natural32(nb_var),one,spf,sjf,ls.v,dx);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    QuadDobl_Complex_Poly_SysFun.Clear(spf);
    QuadDobl_Complex_Jaco_Matrices.Clear(sjm);
    QuadDobl_Complex_Jaco_Matrices.Clear(sjf);
    if isreal then
      Quad_Double_Poly_Systems.Clear(rp);
      Quad_Double_Poly_SysFun.Clear(rpf);
      Quad_Double_Jaco_Matrices.Clear(rjm);
      Quad_Double_Jaco_Matrices.Clear(rjf);
    end if;
    return 0;
  end QuadDobl_Real_Sweep;

  function Job11 return integer32 is -- real natural-parameter sweep

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    precision : constant natural32 := natural32(v_a(v_a'first));

  begin
    case precision is
      when 0 => return Standard_Real_Sweep;
      when 1 => return DoblDobl_Real_Sweep;
      when 2 => return QuadDobl_Real_Sweep;
      when others => null;
    end case;
    return 0;
  exception
    when others => return 11;
  end Job11;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0;   -- define parameters numerically
      when 1 => return Job1;   -- define parameters symbolically
      when 2 => return Job2;   -- return the number of equations
      when 3 => return Job3;   -- return the number of variables
      when 4 => return Job4;   -- return the number of parameters
      when 5 => return Job5;   -- return parameters numerically
      when 6 => return Job6;   -- return parameters symbolically
      when 7 => return Job7;   -- clear parameter definitions
      when 8 => return Job8;   -- set start or target parameter values
      when 9 => return Job9;   -- get start or target parameter values
      when 10 => return Job10; -- complex convex-parameter sweep
      when 11 => return Job11; -- real natural-parameter sweep
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_sweep;
