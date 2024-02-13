with text_io;                           use text_io;
with Interfaces.C;
with String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;
with Characters_and_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_SysFun;
with Multprec_Complex_Laur_JacoMats;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;
with Standard_Root_Refiners;
with DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;
with Multprec_Root_Refiners;
with Verification_of_Solutions;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with QuadDobl_Solutions_Container;
with Multprec_PolySys_Container;
with Multprec_LaurSys_Container;
with Multprec_Solutions_Container;
with PHCpack_Operations;
with Standard_Systems_Pool;
with Solutions_Pool;

package body Newton_Interface is

  function Newton_Standard_Polynomial_Verify
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;
    use Standard_Root_Refiners;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    work : Solution_List;
   -- tmp : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_Standard_Polynomial_Verify ...");
    end if;
   -- put_line("refining the roots ...");
   -- put(lp.all);
    if lp = null or Is_Null(sols) then
      return 9;             
    else
    -- put_line("Solutions before the refine :");
    -- put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
      Copy(sols,work);
      declare
        epsxa : constant double_float := 1.0E-12;
        epsfa : constant double_float := 1.0E-12;
        tolsi : constant double_float := 1.0E-8;
        deflate : boolean := false;
        nit : natural32 := 0;
      begin
        if PHCpack_Operations.Is_File_Defined then
          put(PHCpack_Operations.output_file,natural32(lp'last),lp.all);
          Reporting_Root_Refiner
            (PHCpack_Operations.output_file,
             lp.all,work,epsxa,epsfa,tolsi,nit,5,deflate,false);
        else
          Silent_Root_Refiner
            (lp.all,work,epsxa,epsfa,tolsi,nit,5,deflate);
        end if;
        Standard_Solutions_Container.Clear;
        Standard_Solutions_Container.Initialize(work);
       -- put("Refiner Results :");
       -- tmp := work;
       -- while not Is_Null(tmp) loop
       --   put(" ");
       --   put(Head_Of(tmp).m,1);
       --   tmp := Tail_Of(tmp);
       -- end loop;
       -- new_line;
      end;
      return 0;
     -- sols := Standard_Solutions_Container.Retrieve;
     -- put_line("Solutions after the refine :");
     -- put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Standard_Polynomial_Verify.");
      end if;
      return 9;
  end Newton_Standard_Polynomial_Verify;

  function Newton_Standard_Polynomial_Refine
             ( b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Root_Refiners;

    ls : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    ep : constant Link_to_Eval_Poly_Sys
       := Standard_PolySys_Container.Evaluator;
    jf : constant Link_to_Eval_Jaco_Mat
       := Standard_PolySys_Container.Jacobian_Evaluator;
    sol : Standard_Complex_Solutions.Solution(ls'last)
        := Convert_to_Solution(b,c);
    epsxa : constant double_float := 1.0E-14;
    epsfa : constant double_float := 1.0E-14;
    max : constant natural32 := 3;
    numit : natural32 := 0;
    fail : boolean;
   -- x : Standard_Complex_Vectors.Vector(ls'range)
   --   := (ls'range => Create(1.0));
   -- y : Standard_Complex_Vectors.Vector(ls'range);
   -- r : constant integer := Standard_Random_Numbers.Random(0,1000);

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_Standard_Polynomial_Refine ...");
    end if;
   -- put("starting evaluation with id = "); put(r,1); new_line;
   -- for i in 1..1000 loop
   --   y := Eval(ls.all,x);
   --   -- y := Eval(ls.all,sol.v); --y := Eval(ep.all,sol.v);
   -- end loop;
   -- put("ending evaluation with id = "); put(r,1); new_line;
    Silent_Newton(ep.all,jf.all,sol,epsxa,epsfa,numit,max,fail);
    Assign_Solution(sol,b,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Standard_Polynomial_Refine.");
      end if;
      return 149;
  end Newton_Standard_Polynomial_Refine;

  function Newton_Standard_Polynomial_Step
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;
    use Standard_Root_Refiners;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_Standard_Polynomial_Step ...");
    end if;
    if lp = null or Is_Null(sols) then
      return 199;             
    else
      Copy(sols,work);
      declare
        nbequ : constant integer32 := lp'last;
        nbvar : constant integer32 := Head_Of(sols).n;
        epsxa : constant double_float := 1.0E-12;
        epsfa : constant double_float := 1.0E-12;
        tolsi : constant double_float := 1.0E-8;
        deflate : boolean := false;
        nit : natural32 := 0;
      begin
        if vrblvl > 0 then
          put("the number of equations : "); put(nbequ,1); new_line;
          put("the number of variables : "); put(nbvar,1); new_line;
        end if;
        if nbequ = nbvar then
          Silent_Root_Refiner
            (lp.all,work,epsxa,epsfa,tolsi,nit,1,deflate,vrblvl-1);
        else
          Silent_Root_Sharpener
            (lp.all,work,epsxa,epsfa,tolsi,nit,1,deflate); --,vrblvl-1);
        end if;
        Standard_Solutions_Container.Clear;
        Standard_Solutions_Container.Initialize(work);
      exception
        when others => return 199;
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Standard_Polynomial_Step.");
      end if;
      return 199;
  end Newton_Standard_Polynomial_Step;

  function Newton_DoblDobl_Polynomial_Step
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Root_Refiners;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_DoblDobl_Polynomial_Step ...");
    end if;
    if lp = null then
      if vrblvl > 0
       then put_line("No polynomial system in dobldobl container?");
      end if;
      return 198;
    end if;
    if Is_Null(sols) then
      if vrblvl > 0
       then put_line("No solutions in dobldobl container?");
      end if;
      return 198;             
    else
      Copy(sols,work);
      declare
        dim : constant integer32 := Head_Of(sols).n;
        f : Eval_Poly_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,1..dim) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,1..dim) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          DoblDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res,vrblvl-1);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        DoblDobl_Solutions_Container.Clear;
        DoblDobl_Solutions_Container.Initialize(work);
      exception
        when others => 
          if vrblvl > 0 then
            put("Exception raised in newton_interface.");
            put_line("Newton_DoblDobl_Polynomial_Step.");
          end if;
          return 198;
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_DoblDobl_Polynomial_Step.");
      end if;
      return 198;
  end Newton_DoblDobl_Polynomial_Step;

  function Newton_QuadDobl_Polynomial_Step
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Root_Refiners;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_QuadDobl_Polynomial_Step ...");
    end if;
    if lp = null or Is_Null(sols) then
      return 197;             
    else
      Copy(sols,work);
      declare
        dim : constant integer32 := Head_Of(sols).n;
        f : Eval_Poly_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,1..dim) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,1..dim) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          QuadDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res,vrblvl-1);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        QuadDobl_Solutions_Container.Clear;
        QuadDobl_Solutions_Container.Initialize(work);
      exception
        when others => 
          if vrblvl > 0 then
            put("Exception raised in newton_interface.");
            put_line("Newton_QuadDobl_Polynomial_Step.");
          end if;
          return 197;
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_QuadDobl_Polynomial_Step.");
      end if;
      return 197;
  end Newton_QuadDobl_Polynomial_Step;

  function Newton_Multprec_Polynomial_Step
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Complex_Solutions;
    use Multprec_Root_Refiners;

    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;
    work : Solution_List;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    deci : constant natural32 := natural32(v_a(v_a'first));
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(deci);

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_Multprec_Polynomial_Step ...");
    end if;
    if lp = null or Is_Null(sols) then
      return 195;             
    else
      Copy(sols,work);
      Set_Size(work,size);
      declare
        f : Eval_Poly_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,lp'range) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,lp'range) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          Multprec_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        Multprec_Solutions_Container.Clear;
        Multprec_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Multprec_Polynomial_Step.");
      end if;
      return 195;
  end Newton_Multprec_Polynomial_Step;

  function Newton_Standard_Laurent_Step
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;
    use Standard_Root_Refiners;

    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    work,refsols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_Standard_Laurent_Step ...");
    end if;
    if lp = null or Is_Null(sols) then
      return 326;             
    else
      Copy(sols,work); -- must work on a copy!
      declare
        nbequ : constant integer32 := lp'last;
        nbvar : constant integer32 := Head_Of(sols).n;
        epsxa : constant double_float := 1.0E-12;
        epsfa : constant double_float := 1.0E-12;
        tolsi : constant double_float := 1.0E-8;
        nit : natural32 := 0;
      begin
        if nbequ = nbvar then
          Silent_Root_Refiner(lp.all,work,refsols,epsxa,epsfa,tolsi,nit,1);
        else
          Silent_Root_Sharpener(lp.all,work,refsols,epsxa,epsfa,tolsi,nit,1);
        end if;
        Standard_Solutions_Container.Clear;
        Standard_Solutions_Container.Initialize(work);
      exception
        when others => return 326;
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Standard_Laurent_Step.");
      end if;
      return 326;
  end Newton_Standard_Laurent_Step;

  function Newton_DoblDobl_Laurent_Step
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Root_Refiners;

    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_DoblDobl_Laurent_Step ...");
    end if;
    if lp = null or Is_Null(sols) then
      return 327;             
    else
      Copy(sols,work);
      declare
        dim : constant integer32 := Head_Of(sols).n;
        f : Eval_Laur_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,1..dim) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,1..dim) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          DoblDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        DoblDobl_Solutions_Container.Clear;
        DoblDobl_Solutions_Container.Initialize(work);
      exception
        when others => return 327;
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Standard_Laurent_Step.");
      end if;
      return 327;
  end Newton_DoblDobl_Laurent_Step;

  function Newton_QuadDobl_Laurent_Step
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Root_Refiners;

    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_QuadDobl_Laurent_Step ...");
    end if;
    if lp = null or Is_Null(sols) then
      return 328;             
    else
      Copy(sols,work);
      declare
        dim : constant integer32 := Head_Of(sols).n;
        f : Eval_Laur_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,1..dim) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,1..dim) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          QuadDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        QuadDobl_Solutions_Container.Clear;
        QuadDobl_Solutions_Container.Initialize(work);
      exception
        when others => return 328;
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_QuadDobl_Laurent_Step.");
      end if;
      return 328;
  end Newton_QuadDobl_Laurent_Step;

  function Newton_Multprec_Laurent_Step
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Laur_Systems;
    use Multprec_Complex_Laur_SysFun;
    use Multprec_Complex_Laur_JacoMats;
    use Multprec_Complex_Solutions;
    use Multprec_Root_Refiners;

    lp : constant Link_to_Laur_Sys := Multprec_LaurSys_Container.Retrieve;
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;
    work : Solution_List;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    deci : constant natural32 := natural32(v_a(v_a'first));
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(deci);

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_Multprec_Laurent_Step ...");
    end if;
    if lp = null or Is_Null(sols) then
      return 329;             
    else
      Copy(sols,work);
      Set_Size(work,size);
      declare
        f : Eval_Laur_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,lp'range) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,lp'range) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          Multprec_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        Multprec_Solutions_Container.Clear;
        Multprec_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Multprec_Laurent_Step.");
      end if;
      return 329;
  end Newton_Multprec_Laurent_Step;

  function Newton_Varbprec_Step
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use String_Splitters;
    use Verification_of_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(5));
    dim : constant natural32 := natural32(v_a(v_a'first));
    ns : constant natural32 := natural32(v_a(v_a'first+1));
    wanted : constant natural32 := natural32(v_a(v_a'first+2));
    maxitr : constant natural32 := natural32(v_a(v_a'first+3));
    maxprc : constant natural32 := natural32(v_a(v_a'first+4));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(ns-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(ns));
    s : constant String(1..integer(ns)) := C_Integer_Array_to_String(ns,v_b);
    ls : Array_of_Strings(1..integer(dim)) := Split(natural(dim),s,';');
    sols : constant Multprec_Complex_Solutions.Solution_List
         := Multprec_Solutions_Container.Retrieve;
    work : Multprec_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_Varbprec_Step ...");
    end if;
   -- put("The dimension : "); put(integer32(dim),1); new_line;
   -- put("The number of characters in the system : ");
   -- put(integer32(ns),1); new_line;
   -- put("The number of wanted number of decimal places : ");
   -- put(integer32(wanted),1); new_line;
   -- put("The maximum number of Newton steps : ");
   -- put(integer32(maxitr),1); new_line;
   -- put("The maximum number of decimal places to estimate the loss : ");
   -- put(integer32(maxprc),1); new_line;
   -- put_line("The polynomials : "); put_line(s);
   -- for i in ls'range loop
   --   put("Polynomial "); put(integer32(i),1); put_line(" :");
   --   put_line(ls(i).all);
   -- end loop;
    Multprec_Complex_Solutions.Copy(sols,work);
    Verify_Solutions_of_Laurent_Polynomials
      (ls,work,wanted,maxitr,maxprc);
    Multprec_Solutions_Container.Clear;
    Multprec_Solutions_Container.Initialize(work);
    Clear(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Varbprec_Step.");
      end if;
      return 179;
  end Newton_Varbprec_Step;

  function Newton_Standard_SysPool_Refine
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Characters_and_Numbers;
    use Standard_Complex_Poly_Sysfun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Root_Refiners;
    use Solutions_Pool;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
       -- := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant integer32 := integer32(v_a(v_a'first));
   -- n : constant integer32 := integer32(v_a(v_a'first+1));
    f : constant Link_to_Eval_Poly_Sys := Standard_Systems_Pool.Evaluator(k);
    jf : constant Link_to_Eval_Jaco_Mat
       := Standard_Systems_Pool.Jacobian_Evaluator(k);
    sols : constant Solution_List := Solutions_Pool.Retrieve(k);
    len : constant natural32 := Solutions_Pool.Length(k);
    tmp : Solution_List;
   -- sol : Solution(n) := Convert_to_Solution(b,c);
    ls : Link_to_Solution; -- := Convert_to_Solution(b,c);
    epsxa : constant double_float := 1.0E-14;
    epsfa : constant double_float := 1.0E-14;
    max : constant natural32 := 3;
    cnt : natural32 := 0;
    numit : natural32 := 0;
    fail : boolean;
   -- x : Standard_Complex_Vectors.Vector(1..n) := (1..n => Create(1.0));
   -- y : Standard_Complex_Vectors.Vector(f'range);

  begin
    if vrblvl > 0 then
      put("-> in newton_interface.");
      put_line("Newton_Standard_SysPool_Refine ...");
    end if;
   -- GNAT.Float_Control.Reset;
    put_line("Thread " & convert(k) & " starts refining " 
                       & convert(integer32(len)) & " solutions.");
    tmp := sols;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp); cnt := cnt + 1;
      put_line("Thread " & convert(k) & " refines solution "
                         & convert(integer32(cnt)));
      Silent_Newton(f.all,jf.all,ls.all,epsxa,epsfa,numit,max,fail);
     -- y := Eval(f.all,ls.v);
      tmp := Tail_Of(tmp);
    end loop;
   -- Assign_Solution(sol,b,c);
   -- Assign_Solution(ls,b,c);
    put_line(" done");
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Standard_SysPool_Refine.");
      end if;
      return 305;
  end Newton_Standard_SysPool_Refine;

end Newton_Interface;
