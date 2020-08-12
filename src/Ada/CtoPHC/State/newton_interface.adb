with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;
with Standard_Complex_Poly_Systems;
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

package body Newton_Interface is

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
      put_line("Newton_Standard_Laurent_Step ...");
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
        if nbequ = nbvar then
          Silent_Root_Refiner(lp.all,work,epsxa,epsfa,tolsi,nit,1,deflate);
        else
          Silent_Root_Sharpener(lp.all,work,epsxa,epsfa,tolsi,nit,1,deflate);
        end if;
        Standard_Solutions_Container.Clear;
        Standard_Solutions_Container.Initialize(work);
      exception
        when others => return 199;
      end;
      return 0;
    end if;
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
    if lp = null or Is_Null(sols) then
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
          DoblDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        DoblDobl_Solutions_Container.Clear;
        DoblDobl_Solutions_Container.Initialize(work);
      exception
        when others => return 198;
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
          QuadDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        QuadDobl_Solutions_Container.Clear;
        QuadDobl_Solutions_Container.Initialize(work);
      exception
        when others => return 197;
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in newton_interface.");
        put_line("Newton_Standard_Polynomial_Step.");
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

end Newton_Interface;
