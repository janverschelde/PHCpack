with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Integer_Matrices;
with Standard_Complex_VecVecVecs;
with Symbol_Table;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_JacoMats;    use Standard_Complex_Laur_JacoMats;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Double_Linear_Laurent_Solvers;
with Double_Lseries_Newton_Steps;

package body Test_Double_Lseries_Newton is

  procedure Add_Parameter ( tv : in Table_Vector ) is

    eqmons : Standard_Integer_VecVecs.Link_to_VecVec;
    eqcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    mon : Standard_Integer_Vectors.Link_to_Vector;
    cff : Standard_Complex_Vectors.Link_to_Vector;
    allzero : boolean := false;

  begin
    for i in tv.mons'range loop -- search for the constant term
      eqmons := tv.mons(i);     -- monomials of i-th equation
      for j in eqmons'range loop
        mon := eqmons(j);
        allzero := true;
        for k in mon'range loop
          if mon(k) /= 0
           then allzero := false;
          end if;
          exit when not allzero;
        end loop;
        if allzero then
          eqcffs := tv.cffs(i);
          cff := eqcffs(j);
          cff(1) := Standard_Complex_Numbers.Create(1.0);
        end if;
        exit when allzero;
      end loop;
      exit when allzero;
    end loop;
    if not allzero
     then put_line("Warning: no constant term found.");
    end if;
  end Add_Parameter;

  procedure Interactive_Newton_Steps
              ( deg : in integer32;
                tv : in Table_Vector; tva : in Table_Vector_Array;
                xlead : in out Standard_Integer_Vectors.Vector;
                xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := true ) is
             
    neq : constant integer32 := tva'last;
    nvr : constant integer32 := xlead'last;
    ylead,dxlead,rlead : Standard_Integer_Vectors.Vector(1..nvr);
    ycffs,dxcffs,rcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    Alead,Blead : Standard_Integer_Matrices.Matrix(1..neq,1..nvr);
    Acffs,Bcffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    ans : character;
    stepcnt : integer32 := 1;
    dxnrm,pxnrm : double_float;

  begin
    Double_Linear_Laurent_Solvers.Allocate_Series_Coefficients(nvr,deg,ycffs);
    Double_Linear_Laurent_Solvers.Allocate_Series_Coefficients(nvr,deg,dxcffs);
    Double_Linear_Laurent_Solvers.Allocate_Series_Coefficients(nvr,deg,rcffs);
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,neq,1,nvr,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Bcffs,1,neq,1,nvr,0,deg);
    loop
      put("Step "); put(stepcnt,1); put_line(" ...");
      Double_Lseries_Newton_Steps.Newton_Step
        (deg,tv,tva,xlead,xcffs,ylead,ycffs,Alead,Acffs,Blead,Bcffs,
         dxlead,dxcffs,rlead,rcffs,dxnrm,pxnrm,verbose);
      Double_Linear_Laurent_Solvers.Write(xlead,xcffs,"x");
      put("Do another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      stepcnt := stepcnt + 1;
    end loop;
    Standard_Complex_VecVecs.Deep_Clear(ycffs);
    Standard_Complex_VecVecs.Deep_Clear(dxcffs);
    Standard_Complex_VecVecs.Deep_Clear(rcffs);
    Standard_Complex_VecVecVecs.Clear(Acffs);
    Standard_Complex_VecVecVecs.Clear(Bcffs);
  end Interactive_Newton_Steps;

  procedure Test_Regular_Newton
              ( p : in Laur_Sys; sol : in Standard_Complex_Vectors.Vector;
                deg : in integer32; verbose : in boolean := true ) is

    neq : constant integer32 := p'last;
    dim : constant integer32 := neq;
    nvr : constant integer32 := sol'last;
    tdx : constant integer32 := 0;
    tv : constant Table_Vector(neq) := Make_Table_Vector(p,dim,nvr,tdx,deg);
    jp : constant Jaco_Mat(1..neq,1..dim) := Create(p);
    tva : constant Table_Vector_Array(1..neq)
        := Make_Table_Vector_Array(jp,tdx,deg);
    xlead : Standard_Integer_Vectors.Vector(1..nvr);
    xcffs : Standard_Complex_VecVecs.Link_to_VecVec;

  begin
    put_line("The table representation :"); Write(tv);
    Add_Parameter(tv);
    put_line("After adding a parameter :"); Write(tv);
    Double_Lseries_Newton_Steps.Make_Series(sol,deg,xlead,xcffs);
    Double_Lseries_Newton_Steps.Set_Leading_Exponents(xlead);
    put("A "); put(nvr,1); put_line("-vector of Laurent series :");
    Double_Linear_Laurent_Solvers.Write(xlead,xcffs,"x");
    Interactive_Newton_Steps(deg,tv,tva,xlead,xcffs,verbose);
  end Test_Regular_Newton;

  procedure Test_Singular_Newton
              ( p,sol : in Laur_Sys; tdx,deg : in integer32;
                verbose : in boolean := true ) is

    neq : constant integer32 := p'last;
    dim : constant integer32 := neq + 1;
    nvr : constant integer32 := neq; -- number of variables minus t
    tv : constant Table_Vector(neq) := Make_Table_Vector(p,dim,nvr,tdx,deg);
    jp : constant Jaco_Mat(1..neq,1..dim) := Create(p);
    tva : constant Table_Vector_Array(1..neq)
        := Make_Table_Vector_Array(jp,tdx,deg);
    xlead : Standard_Integer_Vectors.Vector(1..nvr);
    xcffs : Standard_Complex_VecVecs.Link_to_VecVec;

  begin
    put_line("The table representation :"); Write(tv);
    Double_Lseries_Newton_Steps.Make_Series(sol,deg,xlead,xcffs);
    Double_Lseries_Newton_Steps.Set_Leading_Exponents(xlead);
    put("A "); put(nvr,1); put_line("-vector of Laurent series :");
    Double_Linear_Laurent_Solvers.Write(xlead,xcffs,"x");
    Interactive_Newton_Steps(deg,tv,tva,xlead,xcffs,verbose);
  end Test_Singular_Newton;

  procedure Test_Isolated_Start is

    lp : Link_to_Laur_Sys;
    sols : Solution_List;
    ls : Link_to_Solution;
    neq,nvr,nbsols,dim,deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    Standard_System_and_Solutions_io.get(lp,sols);
    neq := lp'last;
    nvr := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Read "); put(neq,1); put(" polynomials in ");
    put(nvr,1); put_line(" variables ...");
    nbsols := integer32(Length_Of(sols));
    put("Read "); put(nbsols,1); put_line(" solutions.");
    if nbsols > 0 then
      ls := Head_Of(sols);
      dim := ls.n;
      if dim = nvr then
        put_line("Solution dimension matches the number of variables.");
        new_line;
        put("Give the largest degree of the series parameter t : ");
        get(deg);
        Test_Regular_Newton(lp.all,ls.v,deg);
      else
        put("Solution dimension : "); put(dim,1);
        put(" != "); put(nvr,1); put_line(".");
      end if;
    end if;
  end Test_Isolated_Start;

  procedure Test_Series_Start is

    lp,lsol : Link_to_Laur_Sys;
    neq,nvr,dim,tdx,deg : integer32 := 0;
    ans : character;
    verbose : boolean;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    neq := lp'last;
    nvr := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Read "); put(neq,1); put(" polynomials in ");
    put(nvr,1); put_line(" variables ...");
    tdx := tsymbol_Index;
    put("-> index of t : "); put(tdx,1); new_line;
    if tdx /= 0 then
      new_line;
      put_line("Reading initial terms of a series ...");
      Symbol_Table.Clear;
      get(lsol);
      dim := lsol'last; 
      new_line;
      put("Read "); put(dim,1); put_line(" polynomials ...");
      deg := Degree(lsol(lsol'first));
      put("The degree of the first polynomial : "); put(deg,1); new_line;
      new_line;
      put("Give the degree : "); get(deg);
      new_line;
      put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
      verbose := (ans = 'y');
      Test_Singular_Newton(lp.all,lsol.all,tdx,deg,verbose);
    end if;
  end Test_Series_Start;
 
  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing Laurent series expansions ...");
    put_line("  1. start at a system with isolated solutions; or");
    put_line("  2. give some leading terms of a power series.");
    put("Type 1 or 2 to select a test : "); 
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Test_Isolated_Start;
      when '2' => Test_Series_Start;
      when others => null;
    end case;
  end Main;

end Test_Double_Lseries_Newton;
