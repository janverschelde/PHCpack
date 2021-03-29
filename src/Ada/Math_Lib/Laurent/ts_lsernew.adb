with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Integer_Matrices;
with Standard_Complex_VecVecVecs;
with Symbol_Table;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_JacoMats;    use Standard_Complex_Laur_JacoMats;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Double_Laurent_Series;
with Double_Linear_Laurent_Solvers;     use Double_Linear_Laurent_Solvers;
with Double_Lseries_Polynomials;        use Double_Lseries_Polynomials;
with Test_Double_Lseries_Matrices;

procedure ts_lsernew is

-- DESCRIPTION :
--   Tests the development of Newton's methon on Laurent series.

  procedure Add_Parameter ( tv : in Table_Vector ) is

  -- DESCRIPTION :
  --   To the first constant term in tv, adds t.
  --   Writes a warning if there is no constant term.

  -- REQUIRED :
  --   The degree in the table vector must be at least one.

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

  procedure Make_Series
              ( sol : in Standard_Complex_Vectors.Vector;
                deg : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Returns a regular power series of degree deg,
  --   with as leading coefficients the coordinates in sol.
  --   Allocates all space for cffs.

  -- REQUIRED : 
  --   lead'range = sol'range.

    xcffs : Standard_Complex_VecVecs.VecVec(sol'range);

  begin
    lead := (sol'range => 0);
    for i in sol'range loop
      declare
        cff : Standard_Complex_Vectors.Vector(0..deg);
      begin
        cff(0) := sol(i);
        for k in 1..deg loop
          cff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        xcffs(i) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    cffs := new Standard_Complex_VecVecs.VecVec'(xcffs);
  end Make_Series;

  procedure Make_Series
              ( sol : in Laur_Sys; deg : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Returns a regular power series of degree deg,
  --   with as leading terms the coordinates in sol.
  --   Allocates all space for cffs.

  -- REQUIRED : 
  --   lead'range = sol'range.

    xcffs : Standard_Complex_VecVecs.VecVec(sol'range);
    dg : constant Standard_Integer_Vectors.Vector(1..1) := (1..1 => 0);
    ldg : Standard_Complex_Laurentials.Degrees
        := new Standard_Integer_Vectors.Vector'(dg);

  begin
    lead := (sol'range => 0);
    for i in sol'range loop
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(sol(i))),1); new_line;
      declare
        cff : Standard_Complex_Vectors.Vector(0..deg);
      begin
        for k in 0..deg loop
          ldg(1) := k;
          cff(k) := Standard_Complex_Laurentials.Coeff(sol(i),ldg);
        end loop;
        xcffs(i) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    cffs := new Standard_Complex_VecVecs.VecVec'(xcffs);
    Standard_Complex_Laurentials.Clear(ldg);
  end Make_Series;

  procedure Newton_Step
              ( deg : in integer32;
                p : in Table_Vector; jp : in Table_Vector_Array;
                xlead : in out Standard_Integer_Vectors.Vector;
                xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                ylead : in out Standard_Integer_Vectors.Vector;
                ycffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                Alead : in out Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : in out Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                dxlead : in out Standard_Integer_Vectors.Vector;
                dxcffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                rlead : in out Standard_Integer_Vectors.Vector;
                rcffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Does one step with Newton's method.

  -- REQUIRED : ycffs, dxcffs, Acffs, and Bcffs are allocated.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   p        table representation of a series system;
  --   jp       table representation of the Jacobian matrix;
  --   xlead    leading exponents of the current solution series;
  --   xcffs    coefficient vector of the current solution series;
  --   verbose  flag for the verbosity.

  -- ON RETURN :
  --   dxlead   leading exponents of the update to the solution series;
  --   dxcffs   coefficient vector of the update to the solution series;
  --   Alead    leading exponents of the LU factorization
  --   Acffs    factors of the Jacobian matrix evaluated at (xlead, xcffs);
  --   Blead    leading exponents of the evaluated Jacobian matrix;
  --   Bcffs    Jacobian matrix evaluated at (xlead, xcffs).

    pivots : Standard_Integer_Vectors.Vector(xlead'range);
    acff,bcff : Standard_Complex_Vectors.Link_to_Vector;

    use Standard_Complex_Numbers;

  begin
    put_line("Evaluating the table vector ...");
    Eval(deg,p,xlead,xcffs,ylead,ycffs.all);
    if verbose
     then Test_Double_Lseries_Matrices.Write(ylead,ycffs,"y");
    end if;
    put_line("Evaluating the table vector array ...");
    Eval(deg,jp,xlead,xcffs,Alead,Acffs);
    if verbose then
      Test_Double_Lseries_Matrices.Copy
        (p.nbt,p.nbt,deg,Alead,Acffs,Blead,Bcffs);
      Test_Double_Lseries_Matrices.Write(Blead,Bcffs,"B");
    end if;
    LU_Factorization(p.nbt,p.nbt,deg,Alead,Acffs,pivots);
    if verbose
     then put("The pivots : "); put(pivots); new_line;
    end if;
    for i in pivots'range loop
      dxlead(i) := 0;
      acff := ycffs(pivots(i));
      bcff := dxcffs(i);
      for k in 0..deg loop
        bcff(k) := -acff(k);
      end loop;
    end loop;
    Test_Double_Lseries_Matrices.Write(dxlead,dxcffs,"b");
    Forward_Substitution(deg,Alead,Acffs,dxlead,dxcffs,ylead,ycffs);
    Backward_Substitution(deg,Alead,Acffs,ylead,ycffs,dxlead,dxcffs);
    if verbose then
      Test_Double_Lseries_Matrices.Write(dxlead,dxcffs,"dx");
      Matrix_Vector_Product(deg,Blead,Bcffs,dxlead,dxcffs,rlead,rcffs);
      Test_Double_Lseries_Matrices.Write(rlead,rcffs,"r");
    end if;
    for i in dxlead'range loop
      Double_Laurent_Series.Add(deg,
        xlead(i),dxlead(i),xcffs(i).all,dxcffs(i).all,ylead(i),ycffs(i).all);
      xlead(i) := ylead(i);
      acff := ycffs(i);
      bcff := xcffs(i);
      for k in 0..deg loop
        bcff(k) := acff(k);
      end loop;
    end loop;
  end Newton_Step;

  procedure Test_Regular_Newton
              ( p : in Laur_Sys; sol : in Standard_Complex_Vectors.Vector;
                deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the system p,
  --   starting at the leading constants of a regular solution,
  --   on Laurent series where the highest power of t equals deg. 

    neq : constant integer32 := p'last;
    dim : constant integer32 := neq;
    nvr : constant integer32 := sol'last;
    tdx : constant integer32 := 0;
    tv : constant Table_Vector(neq) := Make_Table_Vector(p,dim,nvr,tdx,deg);
    jp : constant Jaco_Mat(1..neq,1..dim) := Create(p);
    tva : constant Table_Vector_Array(1..neq)
        := Make_Table_Vector_Array(jp,tdx,deg);
    xlead,ylead,dxlead,rlead : Standard_Integer_Vectors.Vector(1..nvr);
    xcffs,ycffs,dxcffs,rcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    Alead,Blead : Standard_Integer_Matrices.Matrix(1..neq,1..nvr);
    Acffs,Bcffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    ans : character;

  begin
    put_line("The table representation :"); Write(tv);
    Add_Parameter(tv);
    put_line("After adding a parameter :"); Write(tv);
    Make_Series(sol,deg,xlead,xcffs);
    put("A "); put(nvr,1); put_line("-vector of Laurent series :");
    Test_Double_Lseries_Matrices.Write(xlead,xcffs,"x");
    Test_Double_Lseries_Matrices.Allocate_Series_Coefficients(nvr,deg,ycffs);
    Test_Double_Lseries_Matrices.Allocate_Series_Coefficients(nvr,deg,dxcffs);
    Test_Double_Lseries_Matrices.Allocate_Series_Coefficients(nvr,deg,rcffs);
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,neq,1,nvr,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Bcffs,1,neq,1,nvr,0,deg);
    loop
      Newton_Step(deg,tv,tva,xlead,xcffs,ylead,ycffs,
                  Alead,Acffs,Blead,Bcffs,dxlead,dxcffs,rlead,rcffs);
      Test_Double_Lseries_Matrices.Write(xlead,xcffs,"x");
      put("Do another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Regular_Newton;

  procedure Test_Singular_Newton
              ( p,sol : in Laur_Sys; tdx,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the system p,
  --   starting at the series defined in sol,
  --   on Laurent series where the highest power of t equals deg. 

  -- REQUIRED : tdx /= 0 and 
  --   p has one more variable than the number of equations.

  -- ON ENTRY :
  --   p       a system with one parameter t;
  --   sol     as many univariate polynomials in t as p'length;
  --   tdx     the index of t in p;
  --   deg     precision of the series.

    neq : constant integer32 := p'last;
    dim : constant integer32 := neq + 1;
    nvr : constant integer32 := neq; -- number of variables minus t
    tv : constant Table_Vector(neq) := Make_Table_Vector(p,dim,nvr,tdx,deg);
    jp : constant Jaco_Mat(1..neq,1..dim) := Create(p);
    tva : constant Table_Vector_Array(1..neq)
        := Make_Table_Vector_Array(jp,tdx,deg);
    xlead,ylead,dxlead,rlead : Standard_Integer_Vectors.Vector(1..nvr);
    xcffs,ycffs,dxcffs,rcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    Alead,Blead : Standard_Integer_Matrices.Matrix(1..neq,1..nvr);
    Acffs,Bcffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    ans : character;

  begin
    put_line("The table representation :"); Write(tv);
    Make_Series(sol,deg,xlead,xcffs);
    put("A "); put(nvr,1); put_line("-vector of Laurent series :");
    Test_Double_Lseries_Matrices.Write(xlead,xcffs,"x");
    Test_Double_Lseries_Matrices.Allocate_Series_Coefficients(nvr,deg,ycffs);
    Test_Double_Lseries_Matrices.Allocate_Series_Coefficients(nvr,deg,dxcffs);
    Test_Double_Lseries_Matrices.Allocate_Series_Coefficients(nvr,deg,rcffs);
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,neq,1,nvr,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Bcffs,1,neq,1,nvr,0,deg);
    loop
      Newton_Step(deg,tv,tva,xlead,xcffs,ylead,ycffs,
                  Alead,Acffs,Blead,Bcffs,dxlead,dxcffs,rlead,rcffs);
      Test_Double_Lseries_Matrices.Write(xlead,xcffs,"x");
      put("Do another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Singular_Newton;

  procedure Test_Isolated_Start is

  -- DESCRIPTION :
  --   Prompts for a system with start solutions and a degree of t.

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

  -- DESCRIPTION :
  --   Prompts for a system, a series, and a degree of t.

    lp,lsol : Link_to_Laur_Sys;
    neq,nvr,dim,tdx,deg : integer32 := 0;

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
      Test_Singular_Newton(lp.all,lsol.all,tdx,deg);
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

begin
  Main;
end ts_lsernew;
