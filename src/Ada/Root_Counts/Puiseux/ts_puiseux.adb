with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vector_Norms;      use Standard_Complex_Vector_Norms;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Integer32_Transformations; use Standard_Integer32_Transformations;
with Standard_Integer_Transformations_io;
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Drivers_for_Static_Lifting;         use Drivers_for_Static_Lifting;
with Black_Box_Solvers;                  use Black_Box_Solvers;

procedure ts_puiseux is

-- DESCRIPTION :
--   Development of the Newton-Puiseux algorithm.

  procedure Mixed_Cell_Tropisms
              ( report : in boolean;
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision;
                mv : out natural32 ) is

  -- DESCRIPTION :
  --   Given a system of n Laurent polynomials in n+1 variables,
  --   computes the tropisms, where the last variable is the parameter.
  --   The tropisms are computed via a regular mixed cell configuration,
  --   induced by the lifting defined by the last variable.
  --   As the lifting may differ, even if several supports would be the same,
  --   the type of mixed is assumed to be fully mixed.

  -- ON ENTRY :
  --   report   if intermediate output has to be written to screen;
  --   sup      a list of n supports in n+1 variables.

  -- ON RETURN :
  --   mcc      a mixed cell configuration induced by the lifting
  --            defined by the last exponent in each monomial of p;
  --   mv       the mixed volume of the cells in mcc.

    dim : constant integer32 := sup'last;
    mix : Standard_Integer_Vectors.Vector(1..dim);

  begin
    if report
     then put_line("The supports : "); put(sup);
    end if;
    mix := (mix'range => 1);
    if report then
      Integer_Create_Mixed_Cells(standard_output,dim,mix,false,sup,mcc);
      Integer_Volume_Computation(standard_output,dim,mix,true,sup,mcc,mv);
    else
      Integer_Create_Mixed_Cells(dim,mix,sup,mcc);
      Integer_Volume_Computation(dim,mix,true,sup,mcc,mv);
    end if;
  end Mixed_Cell_Tropisms;

  procedure Mixed_Cell_Tropisms
              ( file : in file_type;
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision;
                mv : out natural32 ) is

  -- DESCRIPTION :
  --   Given a system of n Laurent polynomials in n+1 variables,
  --   computes the tropisms, where the last variable is the parameter.
  --   The tropisms are computed via a regular mixed cell configuration,
  --   induced by the lifting defined by the last variable.
  --   As the lifting may differ, even if several supports would be the same,
  --   the type of mixed is assumed to be fully mixed.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   sup      a list of n supports in n+1 variables.

  -- ON RETURN :
  --   mcc      a mixed cell configuration induced by the lifting
  --            defined by the last exponent in each monomial of p;
  --   mv       the mixed volume of the cells in mcc.

    dim : constant integer32 := sup'last;
    mix : Standard_Integer_Vectors.Vector(1..dim);

  begin
    new_line(file);
    put_line(file,"THE SUPPORTS : ");
    put(file,sup);
    mix := (mix'range => 1);
    Integer_Create_Mixed_Cells(file,dim,mix,false,sup,mcc);
    Integer_Volume_Computation(file,dim,mix,true,sup,mcc,mv);
  end Mixed_Cell_Tropisms;

  procedure Initial_Coefficients
              ( file : in file_type;
                p : in Laur_Sys; mic : in Mixed_Cell;
                psub : out Laur_Sys; sols : out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the blackbox solver on the subsystem of p
  --   supported by the points that span the mixed cell mic.
  --
  -- ON ENTRY :
  --   file     to write output to file;
  --   p        system of n polynomials in n+1 variables;
  --   mic      a mixed cell with inner normal computed with a lifting
  --            equal to the exponent of the last variable in p.
  --
  -- ON RETURN :
  --   psub     a square subsystem of p with the last variable removed;
  --   sols     solutions to psub, as computed by the blackbox solver.

    q : Laur_Sys(p'range) := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(1.0);
    rc : natural32;

  begin
    put_line(file,"The initial form system with t :"); put(file,q);
    psub := Eval(q,one,idx);
    put_line(file,"The initial form system with t = 1 :"); put(file,psub);
    Solve(psub,false,rc,sols);
    put(file,"Computed ");
    put(file,Length_Of(sols),1); put_line(file," solutions.");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Initial_Coefficients;

  procedure Initial_Coefficients
              ( p : in Laur_Sys; mic : in Mixed_Cell;
                psub : out Laur_Sys; sols : out Solution_List;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Calls the blackbox solver on the subsystem of p
  --   supported by the points that span the mixed cell mic.
  --
  -- ON ENTRY :
  --   p        system of n polynomials in n+1 variables;
  --   mic      a mixed cell with inner normal computed with a lifting
  --            equal to the exponent of the last variable in p;
  --   report   if true, then output is written to screen.
  --
  -- ON RETURN :
  --   psub     a square subsystem of p with the last variable removed;
  --   sols     solutions to psub, as computed by the blackbox solver.

    q : Laur_Sys(p'range) := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(1.0);
    rc : natural32;

  begin
    if report then
      put_line("The initial form system with t :"); put(q);
    end if;
    psub := Eval(q,one,idx);
    if report then
      put_line("The initial form system with t = 1 :"); put(psub);
    end if;
    Solve(psub,false,rc,sols);
    if report then
      put("Computed "); put(Length_Of(sols),1); put_line(" solutions.");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Initial_Coefficients;

  function Pivot ( v : Standard_Integer_Vectors.Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns the first index k in v for which v(k) /= 0.

  begin
    for k in v'range loop
      if v(k) /= 0
       then return k;
      end if;
    end loop;
    return v'last+1;
  end Pivot;

  procedure Shift ( p : in out Poly; verbose : in boolean ) is

  -- DESCRIPTION :
  --   Multiplies the monomials in p so that all monomials have 
  --   nonnegative exponents.

    mindeg : constant Degrees := Minimal_Degrees(p);
    t : Term;

  begin
    if verbose then
      put("The minimal degrees : ");
      put(Standard_Integer_Vectors.Vector(mindeg.all)); new_line;
      put_line("The polynomial before the shift :");
      put(p); new_line;
    end if;
    t.cf := Create(1.0);
    t.dg := new Standard_Integer_Vectors.Vector(mindeg'range);
    for i in mindeg'range loop
      t.dg(i) := -mindeg(i);
    end loop;
    Mul(p,t);
    if verbose then
      put_line("The polynomial after the shift :");
      put(p); new_line;
    end if;
  end Shift;

  procedure Shift ( p : in out Laur_Sys; verbose : in boolean ) is

  -- DESCRIPTION :
  --   Multiplies the polynomials in p so that all monomials in p
  --   have nonnegative exponents.

  begin
    for i in p'range loop
      Shift(p(i),verbose);
    end loop;
  end Shift;

  procedure Transform_Coordinates
              ( file : in file_type;
                p : in Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Laur_Sys ) is

  -- DESCRIPTION :
  --   Applies the unimodular coordinate transformation defined by v
  --   to the system p, the result is in the system q.
  --   Output is written to file.

    i : constant integer32 := v'last; -- Pivot(v); -- stupid mistake!
    t : Transfo := Build_Transfo(v,i);

  begin
    q := Transform(t,p);
    Shift(q);
    put(file,"The transformation defined by ");
    put(file,v); put_line(file," :");
    Standard_Integer_Transformations_io.put(file,t);
    put_line(file,"The transformed system : "); put_line(file,q);
  end Transform_Coordinates;

  procedure Transform_Coordinates
              ( p : in Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Laur_Sys; report : in boolean ) is

  -- DESCRIPTION :
  --   Applies the unimodular coordinate transformation defined by v
  --   to the system p, the result is in the system q.

    i : constant integer32 := v'last; -- Pivot(v); -- stupid mistake!
    t : Transfo := Build_Transfo(v,i);

  begin
    q := Transform(t,p);
    Shift(q);
    if report then
      put("The transformation defined by "); put(v); put_line(" :");
      Standard_Integer_Transformations_io.put(t);
      put_line("The transformed system : "); put_line(q);
    end if;
  end Transform_Coordinates;

  function Initial_Residual
              ( p : in Laur_Sys;
                sol : in Standard_Complex_Vectors.Vector )
              return double_float is

  -- DESCRIPTION :
  --   Returns the residual of the solution extended with zero,
  --   evaluated at p.

    res : double_float := 0.0;
    ext : Standard_Complex_Vectors.Vector(sol'first..sol'last+1);
    eva : Standard_Complex_Vectors.Vector(p'range);

  begin
    ext(sol'range) := sol;
    ext(sol'last+1) := Create(0.0);
    eva := Eval(p,ext);
    -- res := Norm(eva);
    return res;
  end Initial_Residual;

  function Initial_Residuals
              ( file : in file_type; p : in Laur_Sys;
                sols : in Solution_List ) return double_float is

  -- DESCRIPTION :
  --   Returns the accumulated sum of all residuals of the solutions,
  --   extended with zero for the last coordinate.
  --   The residuals are written to file for each solution.

    res : double_float := 0.0;
    eva : double_float;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    for k in 1..Length_Of(tmp) loop
      ls := Head_Of(tmp);
      eva := Initial_Residual(p,ls.v);
      put(file,"At solution "); put(file,k,1);
      put(file," : "); put(file,eva); new_line(file);
      res := res + eva;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Initial_Residuals;

  function Initial_Residuals
              ( p : in Laur_Sys;
                sols : in Solution_List;
                report : in boolean ) return double_float is

  -- DESCRIPTION :
  --   Returns the accumulated sum of all residuals of the solutions,
  --   extended with zero for the last coordinate.
  --   If report, then the residuals are written to screen for each solution,
  --   otherwise, this function remains silent.

    res : double_float := 0.0;
    eva : double_float;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    for k in 1..Length_Of(tmp) loop
      ls := Head_Of(tmp);
      eva := Initial_Residual(p,ls.v);
      if report then
        put("At solution "); put(k,1);
        put(" : "); put(eva); new_line;
      end if;
      res := res + eva;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Initial_Residuals;

  procedure Initials
              ( file : in file_type;
                p : in Laur_Sys;
                mcc : in Mixed_Subdivision;
                mv : in natural32 ) is

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the cells in mcc.

  -- ON ENTRY :
  --   file     to write output;
  --   p        a Laurent system of n equations in n+1 variables,
  --            where the last variable is the lifting;
  --   mcc      regular mixed cell configuration defined by the
  --            lower hull of the supports of p;
  --   mv       the mixed volume gives a generically sharp upper
  --            bound on the number of initial coefficients.

    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    tvp : Laur_Sys(p'range);

  begin
    for k in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        psub : Laur_Sys(p'range);
        sols : Solution_List;
        res : double_float;
      begin
        put(file,"Tropism "); put(file,k,1); put(file," is ");
        put(file,mic.nor); new_line(file);
        Initial_Coefficients(file,p,mic,psub,sols);
        Transform_Coordinates(file,p,mic.nor.all,tvp);
        res := Initial_Residuals(file,tvp,sols);
        put(file,"The residual :"); put(file,res,3); new_line(file);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Initials;

  procedure Initials
              ( p : in Laur_Sys;
                mcc : in Mixed_Subdivision;
                mv : in natural32; report : in boolean ) is

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the cells in mcc.

  -- ON ENTRY :
  --   p        a Laurent system of n equations in n+1 variables,
  --            where the last variable is the lifting;
  --   mcc      regular mixed cell configuration defined by the
  --            lower hull of the supports of p;
  --   mv       the mixed volume gives a generically sharp upper
  --            bound on the number of initial coefficients;
  --   report   if true, then output is written to screen.

    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    tvp : Laur_Sys(p'range);

  begin
    for k in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        psub : Laur_Sys(p'range);
        sols : Solution_List;
        res : double_float;
      begin
        if report then
          put("Tropism "); put(k,1); put(" is ");
          put(mic.nor); new_line;
        end if;
        Initial_Coefficients(p,mic,psub,sols,report);
        Transform_Coordinates(p,mic.nor.all,tvp,report);
        res := Initial_Residuals(tvp,sols,report);
        if report then
          put("The residual :"); put(res,3); new_line;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Initials;

  procedure Test ( file : in file_type; p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and initial terms.
  --   The output is written to file.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;

  begin
    Mixed_Cell_Tropisms(file,sup,mcc,mv);
    new_line(file);
    put(file,"The number of pretropisms : ");
    put(file,Length_Of(mcc),1); new_line(file);
    put(file,"The number of series : ");
    put(file,mv,1); new_line(file);
    Initials(file,p,mcc,mv);
  end Test;

  procedure Test ( p : in Laur_Sys; report : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and initial terms.
  --   The output is written to screen if report.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;

  begin
    Mixed_Cell_Tropisms(report,sup,mcc,mv);
    if report then
      new_line;
      put("The number of pretropisms : ");
      put(Length_Of(mcc),1); new_line;
      put("The number of series : "); put(mv,1); new_line;
    end if;
    Initials(p,mcc,mv,report);
  end Test;
  
  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.

    lp : Link_to_Laur_Sys;
    nq,nv : integer32;
    ans : character;
    report : boolean;
    file : file_type;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
    if nv /= nq+1 then
      put(nv,1); put(" /= "); put(nq,1); put(" + 1");
    else
      new_line;
      put("Do you want intermediate output to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      report := (ans = 'y');
      if report then
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(file);
        put(file,natural32(nq),natural32(nv),lp.all);
        Test(file,lp.all);
      else
        new_line;
        put("Do you want intermediate output to screen ? (y/n) ");
        Ask_Yes_or_No(ans);
        report := (ans = 'y');
        Test(lp.all,report);
      end if;
    end if;
  end Main;

begin
  Main;
end ts_puiseux;
