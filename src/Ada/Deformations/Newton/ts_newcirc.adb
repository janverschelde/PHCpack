with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Standard_Coefficient_Circuits;      use Standard_Coefficient_Circuits;
with Standard_Circuit_Makers;            use Standard_Circuit_Makers;
with Standard_Newton_Circuits;           use Standard_Newton_Circuits;

procedure ts_newcirc is

-- DESCRIPTION :
--   Test the development of Newton's method on coefficient circuits.

  procedure LU_Newton_Steps
              ( s : in Link_to_System;
                v : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton steps on the system s starting at the vector v.
 
  -- REQUIRED : The system is square.

    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    info : integer32;
    vxr : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    res,rco,err : double_float;
    ans : character;
    condition : boolean;

  begin
    put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
    condition := (ans = 'y');
    put_line("The vector v :"); put_line(v);
    loop
      if condition
       then LU_Newton_Step(s,v,xr,xi,ipvt,res,rco,err,true);
       else LU_Newton_Step(s,v,xr,xi,ipvt,info,res,err,true);
      end if;
      put_line("The vector v :"); put_line(v);
      put("Another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end LU_Newton_Steps;

  procedure Interactive_Run
              ( s : in Link_to_System; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Interactive application of Newton's method on the system of
  --   coefficient circuits s, starting at the solutions in sols.
  --   Pauses frequently, prompting user to continue.

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    ans : character;
    cnt : integer32 := 0;

  begin
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      put("Running Newton's method on solution ");
      put(cnt,1); put_line(" ...");
      LU_Newton_Steps(s,ls.v);
      put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      ptr := Tail_Of(ptr);
    end loop;
  end Interactive_Run;

  procedure Show_Parameters
              ( maxit : in natural32;
                tolres,tolerr : in double_float;
                condition : in boolean ) is
  begin
    put_line("Parameter Settings :");
    put("  1. maximum number of iterations : "); put(maxit,1); new_line;
    put("  2. tolerance on residual        :"); put(tolres,3); new_line;
    put("  3. tolerance on forward error   :"); put(tolerr,3); new_line;
    put("  4. condition number wanted      : ");
    if condition
     then put_line("yes");
     else put_line("no");
    end if;
  end Show_Parameters;

  procedure Set_Parameters
              ( maxit : out natural32;
                tolres,tolerr : out double_float;
                condition : out boolean ) is

  -- DESCRIPTION :
  --   Sets the parameters to run several steps with Newton's method.

    ans : character;

  begin
    maxit := 3; tolres := 1.0E-12; tolerr := 1.0E-12; condition := true; 
    loop
      Show_Parameters(maxit,tolres,tolerr,condition);
      put("Type 1, 2, 3, or 4 to set parameter, or 0 to exit : ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      case ans is
        when '1' => put("-> maximum number of iterations : "); get(maxit);
        when '2' => put("-> tolerance on residual : "); get(tolres);
        when '3' => put("-> tolerance on forward error : "); get(tolerr);
        when '4' => put("-> condition number wanted ? (y/n) ");
                    Ask_Yes_or_No(ans); condition := (ans = 'y');
        when others => null;
      end case;
    end loop;
  end Set_Parameters;

  procedure Monitored_Run
              ( s : in Link_to_System; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Runs several steps of Newton's method on the system s,
  --   starting at the solutions in sols.
  --   For each solution writes one line to screen.

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    info : integer32;
    vxr : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    res,rco,err,tolres,tolerr : double_float;
    numit,maxit : natural32 := 0;
    fail,condition : boolean;

  begin
    Set_Parameters(maxit,tolres,tolerr,condition);
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1; put(cnt,1); put(" : ");
      if condition then
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        res,rco,err,numit,fail,false);
        put(cnt,1); put(" : ");
        put("  err :"); put(err,3); put("  rco :"); put(rco,3);
        put("  res :"); put(res,3); put("  #steps : "); put(numit);
        if fail
         then put_line("  failure");
         else put_line("  success");
        end if;
      else
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        info,res,err,numit,fail,false);
        put("  err :"); put(err,3);
        put("  res :"); put(res,3); put("  #steps : "); put(numit);
        if fail
         then put_line("  failure");
         else put_line("  success");
        end if;
      end if;
      ptr := Tail_Of(ptr);
    end loop;
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end Monitored_Run;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions,
  --   makes a coefficient circuit and then runs Newton's method.

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    nbq,len,dim : integer32 := 0;
    s : Link_to_System;
    ans : character;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(p,sols,"THE SOLUTIONS");
    nbq := p'last;
    len := integer32(Length_Of(sols));
    if len = 0 then
      put_line("No solutions found on file.");
    else
      dim := Head_Of(sols).n;
      new_line;
      put("Read system of "); put(nbq,1); put(" polynomials and ");
      put(len,1); put(" solutions in dimension "); put(dim,1); put_line(".");
      s := Make_Coefficient_System(p);
      new_line;
      put("Interactive run ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then Interactive_Run(s,sols);
       else Monitored_Run(s,sols);
      end if;
    end if;
  end Main;

begin
  Main;
end ts_newcirc;
