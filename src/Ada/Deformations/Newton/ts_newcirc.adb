with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with Standard_Solution_Diagnostics;
with Standard_Condition_Tables;
with Standard_Coefficient_Circuits;      use Standard_Coefficient_Circuits;
with Standard_Circuit_Makers;            use Standard_Circuit_Makers;
with Standard_Newton_Circuits;           use Standard_Newton_Circuits;
with Standard_Solutions_Heap;

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
                tolres,tolerr,tolsing : in double_float;
                condition : in boolean ) is

  -- DESCRIPTION :
  --   Displays the values of the parameters.

  -- ON ENTRY :
  --   maxit        maximum number of iterations;
  --   tolres       tolerance on the residual;
  --   tolerr       tolerance on the forward error;
  --   tolsing      tolerance on a singularity;
  --   condition    true if condition number is wanted,
  --                false otherwise.

  begin
    put_line("Parameter Settings :");
    put("  1. maximum number of iterations : "); put(maxit,1); new_line;
    put("  2. tolerance on residual        :"); put(tolres,3); new_line;
    put("  3. tolerance on forward error   :"); put(tolerr,3); new_line;
    put("  4. tolerance on singularity     :"); put(tolsing,3); new_line;
    put("  5. condition number wanted      : ");
    if condition
     then put_line("yes");
     else put_line("no");
    end if;
  end Show_Parameters;

  procedure Set_Parameters
              ( maxit : out natural32;
                tolres,tolerr,tolsing : out double_float;
                condition : out boolean ) is

  -- DESCRIPTION :
  --   Sets the parameters to run several steps with Newton's method.

  -- ON RETURN :
  --   maxit        maximum number of iterations;
  --   tolres       tolerance on the residual;
  --   tolerr       tolerance on the forward error;
  --   tolsing      tolerance on a singularity;
  --   condition    true if condition number is wanted,
  --                false otherwise.

    ans : character;

  begin
    maxit := 4; condition := true;
    tolres := 1.0E-8; tolerr := 1.0E-8; tolsing := 1.0E-8;
    loop
      Show_Parameters(maxit,tolres,tolerr,tolsing,condition);
      put("Type 1, 2, 3, 4, or 5 to set parameter, or 0 to exit : ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      case ans is
        when '1' => put("-> maximum number of iterations : "); get(maxit);
        when '2' => put("-> tolerance on residual :"); get(tolres);
        when '3' => put("-> tolerance on forward error :"); get(tolerr);
        when '4' => put("-> tolerance on singularity :"); get(tolsing);
        when '5' => put("-> condition number wanted ? (y/n) ");
                    Ask_Yes_or_No(ans); condition := (ans = 'y');
        when others => null;
      end case;
    end loop;
  end Set_Parameters;

  procedure Monitor_Report
              ( idx : in integer32; fail,isreal : in boolean;
                err,rco,res,wgt,tolsing : in double_float ) is

  -- DESCRIPTION :
  --   Writes one line to screen with a report on a solution,
  --   to monitor the progress on the verification.

  -- ON ENTRY :
  --   idx      index number of the current solution;
  --   fail     true if Newton's method failed, false otherwise;
  --   isreal   true if real, false otherwise (only if not fail);
  --   err      forward error;
  --   rco      estimate for the inverse condition number;
  --   res      residual;
  --   wgt      weight of the coordinates;
  --   tolsing  tolerance to decide whether a solution is singular.

  begin
    put(idx,1); put(" : ");
    if fail then
      put_line("no solution");
    else
      put("err :"); put(err,2);
      put("  rco :"); put(rco,2);
      put("  res :"); put(res,2);
      put("  wgt :"); put(wgt,2);
      if isreal
       then put(" real");
       else put(" complex");
      end if;
      if rco < tolsing
       then put_line(" singular");
       else put_line(" regular");
      end if;
    end if;
  end Monitor_Report;

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
    startres,res,rco,err,tolres,tolerr,tolsing : double_float;
    numit,maxit : natural32 := 0;
    fail,condition : boolean;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 

  begin
    new_line;
    Set_Parameters(maxit,tolres,tolerr,tolsing,condition);
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1; put(cnt,1); put(" : ");
      if condition then
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        startres,res,rco,err,numit,fail,false);
        put("  err :"); put(err,3); put("  rco :"); put(rco,3);
        put("  res :"); put(res,3); put("  #steps : "); put(numit);
        if fail
         then put_line("  failure");
         else put_line("  success");
        end if;
        Standard_Condition_Tables.Update_Corrector(t_err,err);
        Standard_Condition_Tables.Update_Condition(t_rco,rco);
        Standard_Condition_Tables.Update_Residuals(t_res,res);
      else
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        info,startres,res,err,numit,fail,false);
        put("  err :"); put(err,3);
        put("  res :"); put(res,3); put("  #steps : "); put(numit);
        if fail
         then put_line("  failure");
         else put_line("  success");
        end if;
        Standard_Condition_Tables.Update_Corrector(t_err,err);
        Standard_Condition_Tables.Update_Residuals(t_res,res);
      end if;
      ptr := Tail_Of(ptr);
    end loop;
    Standard_Condition_Tables.Write_Tables(standard_output,t_err,t_res,t_rco);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end Monitored_Run;

  procedure Monitored_Run
              ( file : in file_type;
                s : in Link_to_System; sols : in Solution_List;
               -- len : in integer32;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs several steps of Newton's method on the system s,
  --   starting at the solutions in sols.
  --   Writes the output to file.
  --   If verbose, then one line is written for every solution.

    use Standard_Solutions_Heap;

    timer : Timing_Widget;
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
    startres,res,rco,err,tolres,tolerr,tolsing : double_float;
    cntfail,cntreal,cntcmplx,cntregu,cntsing,cntclus : natural32 := 0;
    numit,maxit : natural32 := 0;
    fail,condition,isreal : boolean;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 
    nbr : constant integer32 := 2*s.dim;
    wv1 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Random_Weight_Vector(nbr);
    wv2 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Random_Weight_Vector(nbr);
    val1,val2 : double_float;
    weights : Heap; -- Heap(len);

  begin
    new_line;
    Set_Parameters(maxit,tolres,tolerr,tolsing,condition);
    new_line;
    weights.bottom := -1; -- make sure heap is declared as empty
    new_line;
    put_line("See the output file for results ...");
    new_line;
    tstart(timer);
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      put(file,"Solution "); put(file,cnt,1);
      put(file," :    start residual :");
      if condition then
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        startres,res,rco,err,numit,fail,false);
        put(file,startres,3);
        put(file,"  #iterations : "); put(file,numit,1);
        if fail
         then put_line(file,"  failure");
         else put_line(file,"  success");
        end if;
        Standard_Complex_Solutions_io.put_vector(file,ls.v);
        put(file,"== err :"); put(file,err,3);
        put(file," = rco :"); put(file,rco,3);
        put(file," = res :"); put(file,res,3);
        if fail then
          put_line(file," == no solution"); cntfail := cntfail + 1;
        else
          isreal := Standard_Solution_Diagnostics.Is_Real(ls.all,tolsing);
          val1 := Weight(ls.v,wv1);
          val2 := Weight(ls.v,wv2);
          Push(weights,val1,val2,cnt,ls);
          if isreal
           then put(file," == real");    cntreal := cntreal + 1;
           else put(file," == complex"); cntcmplx := cntcmplx + 1;
          end if;
          if rco < tolsing
           then put_line(file," singular"); cntsing := cntsing + 1;
           else put_line(file," regular");  cntregu := cntregu + 1;
          end if;
        end if;
        if verbose
         then Monitor_Report(cnt,fail,isreal,err,rco,res,val1,tolsing);
        end if;
        Standard_Condition_Tables.Update_Corrector(t_err,err);
        Standard_Condition_Tables.Update_Condition(t_rco,rco);
        Standard_Condition_Tables.Update_Residuals(t_res,res);
        ls.rco := rco;
      else
        LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                        info,startres,res,err,numit,fail,false);
        put(file,startres,3);
        put(file,"  #iterations : "); put(file,numit,1);
        if fail
         then put_line(file,"  failure");
         else put_line(file,"  success");
        end if;
        Standard_Complex_Solutions_io.put_vector(file,ls.v);
        put(file,"== err :"); put(file,err,3);
        put(file," = res :"); put(file,res,3);
        if fail then
          put_line(file," == no solution"); cntfail := cntfail + 1;
        else
          isreal := Standard_Solution_Diagnostics.Is_Real(ls.all,tolsing);
          val1 := Weight(ls.v,wv1);
          val2 := Weight(ls.v,wv2);
          Push(weights,val1,val2,cnt,ls);
          if isreal
           then put(file," == real solution");    cntreal := cntreal + 1;
           else put(file," == complex solution"); cntcmplx := cntcmplx + 1;
          end if;
        end if;
        Standard_Condition_Tables.Update_Corrector(t_err,err);
        Standard_Condition_Tables.Update_Residuals(t_res,res);
      end if;
      ls.err := err; ls.res := res;
      ptr := Tail_Of(ptr);
    end loop;
    tstop(timer);
    Standard_Complex_Solutions_io.put_bar(file);
    Count_Clusters(weights,tolsing,cntclus);
    if condition then
      put(file,"number of regular solutions   : ");
      put(file,cntregu,1); new_line(file);
      put(file,"number of singular solutions  : ");
      put(file,cntsing,1); new_line(file);
    end if;
    put(file,"number of real solutions      : ");
    put(file,cntreal,1); new_line(file);
    put(file,"number of complex solutions   : ");
    put(file,cntcmplx,1); new_line(file);
    put(file,"number of clustered solutions : ");
    put(file,cntclus,1); new_line(file);
    put(file,"number of failures            : ");
    put(file,cntfail,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
    new_line(file);
    print_times(file,timer,"Newton with condition table report");
    Clear(weights);
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
    file : file_type;

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
      if ans = 'y' then
        Interactive_Run(s,sols);
      else
        put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'n' then
          Monitored_Run(s,sols);
        else
          new_line;
          put_line("Reading the name of the output file ...");
          Read_Name_and_Create_File(file);
          Monitored_Run(file,s,sols); --,len);
          Close(file);
        end if;
      end if;
    end if;
  end Main;

begin
  Main;
end ts_newcirc;
