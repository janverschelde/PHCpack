with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_System_and_Solutions_io;
with Standard_Complex_Solutions_io;
with Standard_Solution_Diagnostics;
with Standard_Condition_Tables;
with Standard_Circuit_Makers;
with Standard_Newton_Circuits;
with Standard_Inlined_Newton_Circuits;
with Standard_Solutions_Heap;

package body Standard_Refiner_Circuits is

  procedure Show_Parameters
              ( maxit : in natural32;
                tolres,tolerr,tolsing : in double_float ) is
  begin
    put_line("Parameter Settings :");
    put("  1. maximum number of iterations : "); put(maxit,1); new_line;
    put("  2. tolerance on residual        :"); put(tolres,3); new_line;
    put("  3. tolerance on forward error   :"); put(tolerr,3); new_line;
    put("  4. tolerance on singularity     :"); put(tolsing,3); new_line;
  end Show_Parameters;

  procedure Set_Parameters
              ( maxit : out natural32;
                tolres,tolerr,tolsing : out double_float ) is

    ans : character;

  begin
    maxit := 4; tolres := 1.0E-8; tolerr := 1.0E-8; tolsing := 1.0E-8;
    loop
      Show_Parameters(maxit,tolres,tolerr,tolsing);
      put("Type 1, 2, 3, or 4 to set a parameter, or 0 to exit : ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      case ans is
        when '1' => put("-> maximum number of iterations : "); get(maxit);
        when '2' => put("-> tolerance on residual :"); get(tolres);
        when '3' => put("-> tolerance on forward error :"); get(tolerr);
        when '4' => put("-> tolerance on singularity :"); get(tolsing);
        when others => null;
      end case;
    end loop;
  end Set_Parameters;

  procedure Monitor_Report
              ( idx : in integer32; fail,isreal : in boolean;
                err,rco,res,wgt,tolsing : in double_float ) is
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

  procedure Run ( s : in Link_to_System; sols : in Solution_List;
                  maxit : in natural32;
                  tolres,tolerr,tolsing : in double_float;
                  cntfail,cntreal,cntcmplx : out natural32;
                  cntregu,cntsing,cntclus : out natural32;
                  t_err,t_rco,t_res : out Standard_Natural_Vectors.Vector;
                  verbose : in boolean; vrb : in integer32 := 0 ) is

    use Standard_Newton_Circuits;

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    vxr : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    startres,res,rco,err : double_float;
    numit : natural32 := 0;
    fail,isreal : boolean;
    nbr : constant integer32 := 2*s.dim;
    wv1 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Standard_Solutions_Heap.Random_Weight_Vector(nbr);
    wv2 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Standard_Solutions_Heap.Random_Weight_Vector(nbr);
    val1,val2 : double_float;
    weights : Standard_Solutions_Heap.Heap;

  begin
    cntfail := 0; cntreal := 0; cntcmplx := 0;
    cntregu := 0; cntsing := 0; cntclus := 0;
    t_err := Standard_Condition_Tables.Create(natural32(t_err'last));
    t_rco := Standard_Condition_Tables.Create(natural32(t_err'last));
    t_res := Standard_Condition_Tables.Create(natural32(t_err'last));
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Run 1");
    end if;
    if verbose
     then new_line;
    end if;
    weights.bottom := -1;
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      if verbose then
        put("Solution "); put(cnt,1); put(" :    start residual :");
      end if;
      LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                      startres,res,rco,err,numit,fail);
      if verbose then
        put(startres,3);
        put("  #iterations : "); put(numit,1);
        if fail
         then put_line("  failure");
         else put_line("  success");
        end if;
        put("t :"); put(ls.t); new_line;
        put("m : "); put(ls.m,1); new_line;
        put_line("the solution for t :");
        Standard_Complex_Solutions_io.put_vector(ls.v);
        put("err :"); put(err,3);
        put(" = rco :"); put(rco,3);
        put(" = res :"); put(res,3);
      end if;
      if fail then
        if verbose
         then put_line(" == no solution");
        end if;
        cntfail := cntfail + 1;
      else
        isreal := Standard_Solution_Diagnostics.Is_Real(ls.all,tolsing);
        val1 := Standard_Solutions_Heap.Weight(ls.v,wv1);
        val2 := Standard_Solutions_Heap.Weight(ls.v,wv2);
        Standard_Solutions_Heap.Push(weights,val1,val2,cnt,ls);
        if isreal then
          if verbose
           then put(" == real");
          end if;
          cntreal := cntreal + 1;
        else
          if verbose
           then put(" == complex");
          end if;
          cntcmplx := cntcmplx + 1;
        end if;
        if rco < tolsing then
          if verbose
           then put_line(" singular");
          end if;
          cntsing := cntsing + 1;
        else
          if verbose
           then put_line(" regular");
          end if;
          cntregu := cntregu + 1;
        end if;
      end if;
      Standard_Condition_Tables.Update_Corrector(t_err,err);
      Standard_Condition_Tables.Update_Condition(t_rco,rco);
      Standard_Condition_Tables.Update_Residuals(t_res,res);
      ptr := Tail_Of(ptr);
    end loop;
    if verbose
     then put_line("computing clusters ...");
    end if;
    Standard_Solutions_Heap.Count_Clusters(weights,tolsing,cntclus,verbose);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
    Standard_Solutions_Heap.Clear(weights);
  end Run;

  procedure Inlined_Run
              ( s : in Link_to_System; sols : in Solution_List;
                maxit : in natural32;
                tolres,tolerr,tolsing : in double_float;
                cntfail,cntreal,cntcmplx : out natural32;
                cntregu,cntsing,cntclus : out natural32;
                t_err,t_rco,t_res : out Standard_Natural_Vectors.Vector;
                verbose : in boolean; vrb : in integer32 := 0 ) is

    use Standard_Inlined_Newton_Circuits;

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    vxr : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    startres,res,rco,err : double_float;
    numit : natural32 := 0;
    fail,isreal : boolean;
    nbr : constant integer32 := 2*s.dim;
    wv1 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Standard_Solutions_Heap.Random_Weight_Vector(nbr);
    wv2 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Standard_Solutions_Heap.Random_Weight_Vector(nbr);
    val1,val2 : double_float;
    weights : Standard_Solutions_Heap.Heap;

  begin
    cntfail := 0; cntreal := 0; cntcmplx := 0;
    cntregu := 0; cntsing := 0; cntclus := 0;
    t_err := Standard_Condition_Tables.Create(15); 
    t_rco := Standard_Condition_Tables.Create(15); 
    t_res := Standard_Condition_Tables.Create(15); 
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Inlined_Run 1");
    end if;
    weights.bottom := -1;
    if verbose
     then new_line;
    end if;
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      if verbose then
        put("Solution "); put(cnt,1); put(" :    start residual :");
      end if;
      LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                      startres,res,rco,err,numit,fail);
      if verbose then
        put(startres,3);
        put("  #iterations : "); put(numit,1);
        if fail
         then put_line("  failure");
         else put_line("  success");
        end if;
        put("t :"); put(ls.t); new_line;
        put("m : "); put(ls.m,1); new_line;
        put_line("the solution for t :");
        Standard_Complex_Solutions_io.put_vector(ls.v);
        put("err :"); put(err,3);
        put(" = rco :"); put(rco,3);
        put(" = res :"); put(res,3);
      end if;
      if fail then
        if verbose
         then put_line(" == no solution");
        end if;
        cntfail := cntfail + 1;
      else
        isreal := Standard_Solution_Diagnostics.Is_Real(ls.all,tolsing);
        val1 := Standard_Solutions_Heap.Weight(ls.v,wv1);
        val2 := Standard_Solutions_Heap.Weight(ls.v,wv2);
        Standard_Solutions_Heap.Push(weights,val1,val2,cnt,ls);
        if isreal then
          if verbose
           then put(" == real");
          end if;
          cntreal := cntreal + 1;
        else
          if verbose
           then put(" == complex");
          end if;
          cntcmplx := cntcmplx + 1;
        end if;
        if rco < tolsing then
          if verbose
           then put_line(" singular");
          end if;
          cntsing := cntsing + 1;
        else
          if verbose
           then put_line(" regular");
          end if;
          cntregu := cntregu + 1;
        end if;
      end if;
      Standard_Condition_Tables.Update_Corrector(t_err,err);
      Standard_Condition_Tables.Update_Condition(t_rco,rco);
      Standard_Condition_Tables.Update_Residuals(t_res,res);
      ptr := Tail_Of(ptr);
    end loop;
    if verbose
     then put_line("computing clusters ...");
    end if;
    Standard_Solutions_Heap.Count_Clusters(weights,tolsing,cntclus,verbose);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
    Standard_Solutions_Heap.Clear(weights);
  end Inlined_Run;

  procedure Run ( file : in file_type;
                  s : in Link_to_System; sols : in Solution_List;
                  maxit : in natural32;
                  tolres,tolerr,tolsing : in double_float;
                  cntfail,cntreal,cntcmplx : out natural32;
                  cntregu,cntsing,cntclus : out natural32;
                  t_err,t_rco,t_res : out Standard_Natural_Vectors.Vector;
                  verbose : in boolean; vrb : in integer32 := 0 ) is

    use Standard_Newton_Circuits;

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    vxr : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    startres,res,rco,err : double_float;
    numit : natural32 := 0;
    fail,isreal : boolean;
    nbr : constant integer32 := 2*s.dim;
    wv1 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Standard_Solutions_Heap.Random_Weight_Vector(nbr);
    wv2 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Standard_Solutions_Heap.Random_Weight_Vector(nbr);
    val1,val2 : double_float;
    weights : Standard_Solutions_Heap.Heap; -- Heap(len);

  begin
    cntfail := 0; cntreal := 0; cntcmplx := 0;
    cntregu := 0; cntsing := 0; cntclus := 0;
    t_err := Standard_Condition_Tables.Create(15); 
    t_rco := Standard_Condition_Tables.Create(15); 
    t_res := Standard_Condition_Tables.Create(15); 
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Run 2");
    end if;
    weights.bottom := -1; -- make sure heap is declared as empty
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),1); put(file," "); put(file,s.dim,1);
    new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      put(file,"Solution "); put(file,cnt,1);
      put(file," :    start residual :");
      LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                      startres,res,rco,err,numit,fail);
      put(file,startres,3);
      put(file,"  #iterations : "); put(file,numit,1);
      if fail
       then put_line(file,"  failure");
       else put_line(file,"  success");
      end if;
      put(file,"t :"); put(file,ls.t); new_line(file);
      put(file,"m : "); put(file,ls.m,1); new_line(file);
      put_line(file,"the solution for t :");
      Standard_Complex_Solutions_io.put_vector(file,ls.v);
      put(file,"== err :"); put(file,err,3);
      put(file," = rco :"); put(file,rco,3);
      put(file," = res :"); put(file,res,3);
      if fail then
        put_line(file," == no solution"); cntfail := cntfail + 1;
      else
        isreal := Standard_Solution_Diagnostics.Is_Real(ls.all,tolsing);
        val1 := Standard_Solutions_Heap.Weight(ls.v,wv1);
        val2 := Standard_Solutions_Heap.Weight(ls.v,wv2);
        Standard_Solutions_Heap.Push(weights,val1,val2,cnt,ls);
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
      ls.err := err; ls.res := res;
      ptr := Tail_Of(ptr);
    end loop;
    Standard_Complex_Solutions_io.put_bar(file);
    Standard_Solutions_Heap.Count_Clusters(weights,tolsing,cntclus,verbose);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end Run;

  procedure Inlined_Run
              ( file : in file_type;
                s : in Link_to_System; sols : in Solution_List;
                maxit : in natural32;
                tolres,tolerr,tolsing : in double_float;
                cntfail,cntreal,cntcmplx : out natural32;
                cntregu,cntsing,cntclus : out natural32;
                t_err,t_rco,t_res : out Standard_Natural_Vectors.Vector;
                verbose : in boolean; vrb : in integer32 := 0 ) is

    use Standard_Inlined_Newton_Circuits;

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    vxr : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(1..s.dim)
        := (1..s.dim => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    startres,res,rco,err : double_float;
    numit : natural32 := 0;
    fail,isreal : boolean;
    nbr : constant integer32 := 2*s.dim;
    wv1 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Standard_Solutions_Heap.Random_Weight_Vector(nbr);
    wv2 : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Standard_Solutions_Heap.Random_Weight_Vector(nbr);
    val1,val2 : double_float;
    weights : Standard_Solutions_Heap.Heap; -- Heap(len);

  begin
    cntfail := 0; cntreal := 0; cntcmplx := 0;
    cntregu := 0; cntsing := 0; cntclus := 0;
    t_err := Standard_Condition_Tables.Create(15); 
    t_rco := Standard_Condition_Tables.Create(15); 
    t_res := Standard_Condition_Tables.Create(15); 
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Inlined_Run 2");
    end if;
    weights.bottom := -1; -- make sure heap is declared as empty
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),1); put(file," "); put(file,s.dim,1);
    new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      put(file,"Solution "); put(file,cnt,1);
      put(file," :    start residual :");
      LU_Newton_Steps(s,ls.v,xr,xi,maxit,tolres,tolerr,ipvt,
                      startres,res,rco,err,numit,fail);
      put(file,startres,3);
      put(file,"  #iterations : "); put(file,numit,1);
      if fail
       then put_line(file,"  failure");
       else put_line(file,"  success");
      end if;
      put(file,"t :"); put(file,ls.t); new_line(file);
      put(file,"m : "); put(file,ls.m,1); new_line(file);
      put_line(file,"the solution for t :");
      Standard_Complex_Solutions_io.put_vector(file,ls.v);
      put(file,"== err :"); put(file,err,3);
      put(file," = rco :"); put(file,rco,3);
      put(file," = res :"); put(file,res,3);
      if fail then
        put_line(file," == no solution"); cntfail := cntfail + 1;
      else
        isreal := Standard_Solution_Diagnostics.Is_Real(ls.all,tolsing);
        val1 := Standard_Solutions_Heap.Weight(ls.v,wv1);
        val2 := Standard_Solutions_Heap.Weight(ls.v,wv2);
        Standard_Solutions_Heap.Push(weights,val1,val2,cnt,ls);
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
      ls.err := err; ls.res := res;
      ptr := Tail_Of(ptr);
    end loop;
    Standard_Complex_Solutions_io.put_bar(file);
    Standard_Solutions_Heap.Count_Clusters(weights,tolsing,cntclus,verbose);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end Inlined_Run;

  procedure Run ( s : in Link_to_System; sols : in Solution_List;
                  vrb : in integer32 := 0 ) is

    tolres,tolerr,tolsing : double_float;
    cntfail,cntreal,cntcmplx,cntregu,cntsing,cntclus : natural32 := 0;
    maxit : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15);
    ans : character;
    verbose : boolean;

  begin
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Run 3");
    end if;
    new_line;
    Set_Parameters(maxit,tolres,tolerr,tolsing);
    new_line;
    put("See all solutions ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    Run(s,sols,maxit,tolres,tolerr,tolsing,cntfail,cntreal,cntcmplx,
        cntregu,cntsing,cntclus,t_err,t_rco,t_res,verbose,vrb);
    put("number of regular solutions   : "); put(cntregu,1); new_line;
    put("number of singular solutions  : "); put(cntsing,1); new_line;
    put("number of real solutions      : "); put(cntreal,1); new_line;
    put("number of complex solutions   : "); put(cntcmplx,1); new_line;
    put("number of clustered solutions : "); put(cntclus,1); new_line;
    put("number of failures            : "); put(cntfail,1); new_line;
    Standard_Condition_Tables.Write_Tables(standard_output,t_err,t_res,t_rco);
  end Run;

  procedure Inlined_Run
              ( s : in Link_to_System; sols : in Solution_List;
                vrb : in integer32 := 0 ) is

    tolres,tolerr,tolsing : double_float;
    cntfail,cntreal,cntcmplx,cntregu,cntsing,cntclus : natural32 := 0;
    maxit : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15);
    ans : character;
    verbose : boolean;

  begin
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Inlined_Run 3");
    end if;
    new_line;
    Set_Parameters(maxit,tolres,tolerr,tolsing);
    new_line;
    put("See all solutions ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    Inlined_Run(s,sols,maxit,tolres,tolerr,tolsing,cntfail,cntreal,cntcmplx,
                cntregu,cntsing,cntclus,t_err,t_rco,t_res,verbose,vrb);
    put("number of regular solutions   : "); put(cntregu,1); new_line;
    put("number of singular solutions  : "); put(cntsing,1); new_line;
    put("number of real solutions      : "); put(cntreal,1); new_line;
    put("number of complex solutions   : "); put(cntcmplx,1); new_line;
    put("number of clustered solutions : "); put(cntclus,1); new_line;
    put("number of failures            : "); put(cntfail,1); new_line;
    Standard_Condition_Tables.Write_Tables(standard_output,t_err,t_res,t_rco);
  end Inlined_Run;

  procedure Run ( file : in file_type;
                  s : in Link_to_System; sols : in Solution_List;
                  vrb : in integer32 := 0 ) is

    timer : Timing_Widget;
    tolres,tolerr,tolsing : double_float;
    cntfail,cntreal,cntcmplx,cntregu,cntsing,cntclus : natural32 := 0;
    maxit : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15);
    ans : character;
    verbose : boolean;

  begin
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Run 4");
    end if;
    new_line;
    Set_Parameters(maxit,tolres,tolerr,tolsing);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    tstart(timer);
    Run(file,s,sols,maxit,tolres,tolerr,tolsing,cntfail,cntreal,cntcmplx,
        cntregu,cntsing,cntclus,t_err,t_rco,t_res,verbose,vrb);
    tstop(timer);
    put(file,"number of regular solutions   : ");
    put(file,cntregu,1); new_line(file);
    put(file,"number of singular solutions  : ");
    put(file,cntsing,1); new_line(file);
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
    new_line(file);
    print_times(file,timer,"Newton with condition table report");
  end Run;

  procedure Inlined_Run
              ( file : in file_type;
                s : in Link_to_System; sols : in Solution_List;
                vrb : in integer32 := 0 ) is

    timer : Timing_Widget;
    tolres,tolerr,tolsing : double_float;
    cntfail,cntreal,cntcmplx,cntregu,cntsing,cntclus : natural32 := 0;
    maxit : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15);
    ans : character;
    verbose : boolean;

  begin
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Run 4");
    end if;
    new_line;
    Set_Parameters(maxit,tolres,tolerr,tolsing);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    tstart(timer);
    Inlined_Run
      (file,s,sols,maxit,tolres,tolerr,tolsing,cntfail,cntreal,cntcmplx,
       cntregu,cntsing,cntclus,t_err,t_rco,t_res,verbose,vrb);
    tstop(timer);
    put(file,"number of regular solutions   : ");
    put(file,cntregu,1); new_line(file);
    put(file,"number of singular solutions  : ");
    put(file,cntsing,1); new_line(file);
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
    new_line(file);
    print_times(file,timer,"Newton with condition table report");
  end Inlined_Run;

  procedure Main ( vrb : in integer32 := 0 ) is

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    nbq,len,dim : integer32 := 0;
    s : Link_to_System;
    ans : character;
    file : file_type;

  begin
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Main 1");
    end if;
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
      s := Standard_Circuit_Makers.Make_Coefficient_System(p,false);
      new_line;
      put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'n' then
        Inlined_Run(s,sols,vrb-1);
      else
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(file);
        put(file,natural32(nbq),p.all);
        new_line(file);
        Inlined_Run(file,s,sols,vrb-1);
        Close(file);
      end if;
    end if;
    Clear(p); Clear(sols); Clear(s);
  end Main;

  procedure Main ( infilename,outfilename : in string;
                   vrb : in integer32 := 0 ) is

    infile,outfile : file_type;
    p : Link_to_Poly_Sys;
    s : Link_to_System;
    sols : Solution_List;
    nbq,len,dim : integer32;
    ans : character;

  begin
    if vrb > 0
     then put_line("-> in standard_refiner_circuits.Main 2");
    end if;
    if infilename = "" then
      Main(vrb);
    else
      new_line;
      put_line("Opening the file " & infilename & " ...");
      Open_Input_File(infile,infilename);
      Standard_System_and_Solutions_io.get(infile,p,sols,"THE SOLUTIONS");
      nbq := p'last;
      len := integer32(Length_Of(sols));
      if len = 0 then
        put_line("No solutions found on file.");
      else
        dim := Head_Of(sols).n;
        new_line;
        put("Read system of "); put(nbq,1); put(" polynomials and ");
        put(len,1); put(" solutions in dimension "); put(dim,1); put_line(".");
        s := Standard_Circuit_Makers.Make_Coefficient_System(p,false);
        if outfilename = "" then
          new_line;
          put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
          if ans = 'n' then
            Inlined_Run(s,sols,vrb-1);
          else
            new_line;
            put_line("Reading the name of the output file ...");
            Read_Name_and_Create_File(outfile);
            put(outfile,natural32(nbq),p.all);
            new_line(outfile);
            Inlined_Run(outfile,s,sols,vrb-1); Close(outfile);
          end if;
        else
          new_line;
          put_line("Creating file " & outfilename & "...");
          Create_Output_File(outfile,outfilename);
          put(outfile,natural32(nbq),p.all);
          new_line(outfile);
          Inlined_Run(outfile,s,sols,vrb-1); Close(outfile);
        end if;
        Clear(p); Clear(sols); Clear(s);
      end if;
    end if;
  end Main;

end Standard_Refiner_Circuits;
