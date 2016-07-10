with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Standard_Homotopy;
with Standard_Root_Refiners;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Dense_Series_Vectors;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_Systems;
with Series_and_Polynomials;
with Series_and_Polynomials_io;          use Series_and_Polynomials_io;
with Series_and_Homotopies;
with Series_and_Predictors;

procedure ts_serpath is

-- DESCRIPTION :
--   Developing path tracers with power series.

  procedure Standard_Correct
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                t : in double_float;
                sol : in out Standard_Complex_Vectors.Vector;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Applies a couple of Newton iterations to correct the solution sol
  --   of the system at t.

    p : Standard_Complex_Poly_Systems.Poly_Sys(hom'range)
      := Series_and_Homotopies.Eval(hom,t);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(hom'range,sol'range)
       := Standard_Complex_Jaco_Matrices.Create(p);
    f : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := Standard_Complex_Poly_SysFun.Create(p);
    jf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat(hom'range,sol'range)
       := Standard_Complex_Jaco_Matrices.Create(jm);
    err,rco,res : double_float;

    use Standard_Root_Refiners;

  begin
    for k in 1..3 loop
      Standard_Newton_Step(f,jf,sol,err,rco,res);
      put("  err :"); put(err,3);
      put("  rco :"); put(rco,3);
      put("  res :"); put(res,3); new_line;
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_SysFun.Clear(f);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Standard_Correct;

  procedure Standard_Track
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom,
  --   in standard double precision.

    wrk : Standard_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : integer32 := 4;
    srv : Standard_Dense_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Dense_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0E-12;
    tolres : constant double_float := 1.0E-5;
    t,step : double_float := 0.0;
    max_steps : constant natural32 := 500;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;

  begin
    Standard_Series_Poly_Systems.Copy(hom,wrk);
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(nit,wrk,wrk_sol,srv,eva,verbose);
      new_line;
      put_line("Setting the step size based on the power series ...");
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres,verbose);
      put("The computed step size : "); put(step,3);
      if step > 0.1
       then step := 0.1;
      end if;
      if t + step <= 1.0 then
        t := t + step;
      else
        step := 1.0 - t;
        t := 1.0; 
      end if;
      put(" t = "); put(t,3);  new_line;
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,step);
      put_line("Shifting the polynomial system ...");
      Standard_Series_Poly_Systems.Clear(wrk);
      wrk := Series_and_Homotopies.Shift(hom,t);
      put_line("Correcting the solution ...");
      Standard_Correct(hom,t,wrk_sol);
      exit when (t = 1.0);
    end loop;
    put_line("The solution vector :"); put_line(wrk_sol);
    Standard_Series_Poly_Systems.Clear(wrk);
  end Standard_Track;

  procedure DoblDobl_Track
              ( h : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tracks one path starting at the solution s using the homotopy h,
  --   in double double precision.

  begin
    null;
  end DoblDobl_Track;

  procedure QuadDobl_Track
              ( h : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out Quaddobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tracks one path starting at the solution s using the homotopy h,
  --   in quad double precision.

  begin
    null;
  end QuadDobl_Track;

  procedure Standard_Test
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in standard double precision.

    use Standard_Complex_Solutions;

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,true);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;

  begin
    put_line("The homotopy system :"); put_line(h);
    put_line("The series system :"); put(s,1);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      Standard_Track(s,ls.all);
      put("Continue to the next path ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in double double precision.

    use DoblDobl_Complex_Solutions;

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,true);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;

  begin
    put_line("The homotopy system :"); put_line(h);
    put_line("The series system :"); put(s,1);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      DoblDobl_Track(s,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in quad double precision.

    use QuadDobl_Complex_Solutions;

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,true);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;

  begin
    put_line("The homotopy system :"); put_line(h);
    put_line("The series system :"); put(s,1);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      QuadDobl_Track(s,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := 2;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and the start solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    Standard_Homotopy.Create(target.all,start.all,k,gamma);
    new_line;
    Standard_Test(target'last,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := 2;
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    DoblDobl_Homotopy.Create(target.all,start.all,k,gamma);
    new_line;
    DoblDobl_Test(target'last,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := 2;
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    QuadDobl_Homotopy.Create(target.all,start.all,k,gamma);
    new_line;
    QuadDobl_Test(target'last,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision
  --   and then launches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpath;
