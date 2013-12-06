with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_Residual_Evaluations;      use Multprec_Residual_Evaluations;
with Multprec_Root_Refiners;             use Multprec_Root_Refiners;

procedure ts_rootrefi is

-- DESCRIPTION :
--   This routine facilitates interactive testing of the root refiners.

  procedure Call_Standard_Newton
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

    epsxa,epsfa : double_float;
    p_eval : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
           := Standard_Complex_Poly_SysFun.Create(p); 
    jac : Standard_Complex_Jaco_Matrices.Jaco_Mat(p'range,p'range)
        := Standard_Complex_Jaco_Matrices.Create(p);
    jac_eval : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat(p'range,p'range)
             := Standard_Complex_Jaco_Matrices.Create(jac);
    numit : natural32 := 0;
    max : constant natural32 := 5;
    fail : boolean;
    ls : constant Standard_Complex_Solutions.Link_to_Solution := Head_Of(sols);
    zero : Standard_Complex_Solutions.Solution(p'last) := ls.all;

  begin
    epsxa := 1.0E-14;
    epsfa := 1.0E-14;
    Reporting_Newton
      (file,p_eval,jac_eval,zero,epsxa,epsfa,numit,max,fail);
    put_line("The refined solution : "); put(zero); new_line;
    if fail
     then put_line("failed to reach desired precision");
     else put_line("converged to desired precision");
    end if;
    Standard_Complex_Poly_SysFun.Clear(p_eval);
    Standard_Complex_Jaco_Matrices.Clear(jac);
    Standard_Complex_Jaco_Matrices.Clear(jac_eval);
  end Call_Standard_Newton;

  procedure Call_Multprec_Newton
               ( file : in file_type;
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Multprec_Complex_Solutions.Solution_List ) is

    epsxa,epsfa : Floating_Number;
    p_eval : Multprec_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
           := Multprec_Complex_Poly_SysFun.Create(p); 
    jac : Multprec_Complex_Jaco_Matrices.Jaco_Mat(p'range,p'range)
        := Multprec_Complex_Jaco_Matrices.Create(p);
    jac_eval : Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat(p'range,p'range)
             := Multprec_Complex_Jaco_Matrices.Create(jac);
    numit,deci,size : natural32 := 0;
    max : constant natural32 := 5;
    fail : boolean;
    ls : Multprec_Complex_Solutions.Link_to_Solution := Head_Of(sols);
    zero : Multprec_Complex_Solutions.Solution(p'last);

  begin
    Copy(ls.all,zero); -- the "zero := ls.all" left empty zero.v
    put("Give the number of decimal places : "); get(deci);
    size := Decimal_to_Size(deci);
    put("The size of the numbers : "); put(size,1); new_line;
    Set_Size(sols,size);
    put("Give tolerance for error : "); get(epsxa);
    put("Give tolerance for residual : "); get(epsfa);
    Reporting_Newton
      (file,p_eval,jac_eval,zero,epsxa,epsfa,numit,max,fail);
    put_line("The refined solution : "); put(zero); new_line;
    if fail
     then put_line("failed to reach desired precision");
     else put_line("converged to desired precision");
    end if;
    Multprec_Complex_Poly_SysFun.Clear(p_eval);
    Multprec_Complex_Jaco_Matrices.Clear(jac);
    Multprec_Complex_Jaco_Matrices.Clear(jac_eval);
  end Call_Multprec_Newton;

  procedure Call_Standard_Root_Refiner
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    numit : natural32 := 0;
    max : constant natural32 := 5;
    deflate : boolean := true;

  begin
    epsxa := 1.0E-14;
    epsfa := 1.0E-14;
    tolsing := 1.0E-08;
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,deflate,true);
  end Call_Standard_Root_Refiner;

  procedure Call_Multprec_Root_Refiner
               ( file : in file_type;
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Multprec_Complex_Solutions.Solution_List ) is

    epsxa,epsfa,tolsing : Floating_Number;
    numit,deci,size : natural32 := 0;
    max : constant natural32 := 5;

  begin
    put("Give the number of decimal places : "); get(deci);
    size := Decimal_to_Size(deci);
    put("The size of the numbers : "); put(size,1); new_line;
    Set_Size(sols,size);
    put("Give tolerance for error : "); get(epsxa);
    put("Give tolerance for residual : "); get(epsfa);
    tolsing := Create(1.0E-08);
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,true,true);
  end Call_Multprec_Root_Refiner;

  procedure Test_Standard_Root_Refiner is

  -- DESCRIPTION :
  --   Test of root refining on list of solutions as standard vectors.

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Test on refining roots as standard complex vectors.");
    new_line;
    get(lp);
    put_line("The system : "); put(lp.all);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    put(file,natural32(lp'last),lp.all);
    new_line;
    Read(sols);
    new_line;
    new_line(file);
    put_line(file,"THE INITIAL SOLUTIONS : ");
    put(file,sols);
    Call_Standard_Root_Refiner(file,lp.all,sols);
  end Test_Standard_Root_Refiner;

  function Multiple_Root_Test_Polynomial
             ( n,i,k : natural32 ) return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns x_i^k = 0 as an n-variable polynomial.

    use Standard_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := k;
    res := Create(t);
    return res;
  end Multiple_Root_Test_Polynomial;

  function Multiple_Root_Test_System
              ( n : natural32; k : Standard_Natural_Vectors.Vector ) 
              return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system x_i^k(i) = 0, for i in 1..n.

    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));

  begin
    for i in 1..integer32(n) loop
      res(i) := Multiple_Root_Test_Polynomial(n,natural32(i),k(i));
    end loop;
    return res;
  end Multiple_Root_Test_System;

  function Read_Start_Solution ( n : natural32 )
             return Standard_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Reads a complex vector of length n and returns the list of
  --   solutions containing this vector.

    res : Standard_Complex_Solutions.Solution_List;
    sol : Standard_Complex_Solutions.Solution(integer32(n));

  begin
    put("Give "); put(n,1);
    put_line(" complex numbers as components of the start solution :");
    sol.t := Create(0.0);
    sol.m := 1;
    sol.err := 0.0;
    sol.rco := 0.0;
    sol.res := 0.0;
    get(sol.v);
    Standard_Complex_Solutions.Add(res,sol);
    return res;
  end Read_Start_Solution;

  procedure Test_Multprec_Root_Refiner is

  -- DESCRIPTION :
  --   Test of root refining on list of solutions as standard vectors.

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    stsols : Standard_Complex_Solutions.Solution_List;
    mpsols : Multprec_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Test on refining roots as multi-precision complex vectors.");
    new_line;
    get(lp);
    put_line("The system : "); put(lp.all);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    put(file,natural32(lp'last),lp.all);
    new_line;
   -- Read(stsols);
   -- new_line(file);
   -- put_line(file,"THE SOLUTION IN STANDARD PRECISION : ");
   -- put(file,stsols);
   -- mpsols := Create(stsols);
    Read(mpsols);
    new_line;
    new_line(file);
    put_line(file,"THE INITIAL SOLUTIONS : ");
    put(file,Length_Of(mpsols),natural32(lp'last),mpsols);
    declare
      mp : Multprec_Complex_Poly_Systems.Poly_Sys(lp'range) := Convert(lp.all);
    begin
      Call_Multprec_Root_Refiner(file,mp,mpsols);
    end;
  end Test_Multprec_Root_Refiner;

  procedure Test_Multprec_Residual_Evaluator is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mpsols : Multprec_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Test on evaluating residuals with multi-precision arithmetic.");
    new_line;
    get(lp);
    put_line("The system : "); put(lp.all);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    put(file,natural32(lp'last),lp.all);
    new_line;
    Read(mpsols);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(mpsols),natural32(lp'last),mpsols);
    declare
      mp : Multprec_Complex_Poly_Systems.Poly_Sys(lp'range) := Convert(lp.all);
      mp_eval : Multprec_Complex_Poly_SysFun.Eval_Poly_Sys(mp'range)
              := Multprec_Complex_Poly_SysFun.Create(mp);
      deci,size : natural32;
    begin
      put("Give the number of decimal places : "); get(deci);
      size := Decimal_to_Size(deci);
      put("The size of the numbers : "); put(size,1); new_line;
      Set_Size(mpsols,size); 
      put_line(file,"THE RESIDUALS :");
      Residuals(file,mp_eval,mpsols);
      Multprec_Complex_Poly_Systems.Clear(mp);
      Multprec_Complex_Poly_SysFun.Clear(mp_eval);
    end;
  end Test_Multprec_Residual_Evaluator;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of root refiners.");
    loop
      new_line;
      put_line("Choose one of the following :                              ");
      put_line("  0. exit this program.                                    ");
      put_line("  1. Test root refiner for standard complex numbers.       ");
      put_line("  2. Test root refiner for multi-precision complex numbers.");
      put_line("  3. Evaluate residuals with multi-precision arithmetic.   ");
      put("Type 1, 2, or 3 to select : "); Ask_Alternative(ans,"0123");
      exit when ans = '0';
      case ans is
        when '1' => Test_Standard_Root_Refiner;
        when '2' => Test_Multprec_Root_Refiner;
        when '3' => Test_Multprec_Residual_Evaluator;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_rootrefi;
