with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Jaco_Matrices;
with Varbprec_Complex_Newton_Steps;      use Varbprec_Complex_Newton_Steps;

procedure ts_vmpnewt is

-- DESCRIPTION :
--   Development of variable precision Newton's method.

  function Standard_Initial_Approximation
              ( n : integer32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user for n coordinates for an initial approximation
  --   to run Newton's method or generates a random vector.

    res : Standard_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    new_line;
    put("Use random vector in Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      res := Standard_Random_Vectors.Random_Vector(1,n);
    else
      put("Reading "); put(n,1); put(" complex numbers ...");
      for i in 1..n loop
        put("x("); put(i,1); put(") : "); get(res(i));
      end loop;
    end if;
    return res;
  end Standard_Initial_Approximation;

  function DoblDobl_Initial_Approximation
              ( n : integer32 ) return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user for n coordinates for an initial approximation
  --   to run Newton's method or generates a random vector.

    res : DoblDobl_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    new_line;
    put("Use random vector in Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      res := DoblDobl_Random_Vectors.Random_Vector(1,n);
    else
      put("Reading "); put(n,1); put(" complex numbers ...");
      for i in 1..n loop
        put("x("); put(i,1); put(") : "); get(res(i));
      end loop;
    end if;
    return res;
  end DoblDobl_Initial_Approximation;

  function QuadDobl_Initial_Approximation
              ( n : integer32 ) return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user for n coordinates for an initial approximation
  --   to run Newton's method or generates a random vector.

    res : QuadDobl_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    new_line;
    put("Use random vector in Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      res := QuadDobl_Random_Vectors.Random_Vector(1,n);
    else
      put("Reading "); put(n,1); put(" complex numbers ...");
      for i in 1..n loop
        put("x("); put(i,1); put(") : "); get(res(i));
      end loop;
    end if;
    return res;
  end QuadDobl_Initial_Approximation;

  procedure Standard_Test
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps on p,
  --   in standard double precision.

    use Standard_Complex_Jaco_Matrices;

    jp : Jaco_Mat(p'range,p'range) := Create(p);
    z : Standard_Complex_Vectors.Vector(p'range)
      := Standard_Initial_Approximation(p'last);
    fz : Standard_Complex_Vectors.Vector(p'range);
    jpz : Standard_Complex_Matrices.Matrix(p'range,p'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    sysrco,evarco,err : double_float;
    sysloss,evaloss : integer32;
    ans : character;

  begin
    loop
      Estimate_Loss_in_Newton_Step
        (p,jp,z,jpz,piv,fz,sysrco,evarco,sysloss,evaloss);
      put_line("The system evaluated at the current solution :");
      put_line(fz);
      put("linear system rco : "); put(sysrco,3); new_line;
      put("   evaluation rco : "); put(evarco,3); new_line;
      put("estimated loss of linear system solving : ");
      put(sysloss,1); new_line;
      put("estimated loss of polynomial evaluation : ");
      put(evaloss,1); new_line;
      do_Newton_Step(z,jpz,piv,fz,err);
      put_line("The current solution vector : "); put_line(z);
      put("magnitude of correction on root : "); put(err,3); new_line;
      put("Continue ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in double double precision.

    use DoblDobl_Complex_Jaco_Matrices;

    jp : Jaco_Mat(p'range,p'range) := Create(p);
    z : DoblDobl_Complex_Vectors.Vector(p'range)
      := DoblDobl_Initial_Approximation(p'last);
    fz : DoblDobl_Complex_Vectors.Vector(p'range);
    jpz : DoblDobl_Complex_Matrices.Matrix(p'range,p'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    sysrco,evarco,err : double_double;
    sysloss,evaloss : integer32;
    ans : character;

  begin
    loop
      Estimate_Loss_in_Newton_Step
        (p,jp,z,jpz,piv,fz,sysrco,evarco,sysloss,evaloss);
      put_line("The system evaluated at the current solution :");
      put_line(fz);
      put("linear system rco : "); put(sysrco,3); new_line;
      put("   evaluation rco : "); put(evarco,3); new_line;
      put("estimated loss of linear system solving : ");
      put(sysloss,1); new_line;
      put("estimated loss of polynomial evaluation : ");
      put(evaloss,1); new_line;
      do_Newton_Step(z,jpz,piv,fz,err);
      put_line("The current solution vector : "); put_line(z);
      put("magnitude of correction on root : "); put(err,3); new_line;
      put("Continue ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in double double precision.

    use QuadDobl_Complex_Jaco_Matrices;

    jp : Jaco_Mat(p'range,p'range) := Create(p);
    z : QuadDobl_Complex_Vectors.Vector(p'range)
      := QuadDobl_Initial_Approximation(p'last);
    fz : QuadDobl_Complex_Vectors.Vector(p'range);
    jpz : QuadDobl_Complex_Matrices.Matrix(p'range,p'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    sysrco,evarco,err : quad_double;
    sysloss,evaloss : integer32;
    ans : character;

  begin
    loop
      Estimate_Loss_in_Newton_Step
        (p,jp,z,jpz,piv,fz,sysrco,evarco,sysloss,evaloss);
      put_line("The system evaluated at the current solution :");
      put_line(fz);
      put("linear system rco : "); put(sysrco,3); new_line;
      put("   evaluation rco : "); put(evarco,3); new_line;
      put("estimated loss of linear system solving : ");
      put(sysloss,1); new_line;
      put("estimated loss of polynomial evaluation : ");
      put(evaloss,1); new_line;
      do_Newton_Step(z,jpz,piv,fz,err);
      put_line("The current solution vector : "); put_line(z);
      put("magnitude of correction on root : "); put(err,3); new_line;
      put("Continue ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Test;

  procedure Standard_Test_on_Given_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and runs the variable precision Newton steps
  --   in standard double precision.

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    Standard_Test(lp.all);
  end Standard_Test_on_Given_System;

  procedure DoblDobl_Test_on_Given_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and runs the variable precision Newton steps
  --   in double double precision.

    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    DoblDobl_Test(lp.all);
  end DoblDobl_Test_on_Given_System;

  procedure QuadDobl_Test_on_Given_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and runs the variable precision Newton steps
  --   in double double precision.

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    QuadDobl_Test(lp.all);
  end QuadDobl_Test_on_Given_System;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision.

    ans : character;

  begin
    new_line;
    put_line("Testing variable precision Newton steps ...");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision.");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is 
      when '0' => Standard_Test_on_Given_System;
      when '1' => DoblDobl_Test_on_Given_System;
      when '2' => QuadDobl_Test_on_Given_System;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_vmpnewt;
