with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
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
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Vector_Tools;
with Multprec_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;
with Random_Conditioned_Matrices;        use Random_Conditioned_Matrices;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Jaco_Matrices;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Strings;      use Multprec_Complex_Poly_Strings;
with Multprec_Complex_Jaco_Matrices;
with Random_Conditioned_Evaluations;     use Random_Conditioned_Evaluations;
with Varbprec_Complex_Linear_Solvers;    use Varbprec_Complex_Linear_Solvers;
with Varbprec_Polynomial_Evaluations;    use Varbprec_Polynomial_Evaluations;
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

  function Multprec_Initial_Approximation
              ( n : integer32; size : natural32 )
              return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user for n coordinates for an initial approximation
  --   to run Newton's method or generates a random vector.

    res : Multprec_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    new_line;
    put("Use random vector in Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      res := Multprec_Random_Vectors.Random_Vector(1,n,size);
    else
      put("Reading "); put(n,1); put(" complex numbers ...");
      for i in 1..n loop
        put("x("); put(i,1); put(") : "); get(res(i));
      end loop;
    end if;
    Multprec_Complex_Vector_Tools.Set_Size(res,size);
    return res;
  end Multprec_Initial_Approximation;

  procedure Standard_Test
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                z : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps on p,
  --   in standard double precision, starting at z
  --   as the initial approximation for a solution of p.

    use Standard_Complex_Jaco_Matrices;

    jp : Jaco_Mat(p'range,p'range) := Create(p);
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

  procedure Standard_Test
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps on p,
  --   in standard double precision, prompting the user for
  --   an initial approximation for a root, 
  --   or generating just a random vector.

    z : Standard_Complex_Vectors.Vector(p'range)
      := Standard_Initial_Approximation(p'last);

  begin
    Standard_Test(p,z);
  end Standard_Test;

  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in double double precision, starting at z 
  --   as an initial approximation for a solution of p.

    use DoblDobl_Complex_Jaco_Matrices;

    jp : Jaco_Mat(p'range,p'range) := Create(p);
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

  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in double double precision, prompting the user for an initial
  --   approximation for a solution, or generating a random vector.

    z : DoblDobl_Complex_Vectors.Vector(p'range)
      := DoblDobl_Initial_Approximation(p'last);

  begin
    DoblDobl_Test(p,z);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in quad double precision, starting at z
  --   as an initial approximation for a solution.

    use QuadDobl_Complex_Jaco_Matrices;

    jp : Jaco_Mat(p'range,p'range) := Create(p);
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

  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in quad double precision, prompting the user for an initial
  --   approximation for a solution or generating a random vector.

    z : QuadDobl_Complex_Vectors.Vector(p'range)
      := QuadDobl_Initial_Approximation(p'last);

  begin
    QuadDobl_Test(p,z);
  end QuadDobl_Test;

  procedure Multprec_Test
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                size : in natural32;
                z : in out Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in arbitrary multiprecision with numbers of the given size,
  --   using z as initial approximation for a solution of p.

    use Multprec_Complex_Jaco_Matrices;

    jp : Jaco_Mat(p'range,p'range) := Create(p);
    fz : Multprec_Complex_Vectors.Vector(p'range);
    jpz : Multprec_Complex_Matrices.Matrix(p'range,p'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    sysrco,evarco,err : Floating_Number;
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
      Multprec_Complex_Matrices.Clear(jpz);
      Multprec_Complex_Vectors.Clear(fz);
      Clear(sysrco); Clear(evarco); Clear(err);
    end loop;
  end Multprec_Test;

  procedure Multprec_Test
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                size : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in arbitrary multiprecision with numbers of the given size,
  --   prompting the user to provide an initial approximation for
  --   a solution, or generating a random vector.

    z : Multprec_Complex_Vectors.Vector(p'range)
      := Multprec_Initial_Approximation(p'last,size);

  begin
    Multprec_Test(p,size,z);
  end Multprec_Test;

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
  --   in quad double precision.

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    QuadDobl_Test(lp.all);
  end QuadDobl_Test_on_Given_System;

  procedure Multprec_Test_on_Given_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and runs the variable precision Newton steps
  --   in arbitrary multiprecision.

    use Multprec_Complex_Poly_Systems;
    file : file_type;
    n,m : natural32;
    ls : Link_to_Array_of_Strings;
    deci,size : natural32 := 0;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci); skip_line;
    size := Decimal_to_Size(deci);
    new_line;
    put_line("Reading a polynomial system ...");
    Read_Name_and_Open_File(file);
    get(file,integer(n),integer(m),ls);
    if Symbol_Table.Number < m
     then Symbol_Table.Init(m);
    end if;
    declare
      p : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(n))
        := Parse(m,size,ls.all);
    begin
      Multprec_Test(p,size);
    end;
  end Multprec_Test_on_Given_System;

  function Minimum ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the minimum of a and b.

  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Minimum;

  procedure Standard_Conditioned_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+16.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : Standard_Complex_Vectors.Vector(p'range);
    jm : Standard_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    rco : double_float;
    precision : constant integer32 := 16;
    loss_jac,loss_eva,loss,want : integer32 := 0;

  begin
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if rco = 0.0
     then loss_eva := -2**30;
     else loss_eva := integer32(log10(rco));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    loss := Minimum(loss_jac,loss_eva);
    put("-> Estimated total loss : "); put(loss,1); new_line;
    put("Give the wanted number of decimal places : "); get(want);
    if precision + loss >= want then
      put("Double precision suffices to meet "); put(want,1);
      put_line(" accurate decimal places.");
      Standard_Test(p,x);
    else
      put("Double precision does not suffice for "); put(want,1);
      put_line(" accurate decimal places.");
    end if;
  end Standard_Conditioned_Test;

  procedure DoblDobl_Conditioned_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with double double complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+32.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : DoblDobl_Complex_Vectors.Vector(p'range);
    jm : DoblDobl_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    rco : double_double;
    precision : constant integer32 := 32;
    loss_jac,loss_eva,loss,want : integer32 := 0;

  begin
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Is_Zero(rco)
     then loss_eva := -2**30;
     else loss_eva := integer32(to_double(log10(rco)));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation : ");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    loss := Minimum(loss_jac,loss_eva);
    put("-> Estimated total loss : "); put(loss,1); new_line;
    put("Give the wanted number of decimal places : "); get(want);
    if precision + loss >= want then
      put("Double double precision suffices to meet "); put(want,1);
      put_line(" accurate decimal places.");
      DoblDobl_Test(p,x);
    else
      put("Double double precision does not suffice for "); put(want,1);
      put_line(" accurate decimal places.");
    end if;
  end DoblDobl_Conditioned_Test;

  procedure QuadDobl_Conditioned_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with quad double complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+64.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : QuadDobl_Complex_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    rco : quad_double;
    precision : constant integer32 := 64;
    loss_jac,loss_eva,loss,want : integer32 := 0;

  begin
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Is_Zero(rco)
     then loss_eva := -2**30;
     else loss_eva := integer32(to_double(log10(rco)));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation : ");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    loss := Minimum(loss_jac,loss_eva);
    put("-> Estimated total loss : "); put(loss,1); new_line;
    put("Give the wanted number of decimal places : "); get(want);
    if precision + loss >= want then
      put("Quad double precision suffices to meet "); put(want,1);
      put_line(" accurate decimal places.");
      QuadDobl_Test(p,x);
    else
      put("Quad double precision does not suffice for "); put(want,1);
      put_line(" accurate decimal places.");
    end if;
  end QuadDobl_Conditioned_Test;

  procedure Multprec_Conditioned_Test
              ( n,d,m,c,sz : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with multiprecision complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed, provided the value sz for the precision
  --   is sufficiently high.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   sz       size of the numbers in the working precision;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : Multprec_Complex_Vectors.Vector(p'range);
    jm : Multprec_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    rco : Floating_Number;
    precision : constant integer32 := integer32(Size_to_Decimal(sz));
    loss_jac,loss_eva,loss,want : integer32 := 0;

  begin
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,sz,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Equal(rco,0.0) then
      loss_eva := -2**30;
    else
      declare
        mp_log10rco : Floating_Number := log10(rco);
        st_log10rco : double_float := Round(mp_log10rco);
      begin
        loss_eva := integer32(st_log10rco);
      end;
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    loss := Minimum(loss_jac,loss_eva);
    put("-> Estimated total loss : "); put(loss,1); new_line;
    put("Give the wanted number of decimal places : "); get(want);
    if precision + loss >= want then
      put("Current multiprecision suffices to meet "); put(want,1);
      put_line(" accurate decimal places.");
      Multprec_Test(p,sz,x);
    else
      put("Current multiprecision does not suffice for "); put(want,1);
      put_line(" accurate decimal places.");
    end if;
  end Multprec_Conditioned_Test;

  procedure Random_Conditioned_Root_Problem ( preclvl : in character ) is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of the numerical problem
  --   to set the condition of the Jacobian matrix and
  --   the condition of the polynomial evaluation problem.
  --   The level of precision is indicated by the character preclvl.

    n,d,m,deci,size : natural32 := 0;
    condjm,cffsize,pntsize,close,cond : double_float := 0.0;

  begin
    new_line;
    put_line("First part, condition number of Jacobian matrix :");
    put("Give the condition of the Jacobian matrix : "); get(condjm);
    new_line;
    put_line("Second part, parameters of the random polynomials :");
    put("Give number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    new_line;
    put_line("Third part, factors in numerical condition of evaluation :");
    put("Give magnitude of the coefficients : "); get(cffsize);
    put("Give magnitude of the coordinates of the point : "); get(pntsize);
    put("Give closeness to a root : "); get(close);
    cond := cffsize*(pntsize**integer(d))/close;
    put("Predicted condition number : "); put(cond,3); new_line;
    new_line;
    case preclvl is
      when '0' => 
        Standard_Conditioned_Test(n,d,m,0,cffsize,pntsize,close,condjm);
      when '1' =>
        DoblDobl_Conditioned_Test(n,d,m,0,cffsize,pntsize,close,condjm);
      when '2' =>
        QuadDobl_Conditioned_Test(n,d,m,0,cffsize,pntsize,close,condjm);
      when '3' =>
        declare
          deci,size : natural32 := 0;
        begin
          put("Give the number of decimal places : "); get(deci);
          size := Decimal_to_Size(deci);
          Multprec_Conditioned_Test(n,d,m,0,size,cffsize,pntsize,close,condjm);
        end;
      when others => null;
    end case;
  end Random_Conditioned_Root_Problem;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then either calls the procedure to generate a random
  --   conditioned root problems or calls the procedure that will 
  --   test variable precision Newton's method on a given system.

    precision,random : character;

  begin
    new_line;
    put_line("Testing variable precision Newton steps ...");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision'");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(precision,"0123");
    new_line;
    put("Generate random conditioned root problem ? (y/n) ");
    Ask_Yes_or_No(random);
    if random = 'y' then
      new_line;
      put_line("Generating a random conditioned root problem ...");
      Random_Conditioned_Root_Problem(precision);
    else
      case precision is 
        when '0' => Standard_Test_on_Given_System;
        when '1' => DoblDobl_Test_on_Given_System;
        when '2' => QuadDobl_Test_on_Given_System;
        when '3' => Multprec_Test_on_Given_System;
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_vmpnewt;
