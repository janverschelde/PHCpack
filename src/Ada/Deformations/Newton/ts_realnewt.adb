with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Multprec_Floating_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Multprec_Floating_Vectors;
with Multprec_Floating_Vectors_io;       use Multprec_Floating_Vectors_io;
with Multprec_Floating_Vector_Tools;     use Multprec_Floating_Vector_Tools;
with Standard_Floating_Poly_Systems;
with Standard_Floating_Poly_Systems_io;  use Standard_Floating_Poly_Systems_io;
with Standard_Floating_Poly_SysFun;
with Standard_Floating_Jaco_Matrices;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Floating_Poly_Systems;
with Multprec_Floating_Poly_Systems_io;  use Multprec_Floating_Poly_Systems_io;
with Multprec_Floating_Poly_SysFun;

procedure ts_realnewt is

-- DESCRIPTION :
--   Interactive development of Newton's method on polynomial systems
--   with real coefficients.

  procedure Standard_Newton
              ( p : in Standard_Floating_Poly_Systems.Poly_Sys;
                x : in out Standard_Floating_Vectors.Vector ) is

    use Standard_Floating_Vectors;
    use Standard_Floating_Matrices;
    use Standard_Floating_Poly_SysFun;
    use Standard_Floating_Jaco_Matrices;

    ep : constant Eval_Poly_Sys(p'range) := Create(p);
    jm : constant Jaco_Mat(p'range,x'range) := Create(p);
    ejm : constant Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    dx : Vector(x'range);
    A : Matrix(ep'range,x'range);
    piv : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    ans : character;

  begin
    loop
      dx := Eval(ep,x); 
      put_line("current value = "); put_line(dx);
      Min(dx);
      A := Eval(ejm,x);
      lufac(A,A'last(1),piv,info);
      lusolve(A,A'last(1),piv,dx);
      put_line("increment dx = "); put_line(dx);
      Add(x,dx);
      put_line("new solution = "); put_line(x);
      put("Do you want another Newton step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Newton;

  procedure Call_Standard_Newton
              ( p : in Standard_Floating_Poly_Systems.Poly_Sys ) is

    n : constant integer32 := p'last;
    x,y : Standard_Floating_Vectors.Vector(1..n);

  begin
    new_line;
    put("Give "); put(n,1);
    put_line(" real numbers for the initial guess : "); get(x);
    y := Standard_Floating_Poly_SysFun.Eval(p,x);
    put_line("The value of p at x :"); put_line(y);
    Standard_Newton(p,x);
  end Call_Standard_Newton;

  procedure Call_Multprec_Newton
              ( p : in Multprec_Floating_Poly_Systems.Poly_Sys;
                size : in natural32 ) is

    n : constant integer32 := p'last;
    x,y : Multprec_Floating_Vectors.Vector(1..n);

  begin
    new_line;
    put("Give "); put(n,1);
    put_line(" real numbers for the initial guess : "); get(x);
    Set_Size(x,size);
    y := Multprec_Floating_Poly_SysFun.Eval(p,x);
    put_line("The value of p at x :"); put_line(y);
  end Call_Multprec_Newton;

  procedure Main is

    sp : Standard_Floating_Poly_Systems.Link_to_Poly_Sys;
    mp : Multprec_Floating_Poly_Systems.Link_to_Poly_Sys;
    deci,size : natural32 := 0;

  begin
    new_line;
    put_line("Newton's method on polynomial systems with real coefficients.");
    new_line;
    put("Give the number of decimal places : "); get(deci);
    if deci < 16 then
      size := 0;
      new_line;
      put_line("Reading a system with standard real coefficients ...");
      get(sp);
      Call_Standard_Newton(sp.all);
    else
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      new_line;
      put_line("Reading a system with multiprecision real coefficients ...");
      get(mp);
      Set_Size(mp.all,size);
      Call_Multprec_Newton(mp.all,size);
    end if;
  end Main;

begin
  Main;
end ts_realnewt;
