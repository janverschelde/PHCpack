with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Vector_Strings;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;
with Random_Conditioned_Matrices;        use Random_Conditioned_Matrices;
with Varbprec_Complex_Linear_Solvers;    use Varbprec_Complex_Linear_Solvers;
with Symbol_Table;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Jaco_Matrices;
with Multprec_Complex_Jaco_Matrices;
with Multprec_Complex_Poly_Strings;
with Random_Conditioned_Evaluations;     use Random_Conditioned_Evaluations;
with Varbprec_Polynomial_Evaluations;    use Varbprec_Polynomial_Evaluations;
with Varbprec_Complex_Newton_Steps;

package body Random_Conditioned_Root_Problems is

  procedure Standard_Test
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                z : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Jaco_Matrices;
    use Varbprec_Complex_Newton_Steps;

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
    Clear(jp);
  end Standard_Test;

  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Jaco_Matrices;
    use Varbprec_Complex_Newton_Steps;

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
    Clear(jp);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Jaco_Matrices;
    use Varbprec_Complex_Newton_Steps;

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
    Clear(jp);
  end QuadDobl_Test;

  procedure Multprec_Test
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                z : in out Multprec_Complex_Vectors.Vector ) is

    use Multprec_Complex_Jaco_Matrices;
    use Varbprec_Complex_Newton_Steps;

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
    Clear(jp);
  end Multprec_Test;

  procedure Standard_Conditioned_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

    p : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : Standard_Complex_Vectors.Vector(p'range);
    jm : constant Standard_Complex_Matrices.Matrix(p'range,x'range)
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
   -- put_line("The polynomial system : "); put_line(p);
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
    loss := Varbprec_Complex_Newton_Steps.Minimum(loss_jac,loss_eva);
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

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : DoblDobl_Complex_Vectors.Vector(p'range);
    jm : constant DoblDobl_Complex_Matrices.Matrix(p'range,x'range)
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
   -- put_line("The polynomial system : "); put_line(p);
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
    loss := Varbprec_Complex_Newton_Steps.Minimum(loss_jac,loss_eva);
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

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : QuadDobl_Complex_Vectors.Vector(p'range);
    jm : constant QuadDobl_Complex_Matrices.Matrix(p'range,x'range)
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
   -- put_line("The polynomial system :"); put_line(p);
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
    loss := Varbprec_Complex_Newton_Steps.Minimum(loss_jac,loss_eva);
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
   -- put_line("The polynomial system :"); put_line(p);
    rco := Inverse_Condition_Number(p,x);
    if Equal(rco,0.0) then
      loss_eva := -2**30;
    else
      declare
        mp_log10rco : Floating_Number := log10(rco);
        st_log10rco : constant double_float := Round(mp_log10rco);
      begin
        Clear(mp_log10rco);
        loss_eva := integer32(st_log10rco);
      end;
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    loss := Varbprec_Complex_Newton_Steps.Minimum(loss_jac,loss_eva);
    put("-> Estimated total loss : "); put(loss,1); new_line;
    put("Give the wanted number of decimal places : "); get(want);
    if precision + loss >= want then
      put("Current multiprecision suffices to meet "); put(want,1);
      put_line(" accurate decimal places.");
      Multprec_Test(p,x);
    else
      put("Current multiprecision does not suffice for "); put(want,1);
      put_line(" accurate decimal places.");
    end if;
    Multprec_Complex_Matrices.Clear(jm);
  end Multprec_Conditioned_Test;

  procedure Multprec_Conditioned_Root_Problem
              ( n,d,m,c,prcn : in natural32;
                cffsz,pntsz,close,condjm : in double_float;
                f : out Link_to_Array_of_Strings; z : out Link_to_String ) is

    p : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : Multprec_Complex_Vectors.Vector(p'range);
    jm : Multprec_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    sz : constant natural32 := Decimal_to_Size(prcn);

  begin
    Random_Conditioned_Jacobian_Evaluation
      (n,d,m,c,sz,cffsz,pntsz,close,jm,p,x);
   -- put_line("The polynomials :"); put_line(p);
   -- put_line("The initial approximation : "); put_line(x);
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(integer32(n)));
    declare
      strx : constant string := Multprec_Complex_Vector_Strings.Write(x);
      strp : Array_of_Strings(1..integer(p'last));
    begin
      z := new string'(strx);
      for i in p'range loop
        declare
          strpi : constant string
                := Multprec_Complex_Poly_Strings.Write(p(i));
        begin
          strp(integer(i)) := new string'(strpi);
        end;
      end loop;
      f := new Array_of_Strings'(strp);
    end;
    Multprec_Complex_Matrices.Clear(jm);
  end Multprec_Conditioned_Root_Problem;

  procedure Random_Conditioned_Parameters
              ( n,d,m : out natural32;
                condjm,cffsize,pntsize,close,condfz : out double_float ) is
  begin
    condjm := 0.0; cffsize := 0.0; pntsize := 0.0; close := 0.0;
    new_line;
    put_line("First part, condition number of Jacobian matrix :");
    put("Give the condition of the Jacobian matrix : "); get(condjm);
    new_line;
    n := 0; d:= 0; m := 0;
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
    condfz := cffsize*(pntsize**integer(d))/close;
    put("Predicted condition number : "); put(condfz,3); new_line;
  end Random_Conditioned_Parameters;

  procedure Random_Conditioned_Root_Problem ( preclvl : in character ) is

    n,d,m : natural32 := 0;
    condjm,cffsize,pntsize,close,condfz : double_float := 0.0;

  begin
    Random_Conditioned_Parameters(n,d,m,condjm,cffsize,pntsize,close,condfz);
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

  function Maximum ( a,b : integer32 ) return integer32 is
  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Maximum;

  procedure Random_Conditioned_Root_Problem
              ( f : out Link_to_Array_of_Strings;
                z : out Link_to_String ) is

    n,d,m : natural32 := 0;
    condjm,cffsize,pntsize,close,condfz : double_float := 0.0;
    jmloss,fzloss,maxloss : integer32;
    precision : natural32;

  begin
    Random_Conditioned_Parameters(n,d,m,condjm,cffsize,pntsize,close,condfz);
    jmloss := integer32(log10(condjm));
    put("-> expected loss of linear system solving : ");
    put(jmloss,1); new_line;
    fzloss := integer32(log10(condfz));
    put("-> expected loss of polynomial evaluation : ");
    put(fzloss,1); new_line;
    maxloss := Maximum(jmloss,fzloss);
    put("-> maximum loss of accuracy : "); put(maxloss,1); new_line;
    precision := 2*natural32(maxloss);
    put("Precision of generating conditioned root problem : "); 
    put(precision,1); new_line;
    Multprec_Conditioned_Root_Problem
      (n,d,m,0,precision,cffsize,pntsize,close,condjm,f,z);
  end Random_Conditioned_Root_Problem;

end Random_Conditioned_Root_Problems;
