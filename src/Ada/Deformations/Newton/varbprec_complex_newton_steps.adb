with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Multprec_DoblDobl_Convertors;
with Multprec_QuadDobl_Convertors;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Vector_Tools;
with Standard_Complex_Vector_Strings;
with DoblDobl_Complex_Vector_Strings;
with QuadDobl_Complex_Vector_Strings;
with Multprec_Complex_Vector_Strings;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;
with Symbol_Table;
with Standard_Complex_Poly_Strings;
with Standard_Complex_Laur_Strings;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Laur_SysFun;
with DoblDobl_Complex_Poly_Strings;
with DoblDobl_Complex_Laur_Strings;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Laur_Strings;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Laur_SysFun;
with Multprec_Complex_Poly_Strings;
with Multprec_Complex_Laur_Strings;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Laur_SysFun;
with Varbprec_Complex_Linear_Solvers;    use Varbprec_Complex_Linear_Solvers;
with Varbprec_Polynomial_Evaluations;    use Varbprec_Polynomial_Evaluations;
-- for debugging :
--with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;

package body Varbprec_Complex_Newton_Steps is

  procedure Estimate_Loss_in_Newton_Step
              ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                z : in Standard_Complex_Vectors.Vector;
                jfz : out Standard_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Standard_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_float;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : double_float;

  begin
    jfz := Standard_Complex_Jaco_Matrices.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,absfz,denrco,fzrco,fz);
    if absfz > 0.1 then
      if fzrco = 0.0
       then fzloss := -2**30;
       else fzloss := integer32(log10(fzrco));
      end if;
    else
      fzloss := -integer32(log10(denrco));
    end if;
  end Estimate_Loss_in_Newton_Step;

  procedure Estimate_Loss_in_Newton_Step
              ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
                z : in DoblDobl_Complex_Vectors.Vector;
                jfz : out DoblDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out DoblDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_double;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : double_double;

  begin
    jfz := DoblDobl_Complex_Jaco_Matrices.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,absfz,denrco,fzrco,fz);
    if absfz > 0.1 then
      if Is_Zero(fzrco)
       then fzloss := -2**30;
       else fzloss := integer32(to_double(log10(fzrco)));
      end if;
    else
      fzloss := -integer32(to_double(log10(denrco)));
    end if;
  end Estimate_Loss_in_Newton_Step;

  procedure Estimate_Loss_in_Newton_Step
              ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
                z : in QuadDobl_Complex_Vectors.Vector;
                jfz : out QuadDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out QuadDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out quad_double;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : quad_double;

  begin
    jfz := QuadDobl_Complex_Jaco_Matrices.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,absfz,denrco,fzrco,fz);
    if absfz > 0.1 then
      if Is_Zero(fzrco)
       then fzloss := -2**30;
       else fzloss := integer32(to_double(log10(fzrco)));
      end if;
    else
      fzloss := -integer32(to_double(log10(fzrco)));
    end if;
  end Estimate_Loss_in_Newton_Step;

  procedure Estimate_Loss_in_Newton_Step
              ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
                jf : in Multprec_Complex_Jaco_Matrices.Jaco_Mat;
                z : in Multprec_Complex_Vectors.Vector;
                jfz : out Multprec_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Multprec_Complex_Vectors.Vector;
                jfzrco,fzrco : out Floating_Number;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : Floating_Number;

  begin
    jfz := Multprec_Complex_Jaco_Matrices.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,absfz,denrco,fzrco,fz);
    if absfz > 0.1 then
      if Equal(fzrco,0.0) then
       fzloss := -2**30;
      else 
        declare
          log10fzrco : Floating_Number := log10(fzrco);
        begin
          fzloss := integer32(Round(log10(fzrco)));
          Clear(log10fzrco);
        end;
      end if;
    else
      declare
        log10denrco : Floating_Number := log10(denrco);
      begin
        fzloss := -integer32(Round(log10(denrco)));
        Clear(log10denrco);
      end;
    end if;
    Clear(absfz); Clear(denrco);
  end Estimate_Loss_in_Newton_Step;

  procedure Estimate_Loss_in_Newton_Step
              ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
                jf : in Standard_Complex_Laur_JacoMats.Jaco_Mat;
                z : in Standard_Complex_Vectors.Vector;
                jfz : out Standard_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Standard_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_float;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : double_float;

  begin
    jfz := Standard_Complex_Laur_JacoMats.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,absfz,denrco,fzrco,fz);
    if absfz > 0.1 then
      if fzrco = 0.0
       then fzloss := -2**30;
       else fzloss := integer32(log10(fzrco));
      end if;
    else
      fzloss := -integer32(log10(denrco));
    end if;
  end Estimate_Loss_in_Newton_Step;

  procedure Estimate_Loss_in_Newton_Step
              ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                jf : in DoblDobl_Complex_Laur_JacoMats.Jaco_Mat;
                z : in DoblDobl_Complex_Vectors.Vector;
                jfz : out DoblDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out DoblDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_double;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : double_double;

  begin
    jfz := DoblDobl_Complex_Laur_JacoMats.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,absfz,denrco,fzrco,fz);
    if absfz > 0.1 then
      if Is_Zero(fzrco)
       then fzloss := -2**30;
       else fzloss := integer32(to_double(log10(fzrco)));
      end if;
    else
      fzloss := -integer32(to_double(log10(denrco)));
    end if;
  end Estimate_Loss_in_Newton_Step;

  procedure Estimate_Loss_in_Newton_Step
              ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Jaco_Mat;
                z : in QuadDobl_Complex_Vectors.Vector;
                jfz : out QuadDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out QuadDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out quad_double;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : quad_double;

  begin
    jfz := QuadDobl_Complex_Laur_JacoMats.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,absfz,denrco,fzrco,fz);
    if absfz > 0.1 then
      if Is_Zero(fzrco)
       then fzloss := -2**30;
       else fzloss := integer32(to_double(log10(fzrco)));
      end if;
    else
      fzloss := -integer32(to_double(log10(fzrco)));
    end if;
  end Estimate_Loss_in_Newton_Step;

  procedure Estimate_Loss_in_Newton_Step
              ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
                jf : in Multprec_Complex_Laur_JacoMats.Jaco_Mat;
                z : in Multprec_Complex_Vectors.Vector;
                jfz : out Multprec_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Multprec_Complex_Vectors.Vector;
                jfzrco,fzrco : out Floating_Number;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : Floating_Number;

  begin
    jfz := Multprec_Complex_Laur_JacoMats.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,absfz,denrco,fzrco,fz);
    if absfz > 0.1 then
      if Equal(fzrco,0.0) then
       fzloss := -2**30;
      else 
        declare
          log10fzrco : Floating_Number := log10(fzrco);
        begin
          fzloss := integer32(Round(log10(fzrco)));
          Clear(log10fzrco);
        end;
      end if;
    else
      declare
        log10denrco : Floating_Number := log10(denrco);
      begin
        fzloss := -integer32(Round(log10(denrco)));
        Clear(log10denrco);
      end;
    end if;
    Clear(absfz); Clear(denrco);
  end Estimate_Loss_in_Newton_Step;

  procedure do_Newton_Step
              ( z : in out Standard_Complex_Vectors.Vector;
                jfz : in Standard_Complex_Matrices.Matrix;
                piv : in Standard_Integer_Vectors.Vector;
                fz : in Standard_Complex_Vectors.Vector;
                err : out double_float ) is

    use Standard_Complex_Numbers;

    dz : Standard_Complex_Vectors.Vector(z'range);
    avl : double_float;

  begin
    for i in dz'range loop
      dz(i) := -fz(i);
    end loop;
    lusolve(jfz,jfz'last(1),piv,dz);
    err := 0.0;
    for i in dz'range loop
      z(i) := z(i) + dz(i);
      avl := AbsVal(dz(i));
      if avl > err
       then err := avl;
      end if;
    end loop;
  end do_Newton_Step;

  procedure do_Newton_Step
              ( z : in out DoblDobl_Complex_Vectors.Vector;
                jfz : in DoblDobl_Complex_Matrices.Matrix;
                piv : in Standard_Integer_Vectors.Vector;
                fz : in DoblDobl_Complex_Vectors.Vector;
                err : out double_double ) is

    use DoblDobl_Complex_Numbers;

    dz : DoblDobl_Complex_Vectors.Vector(z'range);
    avl : double_double;

  begin
    for i in dz'range loop
      dz(i) := -fz(i);
    end loop;
    lusolve(jfz,jfz'last(1),piv,dz);
    err := create(0.0);
    for i in dz'range loop
      z(i) := z(i) + dz(i);
      avl := AbsVal(dz(i));
      if avl > err
       then err := avl;
      end if;
    end loop;
  end do_Newton_Step;

  procedure do_Newton_Step
              ( z : in out QuadDobl_Complex_Vectors.Vector;
                jfz : in QuadDobl_Complex_Matrices.Matrix;
                piv : in Standard_Integer_Vectors.Vector;
                fz : in QuadDobl_Complex_Vectors.Vector;
                err : out quad_double ) is

    use QuadDobl_Complex_Numbers;

    dz : QuadDobl_Complex_Vectors.Vector(z'range);
    avl : quad_double;

  begin
    for i in dz'range loop
      dz(i) := -fz(i);
    end loop;
    lusolve(jfz,jfz'last(1),piv,dz);
    err := create(0.0);
    for i in dz'range loop
      z(i) := z(i) + dz(i);
      avl := AbsVal(dz(i));
      if avl > err
       then err := avl;
      end if;
    end loop;
  end do_Newton_Step;

  procedure do_Newton_Step
              ( z : in out Multprec_Complex_Vectors.Vector;
                jfz : in Multprec_Complex_Matrices.Matrix;
                piv : in Standard_Integer_Vectors.Vector;
                fz : in Multprec_Complex_Vectors.Vector;
                err : out Floating_Number ) is

    use Multprec_Complex_Numbers;

    dz : Multprec_Complex_Vectors.Vector(z'range);
    avl : Floating_Number;

  begin
    for i in dz'range loop
      Copy(fz(i),dz(i));
      Min(dz(i));
    end loop;
    lusolve(jfz,jfz'last(1),piv,dz);
    err := create(0.0);
    for i in dz'range loop
      Add(z(i),dz(i));
      avl := AbsVal(dz(i));
      if avl > err
       then Copy(avl,err);
      end if;
      Clear(avl);
    end loop;
    Multprec_Complex_Vectors.Clear(dz);
  end do_Newton_Step;

  function Minimum ( a,b : integer32 ) return integer32 is
  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Minimum;

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                z : in out Standard_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_float;
                fz : out Standard_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 16;
    jfz : Standard_Complex_Matrices.Matrix(f'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(f,jf,z,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
                z : in out DoblDobl_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_double;
                fz : out DoblDobl_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 32;
    jfz : DoblDobl_Complex_Matrices.Matrix(f'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(f,jf,z,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
                z : in out QuadDobl_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out quad_double;
                fz : out QuadDobl_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 64;
    jfz : QuadDobl_Complex_Matrices.Matrix(f'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(f,jf,z,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
                jf : in Multprec_Complex_Jaco_Matrices.Jaco_Mat;
                z : in out Multprec_Complex_Vectors.Vector;
                prec,want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out Floating_Number;
                fz : out Multprec_Complex_Vectors.Vector;
                fail : out boolean ) is

    jfz : Multprec_Complex_Matrices.Matrix(f'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(f,jf,z,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    if prec + loss < want then
      fail := true;
    else
      fail := false;
      do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
                jf : in Standard_Complex_Laur_JacoMats.Jaco_Mat;
                z : in out Standard_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_float;
                fz : out Standard_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 16;
    jfz : Standard_Complex_Matrices.Matrix(f'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(f,jf,z,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                jf : in DoblDobl_Complex_Laur_JacoMats.Jaco_Mat;
                z : in out DoblDobl_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_double;
                fz : out DoblDobl_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 32;
    jfz : DoblDobl_Complex_Matrices.Matrix(f'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(f,jf,z,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Jaco_Mat;
                z : in out QuadDobl_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out quad_double;
                fz : out QuadDobl_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 64;
    jfz : QuadDobl_Complex_Matrices.Matrix(f'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(f,jf,z,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
                jf : in Multprec_Complex_Laur_JacoMats.Jaco_Mat;
                z : in out Multprec_Complex_Vectors.Vector;
                prec,want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out Floating_Number;
                fz : out Multprec_Complex_Vectors.Vector;
                fail : out boolean ) is

    jfz : Multprec_Complex_Matrices.Matrix(f'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(f,jf,z,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    if prec + loss < want then
      fail := true;
    else
      fail := false;
      do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Standard_Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out double_float; loss : out integer32 ) is
  
    pf : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(f'last))
       := Standard_Complex_Poly_Strings.Parse(natural32(f'last),f);
    cz : constant Standard_Complex_Vectors.Vector(pf'range)
       := Standard_Complex_Vector_Strings.Parse(z);
    jf : Standard_Complex_Jaco_Matrices.Jaco_Mat(pf'range,cz'range)
       := Standard_Complex_Jaco_Matrices.Create(pf);
    fz : Standard_Complex_Vectors.Vector(pf'range);
    jfz : Standard_Complex_Matrices.Matrix(pf'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(pf,jf,cz,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    Standard_Complex_Poly_Systems.Clear(pf);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Standard_Estimate_Loss_for_Polynomial_System;

  procedure DoblDobl_Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out double_double; loss : out integer32 ) is

    pf : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(f'last))
       := DoblDobl_Complex_Poly_Strings.Parse(natural32(f'last),f);
    cz : constant DoblDobl_Complex_Vectors.Vector(pf'range)
       := DoblDobl_Complex_Vector_Strings.Parse(z);
    jf : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(pf'range,cz'range)
       := DoblDobl_Complex_Jaco_Matrices.Create(pf);
    fz : DoblDobl_Complex_Vectors.Vector(pf'range);
    jfz : DoblDobl_Complex_Matrices.Matrix(pf'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(pf,jf,cz,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    DoblDobl_Complex_Poly_Systems.Clear(pf);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end DoblDobl_Estimate_Loss_for_Polynomial_System;

  procedure QuadDobl_Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out quad_double; loss : out integer32 ) is

    pf : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(f'last))
       := QuadDobl_Complex_Poly_Strings.Parse(natural32(f'last),f);
    cz : constant QuadDobl_Complex_Vectors.Vector(pf'range)
       := QuadDobl_Complex_Vector_Strings.Parse(z);
    jf : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(pf'range,cz'range)
       := QuadDobl_Complex_Jaco_Matrices.Create(pf);
    fz : QuadDobl_Complex_Vectors.Vector(pf'range);
    jfz : QuadDobl_Complex_Matrices.Matrix(pf'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(pf,jf,cz,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    QuadDobl_Complex_Poly_Systems.Clear(pf);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end QuadDobl_Estimate_Loss_for_Polynomial_System;

  procedure Multprec_Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                prec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 ) is

    size : constant natural32 := Decimal_to_Size(prec);
    pf : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(f'last))
       := Multprec_Complex_Poly_Strings.Parse(natural32(f'last),size,f);
    cz : Multprec_Complex_Vectors.Vector(pf'range)
       := Multprec_Complex_Vector_Strings.Parse(z);
    jf : Multprec_Complex_Jaco_Matrices.Jaco_Mat(pf'range,cz'range)
       := Multprec_Complex_Jaco_Matrices.Create(pf);
    fz : Multprec_Complex_Vectors.Vector(pf'range);
    jfz : Multprec_Complex_Matrices.Matrix(pf'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(pf,jf,cz,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    Multprec_Complex_Poly_Systems.Clear(pf);
    Multprec_Complex_Jaco_Matrices.Clear(jf);
    Multprec_Complex_Matrices.Clear(jfz);
    Multprec_Complex_Vectors.Clear(cz);
    Multprec_Complex_Vectors.Clear(fz);
  end Multprec_Estimate_Loss_for_Polynomial_System;

  procedure Standard_Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out double_float; loss : out integer32 ) is
  
    pf : Standard_Complex_Laur_Systems.Laur_Sys(1..integer32(f'last))
       := Standard_Complex_Laur_Strings.Parse(natural32(f'last),f);
    cz : constant Standard_Complex_Vectors.Vector(pf'range)
       := Standard_Complex_Vector_Strings.Parse(z);
    jf : Standard_Complex_Laur_JacoMats.Jaco_Mat(pf'range,cz'range)
       := Standard_Complex_Laur_JacoMats.Create(pf);
    fz : Standard_Complex_Vectors.Vector(pf'range);
    jfz : Standard_Complex_Matrices.Matrix(pf'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(pf,jf,cz,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    Standard_Complex_Laur_Systems.Clear(pf);
    Standard_Complex_Laur_JacoMats.Clear(jf);
  end Standard_Estimate_Loss_for_Laurent_Polynomials;

  procedure DoblDobl_Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out double_double; loss : out integer32 ) is

    pf : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..integer32(f'last))
       := DoblDobl_Complex_Laur_Strings.Parse(natural32(f'last),f);
    cz : constant DoblDobl_Complex_Vectors.Vector(pf'range)
       := DoblDobl_Complex_Vector_Strings.Parse(z);
    jf : DoblDobl_Complex_Laur_JacoMats.Jaco_Mat(pf'range,cz'range)
       := DoblDobl_Complex_Laur_JacoMats.Create(pf);
    fz : DoblDobl_Complex_Vectors.Vector(pf'range);
    jfz : DoblDobl_Complex_Matrices.Matrix(pf'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(pf,jf,cz,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    DoblDobl_Complex_Laur_Systems.Clear(pf);
    DoblDobl_Complex_Laur_JacoMats.Clear(jf);
  end DoblDobl_Estimate_Loss_for_Laurent_Polynomials;

  procedure QuadDobl_Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out quad_double; loss : out integer32 ) is

    pf : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..integer32(f'last))
       := QuadDobl_Complex_Laur_Strings.Parse(natural32(f'last),f);
    cz : constant QuadDobl_Complex_Vectors.Vector(pf'range)
       := QuadDobl_Complex_Vector_Strings.Parse(z);
    jf : QuadDobl_Complex_Laur_JacoMats.Jaco_Mat(pf'range,cz'range)
       := QuadDobl_Complex_Laur_JacoMats.Create(pf);
    fz : QuadDobl_Complex_Vectors.Vector(pf'range);
    jfz : QuadDobl_Complex_Matrices.Matrix(pf'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(pf,jf,cz,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    QuadDobl_Complex_Laur_Systems.Clear(pf);
    QuadDobl_Complex_Laur_JacoMats.Clear(jf);
  end QuadDobl_Estimate_Loss_for_Laurent_Polynomials;

  procedure Multprec_Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                prec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 ) is

    size : constant natural32 := Decimal_to_Size(prec);
    pf : Multprec_Complex_Laur_Systems.Laur_Sys(1..integer32(f'last))
       := Multprec_Complex_Laur_Strings.Parse(natural32(f'last),size,f);
    cz : Multprec_Complex_Vectors.Vector(pf'range)
       := Multprec_Complex_Vector_Strings.Parse(z);
    jf : Multprec_Complex_Laur_JacoMats.Jaco_Mat(pf'range,cz'range)
       := Multprec_Complex_Laur_JacoMats.Create(pf);
    fz : Multprec_Complex_Vectors.Vector(pf'range);
    jfz : Multprec_Complex_Matrices.Matrix(pf'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(pf,jf,cz,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Minimum(jflss,fzlss);
    Multprec_Complex_Laur_Systems.Clear(pf);
    Multprec_Complex_Laur_JacoMats.Clear(jf);
    Multprec_Complex_Matrices.Clear(jfz);
    Multprec_Complex_Vectors.Clear(cz);
    Multprec_Complex_Vectors.Clear(fz);
  end Multprec_Estimate_Loss_for_Laurent_Polynomials;

  procedure Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                maxprec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 ) is

    dim : constant natural32 := natural32(f'last);
    d_jfrco,d_fzrco : double_float;
    dd_jfrco,dd_fzrco : double_double;
    qd_jfrco,qd_fzrco : quad_double;
    prcn : natural32;

  begin
    if Symbol_Table.Number < dim
     then Symbol_Table.Init(dim);
    end if;
    Standard_Estimate_Loss_for_Polynomial_System(f,z,d_jfrco,d_fzrco,loss);
    if loss > -15 then
      jfrco := create(d_jfrco);
      fzrco := create(d_fzrco);
    else
      DoblDobl_Estimate_Loss_for_Polynomial_System(f,z,dd_jfrco,dd_fzrco,loss);
      if loss > -30 then
        jfrco := Multprec_DoblDobl_Convertors.to_floating_number(dd_jfrco);
        fzrco := Multprec_DoblDobl_Convertors.to_floating_number(dd_fzrco);
      else
        QuadDobl_Estimate_Loss_for_Polynomial_System
          (f,z,qd_jfrco,qd_fzrco,loss);
        if loss > -60 then
          jfrco := Multprec_QuadDobl_Convertors.to_floating_number(qd_jfrco);
          fzrco := Multprec_QuadDobl_Convertors.to_floating_number(qd_fzrco);
        else
          prcn := 80;
          loop
            Multprec_Estimate_Loss_for_Polynomial_System
              (f,z,prcn,jfrco,fzrco,loss);
            exit when (natural32(loss) < prcn);
            Clear(jfrco); Clear(fzrco);
            prcn := prcn + prcn/16;  -- increment precision with 25%
            exit when (prcn > maxprec);
          end loop;
        end if;
      end if;
    end if;
  end Estimate_Loss_for_Polynomial_System;

  function Estimate_Loss_for_Polynomial_System
              ( f : Array_of_Strings; z : string; maxprec : natural32)
              return integer32 is

    res : integer32;
    jfrco,fzrco : Floating_Number;

  begin
    Estimate_Loss_for_Polynomial_System(f,z,maxprec,jfrco,fzrco,res);
    Clear(jfrco); Clear(fzrco);
    return res;
  end Estimate_Loss_for_Polynomial_System;

  procedure Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                maxprec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 ) is

    dim : constant natural32 := natural32(f'last);
    d_jfrco,d_fzrco : double_float;
    dd_jfrco,dd_fzrco : double_double;
    qd_jfrco,qd_fzrco : quad_double;
    prcn : natural32;

  begin
    if Symbol_Table.Number < dim
     then Symbol_Table.Init(dim);
    end if;
    Standard_Estimate_Loss_for_Laurent_Polynomials(f,z,d_jfrco,d_fzrco,loss);
    if loss > -15 then
      jfrco := create(d_jfrco);
      fzrco := create(d_fzrco);
    else
      DoblDobl_Estimate_Loss_for_Laurent_Polynomials
        (f,z,dd_jfrco,dd_fzrco,loss);
      if loss > -30 then
        jfrco := Multprec_DoblDobl_Convertors.to_floating_number(dd_jfrco);
        fzrco := Multprec_DoblDobl_Convertors.to_floating_number(dd_fzrco);
      else
        QuadDobl_Estimate_Loss_for_Laurent_Polynomials
          (f,z,qd_jfrco,qd_fzrco,loss);
        if loss > -60 then
          jfrco := Multprec_QuadDobl_Convertors.to_floating_number(qd_jfrco);
          fzrco := Multprec_QuadDobl_Convertors.to_floating_number(qd_fzrco);
        else
          prcn := 80;
          loop
            Multprec_Estimate_Loss_for_Laurent_Polynomials
              (f,z,prcn,jfrco,fzrco,loss);
            exit when (natural32(loss) < prcn);
            Clear(jfrco); Clear(fzrco);
            prcn := prcn + prcn/16;  -- increment precision with 25%
            exit when (prcn > maxprec);
          end loop;
        end if;
      end if;
    end if;
  end Estimate_Loss_for_Laurent_Polynomials;

  function Estimate_Loss_for_Laurent_Polynomials
              ( f : Array_of_Strings; z : string; maxprec : natural32)
              return integer32 is

    res : integer32;
    jfrco,fzrco : Floating_Number;

  begin
    Estimate_Loss_for_Laurent_Polynomials(f,z,maxprec,jfrco,fzrco,res);
    Clear(jfrco); Clear(fzrco);
    return res;
  end Estimate_Loss_for_Laurent_Polynomials;

  procedure Standard_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float ) is

    n : constant integer32 := integer32(f'last);
    p : Standard_Complex_Poly_Systems.Poly_Sys(1..n)
      := Standard_Complex_Poly_Strings.Parse(natural32(n),f);
    x : Standard_Complex_Vectors.Vector(p'range)
      := Standard_Complex_Vector_Strings.Parse(z.all);
    jf : Standard_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Complex_Jaco_Matrices.Create(p);
    fz : constant Standard_Complex_Vectors.Vector(p'range)
       := Standard_Complex_Poly_SysFun.Eval(p,x);
    jfz : Standard_Complex_Matrices.Matrix(p'range,x'range)
        := Standard_Complex_Jaco_Matrices.Eval(jf,x);
    piv : Standard_Integer_Vectors.Vector(x'range);

  begin
    lufco(jfz,n,piv,rco);
    do_Newton_Step(x,jfz,piv,fz,err);
    res := Max_Norm(fz);
    declare
      nz : constant string
         := Standard_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Standard_Newton_Step_on_Polynomial_System;

  procedure DoblDobl_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float ) is

    n : constant integer32 := integer32(f'last);
    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..n)
      := DoblDobl_Complex_Poly_Strings.Parse(natural32(n),f);
    x : DoblDobl_Complex_Vectors.Vector(p'range)
      := DoblDobl_Complex_Vector_Strings.Parse(z.all);
    jf : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Complex_Jaco_Matrices.Create(p);
    fz : constant DoblDobl_Complex_Vectors.Vector(p'range)
       := DoblDobl_Complex_Poly_SysFun.Eval(p,x);
    jfz : DoblDobl_Complex_Matrices.Matrix(p'range,x'range)
        := DoblDobl_Complex_Jaco_Matrices.Eval(jf,x);
    piv : Standard_Integer_Vectors.Vector(x'range);
    dd_rco,dd_err : double_double;

  begin
    lufco(jfz,n,piv,dd_rco);
    rco := to_double(dd_rco);
    do_Newton_Step(x,jfz,piv,fz,dd_err);
    err := to_double(dd_err);
    res := to_double(Max_Norm(fz));
    declare
      nz : constant string
         := DoblDobl_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    DoblDobl_Complex_Poly_Systems.Clear(p);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end DoblDobl_Newton_Step_on_Polynomial_System;

  procedure QuadDobl_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float ) is

    n : constant integer32 := integer32(f'last);
    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..n)
      := QuadDobl_Complex_Poly_Strings.Parse(natural32(n),f);
    x : QuadDobl_Complex_Vectors.Vector(p'range)
      := QuadDobl_Complex_Vector_Strings.Parse(z.all);
    jf : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Complex_Jaco_Matrices.Create(p);
    fz : constant QuadDobl_Complex_Vectors.Vector(p'range)
       := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    jfz : QuadDobl_Complex_Matrices.Matrix(p'range,x'range)
        := QuadDobl_Complex_Jaco_Matrices.Eval(jf,x);
    piv : Standard_Integer_Vectors.Vector(x'range);
    qd_rco,qd_err : quad_double;

  begin
    lufco(jfz,n,piv,qd_rco);
    rco := to_double(qd_rco);
    do_Newton_Step(x,jfz,piv,fz,qd_err);
    err := to_double(qd_err);
    res := to_double(Max_Norm(fz));
    declare
      nz : constant string
         := QuadDobl_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    QuadDobl_Complex_Poly_Systems.Clear(p);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end QuadDobl_Newton_Step_on_Polynomial_System;

  procedure Multprec_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                prcn : in natural32; err,rco,res : out double_float ) is

    size : constant natural32 := Decimal_to_Size(prcn);
    n : constant integer32 := integer32(f'last);
    p : Multprec_Complex_Poly_Systems.Poly_Sys(1..n)
       := Multprec_Complex_Poly_Strings.Parse(natural32(n),size,f);
    x : Multprec_Complex_Vectors.Vector(p'range)
      := Multprec_Complex_Vector_Strings.Parse(z.all);
    jf : Multprec_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Multprec_Complex_Jaco_Matrices.Create(p);
    fz : Multprec_Complex_Vectors.Vector(p'range);
    jfz : Multprec_Complex_Matrices.Matrix(p'range,x'range);
    piv : Standard_Integer_Vectors.Vector(x'range);
    mp_rco,mp_err,mp_res : Floating_Number;

  begin
    Multprec_Complex_Vector_Tools.Set_Size(x,size);
    fz := Multprec_Complex_Poly_SysFun.Eval(p,x);
    jfz := Multprec_Complex_Jaco_Matrices.Eval(jf,x);
    lufco(jfz,n,piv,mp_rco);
    rco := Round(mp_rco);
    do_Newton_Step(x,jfz,piv,fz,mp_err);
    err := Round(mp_err);
    mp_res := Max_Norm(fz);
    res := Round(mp_res);
    declare
      nz : constant string
         := Multprec_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    Multprec_Complex_Poly_Systems.Clear(p);
    Multprec_Complex_Jaco_Matrices.Clear(jf);
    Multprec_Complex_Matrices.Clear(jfz);
    Multprec_Complex_Vectors.Clear(fz);
    Multprec_Complex_Vectors.Clear(x);
    Clear(mp_rco); Clear(mp_err); Clear(mp_res);
  end Multprec_Newton_Step_on_Polynomial_System;

  procedure Standard_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float ) is

    n : constant integer32 := integer32(f'last);
    p : Standard_Complex_Laur_Systems.Laur_Sys(1..n)
      := Standard_Complex_Laur_Strings.Parse(natural32(n),f);
    x : Standard_Complex_Vectors.Vector(p'range)
      := Standard_Complex_Vector_Strings.Parse(z.all);
    jf : Standard_Complex_Laur_JacoMats.Jaco_Mat(p'range,x'range)
       := Standard_Complex_Laur_JacoMats.Create(p);
    fz : constant Standard_Complex_Vectors.Vector(p'range)
       := Standard_Complex_Laur_SysFun.Eval(p,x);
    jfz : Standard_Complex_Matrices.Matrix(p'range,x'range)
        := Standard_Complex_Laur_JacoMats.Eval(jf,x);
    piv : Standard_Integer_Vectors.Vector(x'range);

  begin
    lufco(jfz,n,piv,rco);
    do_Newton_Step(x,jfz,piv,fz,err);
    res := Max_Norm(fz);
    declare
      nz : constant string
         := Standard_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    Standard_Complex_Laur_Systems.Clear(p);
    Standard_Complex_Laur_JacoMats.Clear(jf);
 -- exception
 --   when others =>
 --     put_line("Exception in Standard_Newton_Step_on_Laurent_Polynomials.");
 --     put_line("The string z :"); put_line(z.all);
 --     put_line("The vector x :"); put_line(x);
 --     raise;
  end Standard_Newton_Step_on_Laurent_Polynomials;

  procedure DoblDobl_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float ) is

    n : constant integer32 := integer32(f'last);
    p : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..n)
      := DoblDobl_Complex_Laur_Strings.Parse(natural32(n),f);
    x : DoblDobl_Complex_Vectors.Vector(p'range)
      := DoblDobl_Complex_Vector_Strings.Parse(z.all);
    jf : DoblDobl_Complex_Laur_JacoMats.Jaco_Mat(p'range,x'range)
       := DoblDobl_Complex_Laur_JacoMats.Create(p);
    fz : constant DoblDobl_Complex_Vectors.Vector(p'range)
       := DoblDobl_Complex_Laur_SysFun.Eval(p,x);
    jfz : DoblDobl_Complex_Matrices.Matrix(p'range,x'range)
        := DoblDobl_Complex_Laur_JacoMats.Eval(jf,x);
    piv : Standard_Integer_Vectors.Vector(x'range);
    dd_rco,dd_err : double_double;

  begin
    lufco(jfz,n,piv,dd_rco);
    rco := to_double(dd_rco);
    do_Newton_Step(x,jfz,piv,fz,dd_err);
    err := to_double(dd_err);
    res := to_double(Max_Norm(fz));
    declare
      nz : constant string
         := DoblDobl_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    DoblDobl_Complex_Laur_Systems.Clear(p);
    DoblDobl_Complex_Laur_JacoMats.Clear(jf);
  end DoblDobl_Newton_Step_on_Laurent_Polynomials;

  procedure QuadDobl_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float ) is

    n : constant integer32 := integer32(f'last);
    p : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..n)
      := QuadDobl_Complex_Laur_Strings.Parse(natural32(n),f);
    x : QuadDobl_Complex_Vectors.Vector(p'range)
      := QuadDobl_Complex_Vector_Strings.Parse(z.all);
    jf : QuadDobl_Complex_Laur_JacoMats.Jaco_Mat(p'range,x'range)
       := QuadDobl_Complex_Laur_JacoMats.Create(p);
    fz : constant QuadDobl_Complex_Vectors.Vector(p'range)
       := QuadDobl_Complex_Laur_SysFun.Eval(p,x);
    jfz : QuadDobl_Complex_Matrices.Matrix(p'range,x'range)
        := QuadDobl_Complex_Laur_JacoMats.Eval(jf,x);
    piv : Standard_Integer_Vectors.Vector(x'range);
    qd_rco,qd_err : quad_double;

  begin
    lufco(jfz,n,piv,qd_rco);
    rco := to_double(qd_rco);
    do_Newton_Step(x,jfz,piv,fz,qd_err);
    err := to_double(qd_err);
    res := to_double(Max_Norm(fz));
    declare
      nz : constant string
         := QuadDobl_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    QuadDobl_Complex_Laur_Systems.Clear(p);
    QuadDobl_Complex_Laur_JacoMats.Clear(jf);
  end QuadDobl_Newton_Step_on_Laurent_Polynomials;

  procedure Multprec_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                prcn : in natural32; err,rco,res : out double_float ) is

    size : constant natural32 := Decimal_to_Size(prcn);
    n : constant integer32 := integer32(f'last);
    p : Multprec_Complex_Laur_Systems.Laur_Sys(1..n)
       := Multprec_Complex_Laur_Strings.Parse(natural32(n),size,f);
    x : Multprec_Complex_Vectors.Vector(p'range)
      := Multprec_Complex_Vector_Strings.Parse(z.all);
    jf : Multprec_Complex_Laur_JacoMats.Jaco_Mat(p'range,x'range)
       := Multprec_Complex_Laur_JacoMats.Create(p);
    fz : Multprec_Complex_Vectors.Vector(p'range);
    jfz : Multprec_Complex_Matrices.Matrix(p'range,x'range);
    piv : Standard_Integer_Vectors.Vector(x'range);
    mp_rco,mp_err,mp_res : Floating_Number;

  begin
    Multprec_Complex_Vector_Tools.Set_Size(x,size);
    fz := Multprec_Complex_Laur_SysFun.Eval(p,x);
    jfz := Multprec_Complex_Laur_JacoMats.Eval(jf,x);
    lufco(jfz,n,piv,mp_rco);
    rco := Round(mp_rco);
    do_Newton_Step(x,jfz,piv,fz,mp_err);
    err := Round(mp_err);
    mp_res := Max_Norm(fz);
    res := Round(mp_res);
    declare
      nz : constant string
         := Multprec_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    Multprec_Complex_Laur_Systems.Clear(p);
    Multprec_Complex_Laur_JacoMats.Clear(jf);
    Multprec_Complex_Matrices.Clear(jfz);
    Multprec_Complex_Vectors.Clear(fz);
    Multprec_Complex_Vectors.Clear(x);
    Clear(mp_rco); Clear(mp_err); Clear(mp_res);
  end Multprec_Newton_Step_on_Laurent_Polynomials;

  procedure do_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                loss,want : in integer32; err,rco,res : out double_float ) is

    precision : constant natural32
              := natural32(-loss) + natural32(want);

  begin
    if precision <= 16 then
      Standard_Newton_Step_on_Polynomial_System(f,z,err,rco,res);
    elsif precision <= 32 then
      DoblDobl_Newton_Step_on_Polynomial_System(f,z,err,rco,res);
    elsif precision <= 64 then
      QuadDobl_Newton_Step_on_Polynomial_System(f,z,err,rco,res);
    else
      Multprec_Newton_Step_on_Polynomial_System(f,z,precision,err,rco,res);
    end if;
  end do_Newton_Step_on_Polynomial_System;

  procedure do_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                loss,want : in integer32; err,rco,res : out double_float ) is

    precision : constant natural32
              := natural32(-loss) + natural32(want);

  begin
    if precision <= 16 then
      Standard_Newton_Step_on_Laurent_Polynomials(f,z,err,rco,res);
    elsif precision <= 32 then
      DoblDobl_Newton_Step_on_Laurent_Polynomials(f,z,err,rco,res);
    elsif precision <= 64 then
      QuadDobl_Newton_Step_on_Laurent_Polynomials(f,z,err,rco,res);
    else
      Multprec_Newton_Step_on_Laurent_Polynomials(f,z,precision,err,rco,res);
    end if;
 -- exception
 --   when others =>
 --     put_line("Exception in do_Newton_Step_on_Laurent_Polynomials ...");
 --     raise;
  end do_Newton_Step_on_Laurent_Polynomials;

  procedure Newton_Steps_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                want : in integer32; maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float ) is

    precision : natural32;
    err_accu,res_accu : integer32;

  begin
    for i in 1..maxitr loop
      loss := Estimate_Loss_for_Polynomial_System(f,z.all,maxprc);
      do_Newton_Step_on_Polynomial_System(f,z,loss,want,err,rco,res);
      precision := natural32(-loss) + natural32(want);
      if err = 0.0
       then err_accu := integer32(precision);
       else err_accu := abs(integer32(log10(err)));
      end if;
      if res = 0.0
       then res_accu := integer32(precision);
       else res_accu := abs(integer32(log10(res)));
      end if;
      exit when ((err_accu >= want) and (res_accu >= want));
    end loop;
  end Newton_Steps_on_Polynomial_System;

  procedure Newton_Steps_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                want : in integer32; maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float ) is

    precision : natural32;
    err_accu,res_accu : integer32;

  begin
    for i in 1..maxitr loop
      loss := Estimate_Loss_for_Laurent_Polynomials(f,z.all,maxprc);
      do_Newton_Step_on_Laurent_Polynomials(f,z,loss,want,err,rco,res);
      precision := natural32(-loss) + natural32(want);
      if err = 0.0 
       then err_accu := integer32(precision);
       else err_accu := abs(integer32(log10(err)));
      end if;
      if res = 0.0
       then res_accu := integer32(precision);
       else res_accu := abs(integer32(log10(res)));
      end if;
      exit when ((err_accu >= want) and (res_accu >= want));
    end loop;
  end Newton_Steps_on_Laurent_Polynomials;

  procedure Newton_Steps_on_Polynomial_System
              ( file : in file_type;
                f : in Array_of_Strings; z : in out Link_to_String;
                want : in integer32; maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float ) is

    precision : natural32;
    err_accu,res_accu : integer32;

  begin
    for i in 1..maxitr loop
      put(file,"estimated loss in step "); put(file,i,1); put(file," : ");
      loss := Estimate_Loss_for_Polynomial_System(f,z.all,maxprc);
      put(file,loss,1);
      precision := natural32(-loss) + natural32(want);
      put(file,", precision : "); put(file,precision,1); new_line(file);
      do_Newton_Step_on_Polynomial_System(f,z,loss,want,err,rco,res);
      put(file,"  err :"); put(file,err,3);
      put(file,"  rco :"); put(file,rco,3);
      put(file,"  res :"); put(file,res,3); new_line(file);
      if err = 0.0
       then err_accu := integer32(precision);
       else err_accu := abs(integer32(log10(err)));
      end if;
      if res = 0.0
       then res_accu := integer32(precision);
       else res_accu := abs(integer32(log10(res)));
      end if;
      exit when ((err_accu >= want) and (res_accu >= want));
    end loop;
  end Newton_Steps_on_Polynomial_System;

  procedure Newton_Steps_on_Laurent_Polynomials
              ( file : in file_type;
                f : in Array_of_Strings; z : in out Link_to_String;
                want : in integer32; maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float ) is

    precision : natural32;
    err_accu,res_accu : integer32;

  begin
    for i in 1..maxitr loop
      put(file,"estimated loss in step "); put(file,i,1); put(file," : ");
      loss := Estimate_Loss_for_Laurent_Polynomials(f,z.all,maxprc);
      put(file,loss,1);
      precision := natural32(-loss) + natural32(want);
      put(file,", precision : "); put(file,precision,1); new_line(file);
      do_Newton_Step_on_Laurent_Polynomials(f,z,loss,want,err,rco,res);
      put(file,"  err :"); put(file,err,3);
      put(file,"  rco :"); put(file,rco,3);
      put(file,"  res :"); put(file,res,3); new_line(file);
      if err = 0.0
       then err_accu := integer32(precision);
       else err_accu := abs(integer32(log10(err)));
      end if;
      if res = 0.0
       then res_accu := integer32(precision);
       else res_accu := abs(integer32(log10(res)));
      end if;
      exit when ((err_accu >= want) and (res_accu >= want));
    end loop;
 -- exception
 --   when others =>
 --     put_line("Exception in Newton_Steps_on_Laurent_Polynomials");
 --     raise;
  end Newton_Steps_on_Laurent_Polynomials;
                
end Varbprec_Complex_Newton_Steps;
