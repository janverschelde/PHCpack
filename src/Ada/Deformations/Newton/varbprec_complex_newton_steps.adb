with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Multprec_DoblDobl_Convertors;
with Multprec_QuadDobl_Convertors;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;
with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Standard_Complex_Vector_Strings;
with DoblDobl_Complex_Vector_Strings;
with QuadDobl_Complex_Vector_Strings;
with Multprec_Complex_Vector_Strings;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;
with Symbol_Table;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Strings;
with DoblDobl_Complex_Poly_Functions;    use DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Poly_Functions;    use QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Strings;
with Multprec_Complex_Poly_Functions;    use Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_Strings;
with Varbprec_Complex_Linear_Solvers;    use Varbprec_Complex_Linear_Solvers;
with Varbprec_Polynomial_Evaluations;    use Varbprec_Polynomial_Evaluations;

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
  begin
    jfz := Standard_Complex_Jaco_Matrices.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,fzrco,fz);
    if fzrco = 0.0
     then fzloss := -2**30;
     else fzloss := integer32(log10(fzrco));
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

  begin
    jfz := DoblDobl_Complex_Jaco_Matrices.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,fzrco,fz);
    if Is_Zero(fzrco)
     then fzloss := -2**30;
     else fzloss := integer32(to_double(log10(fzrco)));
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

  begin
    jfz := QuadDobl_Complex_Jaco_Matrices.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,fzrco,fz);
    if Is_Zero(fzrco)
     then fzloss := -2**30;
     else fzloss := integer32(to_double(log10(fzrco)));
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
  begin
    jfz := Multprec_Complex_Jaco_Matrices.Eval(jf,z);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    Evaluate_with_Inverse_Condition(f,z,fzrco,fz);
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

  procedure Standard_Estimate_Loss_of_Accuracy
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
  end Standard_Estimate_Loss_of_Accuracy;

  procedure DoblDobl_Estimate_Loss_of_Accuracy
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
  end DoblDobl_Estimate_Loss_of_Accuracy;

  procedure QuadDobl_Estimate_Loss_of_Accuracy
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
  end QuadDobl_Estimate_Loss_of_Accuracy;

  procedure Multprec_Estimate_Loss_of_Accuracy
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
  end Multprec_Estimate_Loss_of_Accuracy;

  procedure Estimate_Loss_of_Accuracy
              ( f : in Array_of_Strings; z : in string;
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
    Standard_Estimate_Loss_of_Accuracy(f,z,d_jfrco,d_fzrco,loss);
    if loss > -15 then
      jfrco := create(d_jfrco);
      fzrco := create(d_fzrco);
    else
      DoblDobl_Estimate_Loss_of_Accuracy(f,z,dd_jfrco,dd_fzrco,loss);
      if loss > -30 then
        jfrco := Multprec_DoblDobl_Convertors.to_floating_number(dd_jfrco);
        fzrco := Multprec_DoblDobl_Convertors.to_floating_number(dd_fzrco);
      else
        QuadDobl_Estimate_Loss_of_Accuracy(f,z,qd_jfrco,qd_fzrco,loss);
        if loss > -60 then
          jfrco := Multprec_QuadDobl_Convertors.to_floating_number(qd_jfrco);
          fzrco := Multprec_QuadDobl_Convertors.to_floating_number(qd_fzrco);
        else
          prcn := 80;
          for i in 1..16 loop  -- no infinite loop ...
            Multprec_Estimate_Loss_of_Accuracy(f,z,prcn,jfrco,fzrco,loss);
            exit when (natural32(loss) < prcn);
            Clear(jfrco); Clear(fzrco);
            prcn := prcn + 16; -- increment working precision
          end loop;
        end if;
      end if;
    end if;
  end Estimate_Loss_of_Accuracy;

  function Estimate_Loss_of_Accuracy
              ( f : Array_of_Strings; z : string ) return integer32 is

    res : integer32;
    jfrco,fzrco : Floating_Number;

  begin
    Estimate_Loss_of_Accuracy(f,z,jfrco,fzrco,res);
    Clear(jfrco); Clear(fzrco);
    return res;
  end Estimate_Loss_of_Accuracy;
                
end Varbprec_Complex_Newton_Steps;
