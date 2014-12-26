with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;
with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Multprec_DoblDobl_Convertors;
with Multprec_QuadDobl_Convertors;
with Standard_Complex_Vectors;
with Standard_Complex_Vector_Strings;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Strings;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Strings;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vector_Strings;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Varbprec_Complex_Linear_Solvers;    use Varbprec_Complex_Linear_Solvers;
with Varbprec_Polynomial_Evaluations;    use Varbprec_Polynomial_Evaluations;
with Varbprec_Complex_Newton_Steps;
with Varbprec_Homotopy;

package body Varbprec_Corrector_Steps is

-- PART I : estimators and corrector steps at various level of precision

  procedure Estimate_Loss_in_Newton_Step
              ( z : in Standard_Complex_Vectors.Vector;
                t : in Standard_Complex_Numbers.Complex_Number;
                jfz : out Standard_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Standard_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_float;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : double_float;
    h : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.Standard_Homotopy_System;
    zt : Standard_Complex_Vectors.Vector(z'first..z'last+1);

  begin
    jfz := Varbprec_Homotopy.Diff(z,t);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    zt(z'range) := z;
    zt(zt'last) := t;
    Evaluate_with_Inverse_Condition(h.all,zt,absfz,denrco,fzrco,fz);
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
              ( z : in DoblDobl_Complex_Vectors.Vector;
                t : in DoblDobl_Complex_Numbers.Complex_Number;
                jfz : out DoblDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out DoblDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_double;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : double_double;
    h : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.DoblDobl_Homotopy_System;
    zt : DoblDobl_Complex_Vectors.Vector(z'first..z'last+1);

  begin
    jfz := Varbprec_Homotopy.Diff(z,t);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    zt(z'range) := z;
    zt(zt'last) := t;
    Evaluate_with_Inverse_Condition(h.all,zt,absfz,denrco,fzrco,fz);
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
              ( z : in QuadDobl_Complex_Vectors.Vector;
                t : in QuadDobl_Complex_Numbers.Complex_Number;
                jfz : out QuadDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out QuadDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out quad_double;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : quad_double;
    h : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.QuadDobl_Homotopy_System;
    zt : QuadDobl_Complex_Vectors.Vector(z'first..z'last+1);

  begin
    jfz := Varbprec_Homotopy.Diff(z,t);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    zt(z'range) := z;
    zt(zt'last) := t;
    Evaluate_with_Inverse_Condition(h.all,zt,absfz,denrco,fzrco,fz);
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
              ( z : in Multprec_Complex_Vectors.Vector;
                t : in Multprec_Complex_Numbers.Complex_Number;
                deci : in natural32;
                jfz : out Multprec_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Multprec_Complex_Vectors.Vector;
                jfzrco,fzrco : out Floating_Number;
                jfzloss,fzloss : out integer32 ) is

    absfz,denrco : Floating_Number;
    h : constant Multprec_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.Multprec_Homotopy_System(deci);
    zt : Multprec_Complex_Vectors.Vector(z'first..z'last+1);

  begin
    jfz := Varbprec_Homotopy.Diff(z,t,deci);
    Estimated_Loss_of_Decimal_Places(jfz,piv,jfzrco,jfzloss);
    for i in z'range loop
      Multprec_Complex_Numbers.Copy(z(i),zt(i));
      Multprec_Complex_Numbers.Copy(t,zt(zt'last));
    end loop;
    Evaluate_with_Inverse_Condition(h.all,zt,absfz,denrco,fzrco,fz);
    Multprec_Complex_Vectors.Clear(zt);
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

  procedure Newton_Step_to_Wanted_Accuracy
              ( z : in out Standard_Complex_Vectors.Vector;
                t : in Standard_Complex_Numbers.Complex_Number;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_float;
                fz : out Standard_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 16;
    jfz : Standard_Complex_Matrices.Matrix(z'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(z,t,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Varbprec_Complex_Newton_Steps.Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      Varbprec_Complex_Newton_Steps.do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( z : in out DoblDobl_Complex_Vectors.Vector;
                t : in DoblDobl_Complex_Numbers.Complex_Number;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_double;
                fz : out DoblDobl_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 32;
    jfz : DoblDobl_Complex_Matrices.Matrix(z'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(z,t,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Varbprec_Complex_Newton_Steps.Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      Varbprec_Complex_Newton_Steps.do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( z : in out QuadDobl_Complex_Vectors.Vector;
                t : in QuadDobl_Complex_Numbers.Complex_Number;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out quad_double;
                fz : out QuadDobl_Complex_Vectors.Vector;
                fail : out boolean ) is

    precision : constant integer32 := 64;
    jfz : QuadDobl_Complex_Matrices.Matrix(z'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(z,t,jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Varbprec_Complex_Newton_Steps.Minimum(jflss,fzlss);
    if precision + loss < want then
      fail := true;
    else
      fail := false;
      Varbprec_Complex_Newton_Steps.do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

  procedure Newton_Step_to_Wanted_Accuracy
              ( z : in out Multprec_Complex_Vectors.Vector;
                t : in Multprec_Complex_Numbers.Complex_Number;
                prec,want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out Floating_Number;
                fz : out Multprec_Complex_Vectors.Vector;
                fail : out boolean ) is

    jfz : Multprec_Complex_Matrices.Matrix(z'range,z'range);
    piv : Standard_Integer_Vectors.Vector(z'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step
      (z,t,natural32(prec),jfz,piv,fz,jfzrco,fzrco,jflss,fzlss);
    loss := Varbprec_Complex_Newton_Steps.Minimum(jflss,fzlss);
    if prec + loss < want then
      fail := true;
    else
      fail := false;
      Varbprec_Complex_Newton_Steps.do_Newton_Step(z,jfz,piv,fz,err);
    end if;
  end Newton_Step_to_Wanted_Accuracy;

-- PART II : estimating loss of accuracy

  procedure Standard_Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in Standard_Complex_Numbers.Complex_Number;
                jfrco,fzrco : out double_float; loss : out integer32 ) is

    h : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.Standard_Homotopy_System;
    cz : constant Standard_Complex_Vectors.Vector(h'range)
       := Standard_Complex_Vector_Strings.Parse(z);
    fz : Standard_Complex_Vectors.Vector(h'range);
    jfz : Standard_Complex_Matrices.Matrix(h'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(cz,t,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Varbprec_Complex_Newton_Steps.Minimum(jflss,fzlss);
  end Standard_Estimate_Loss_for_Polynomial_Homotopy;

  procedure DoblDobl_Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in DoblDobl_Complex_Numbers.Complex_Number;
                jfrco,fzrco : out double_double; loss : out integer32 ) is

    h : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.DoblDobl_Homotopy_System;
    cz : constant DoblDobl_Complex_Vectors.Vector(h'range)
       := DoblDobl_Complex_Vector_Strings.Parse(z);
    fz : DoblDobl_Complex_Vectors.Vector(h'range);
    jfz : DoblDobl_Complex_Matrices.Matrix(h'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(cz,t,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Varbprec_Complex_Newton_Steps.Minimum(jflss,fzlss);
  end DoblDobl_Estimate_Loss_for_Polynomial_Homotopy;

  procedure QuadDobl_Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in QuadDobl_Complex_Numbers.Complex_Number;
                jfrco,fzrco : out quad_double; loss : out integer32 ) is

    h : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.QuadDobl_Homotopy_System;
    cz : constant QuadDobl_Complex_Vectors.Vector(h'range)
       := QuadDobl_Complex_Vector_Strings.Parse(z);
    fz : QuadDobl_Complex_Vectors.Vector(h'range);
    jfz : QuadDobl_Complex_Matrices.Matrix(h'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(cz,t,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Varbprec_Complex_Newton_Steps.Minimum(jflss,fzlss);
  end QuadDobl_Estimate_Loss_for_Polynomial_Homotopy;

  procedure Multprec_Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in Multprec_Complex_Numbers.Complex_Number;
                prec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 ) is

    h : constant Multprec_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.Multprec_Homotopy_System(prec);
    cz : Multprec_Complex_Vectors.Vector(h'range)
       := Multprec_Complex_Vector_Strings.Parse(z);
    fz : Multprec_Complex_Vectors.Vector(h'range);
    jfz : Multprec_Complex_Matrices.Matrix(h'range,cz'range);
    piv : Standard_Integer_Vectors.Vector(cz'range);
    jflss,fzlss : integer32;

  begin
    Estimate_Loss_in_Newton_Step(cz,t,prec,jfz,piv,fz,jfrco,fzrco,jflss,fzlss);
    loss := Varbprec_Complex_Newton_Steps.Minimum(jflss,fzlss);
    Multprec_Complex_Matrices.Clear(jfz);
    Multprec_Complex_Vectors.Clear(cz);
    Multprec_Complex_Vectors.Clear(fz);
  end Multprec_Estimate_Loss_for_Polynomial_Homotopy;

  procedure Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in Standard_Complex_Numbers.Complex_Number;
                maxprec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 ) is

    h : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.Standard_Homotopy_System;
    dim : constant natural32 := natural32(h'last);
    d_jfrco,d_fzrco : double_float;
    dd_jfrco,dd_fzrco,dd_re_t,dd_im_t : double_double;
    qd_jfrco,qd_fzrco,qd_re_t,qd_im_t : quad_double;
    mp_re_t,mp_im_t : Floating_Number;
    dd_t : DoblDobl_Complex_Numbers.Complex_Number;
    qd_t : QuadDobl_Complex_Numbers.Complex_Number;
    mp_t : Multprec_Complex_Numbers.Complex_Number;
    prcn : natural32;

  begin
    if Symbol_Table.Number < dim
     then Symbol_Table.Init(dim);
    end if;
    Standard_Estimate_Loss_for_Polynomial_Homotopy(z,t,d_jfrco,d_fzrco,loss);
    if loss > -15 then
      jfrco := create(d_jfrco);
      fzrco := create(d_fzrco);
    else
      dd_re_t := create(Standard_Complex_Numbers.REAL_PART(t));
      dd_im_t := create(Standard_Complex_Numbers.IMAG_PART(t));
      dd_t := DoblDobl_Complex_Numbers.Create(dd_re_t,dd_im_t);
      DoblDobl_Estimate_Loss_for_Polynomial_Homotopy
        (z,dd_t,dd_jfrco,dd_fzrco,loss);
      if loss > -30 then
        jfrco := Multprec_DoblDobl_Convertors.to_floating_number(dd_jfrco);
        fzrco := Multprec_DoblDobl_Convertors.to_floating_number(dd_fzrco);
      else
        qd_re_t := create(Standard_Complex_Numbers.REAL_PART(t));
        qd_im_t := create(Standard_Complex_Numbers.IMAG_PART(t));
        qd_t := QuadDobl_Complex_Numbers.Create(qd_re_t,qd_im_t);
        QuadDobl_Estimate_Loss_for_Polynomial_Homotopy
          (z,qd_t,qd_jfrco,qd_fzrco,loss);
        if loss > -60 then
          jfrco := Multprec_QuadDobl_Convertors.to_floating_number(qd_jfrco);
          fzrco := Multprec_QuadDobl_Convertors.to_floating_number(qd_fzrco);
        else
          mp_re_t := create(Standard_Complex_Numbers.REAL_PART(t));
          mp_im_t := create(Standard_Complex_Numbers.IMAG_PART(t));
          mp_t := Multprec_Complex_Numbers.Create(mp_re_t,mp_im_t);
          prcn := 80;
          loop
            Multprec_Estimate_Loss_for_Polynomial_Homotopy
              (z,mp_t,prcn,jfrco,fzrco,loss);
            exit when (natural32(loss) < prcn);
            Clear(jfrco); Clear(fzrco);
            prcn := prcn + prcn/16;  -- increment precision with 25%
            exit when (prcn > maxprec);
          end loop;
        end if;
      end if;
    end if;
  end Estimate_Loss_for_Polynomial_Homotopy;

  function Estimate_Loss_for_Polynomial_Homotopy
              ( z : string; t : Standard_Complex_Numbers.Complex_Number;
                maxprec : natural32 ) return integer32 is

    res : integer32;
    jfrco,fzrco : Floating_Number;

  begin
    Estimate_Loss_for_Polynomial_Homotopy(z,t,maxprec,jfrco,fzrco,res);
    Clear(jfrco); Clear(fzrco);
    return res;
  end Estimate_Loss_for_Polynomial_Homotopy;

end Varbprec_Corrector_Steps;
