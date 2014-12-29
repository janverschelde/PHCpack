with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;
with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Multprec_DoblDobl_Convertors;
with Multprec_QuadDobl_Convertors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vector_Strings;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Strings;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Strings;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vector_Strings;
with Multprec_Complex_Vector_Tools;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;
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

-- PART III : computing after parsing to the precision

  procedure Standard_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                err,rco,res : out double_float ) is

    h : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.Standard_Homotopy_System;
    n : constant integer32 := h'last;
    x : Standard_Complex_Vectors.Vector(1..n)
      := Standard_Complex_Vector_Strings.Parse(z.all);
    fz : constant Standard_Complex_Vectors.Vector(h'range)
       := Varbprec_Homotopy.Eval(x,t);
    jfz : Standard_Complex_Matrices.Matrix(h'range,x'range)
        := Varbprec_Homotopy.Diff(x,t);
    piv : Standard_Integer_Vectors.Vector(x'range);

  begin
    lufco(jfz,n,piv,rco);
    Varbprec_Complex_Newton_Steps.do_Newton_Step(x,jfz,piv,fz,err);
    res := Max_Norm(fz);
    declare
      nz : constant string
         := Standard_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
  end Standard_Newton_Step_on_Polynomial_Homotopy;

  procedure DoblDobl_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                err,rco,res : out double_float ) is

    h : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.DoblDobl_Homotopy_System;
    n : constant integer32 := h'last;
    x : DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Complex_Vector_Strings.Parse(z.all);
    dd_re : constant double_double
          := Create(Standard_Complex_Numbers.REAL_PART(t));
    dd_im : constant double_double
          := Create(Standard_Complex_Numbers.IMAG_PART(t));
    dd_t : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(dd_re,dd_im);
    fz : constant DoblDobl_Complex_Vectors.Vector(h'range)
       := Varbprec_Homotopy.Eval(x,dd_t);
    jfz : DoblDobl_Complex_Matrices.Matrix(h'range,x'range)
        := Varbprec_Homotopy.Diff(x,dd_t);
    piv : Standard_Integer_Vectors.Vector(x'range);
    dd_err,dd_rco,dd_res : double_double;

  begin
    lufco(jfz,n,piv,dd_rco);
    Varbprec_Complex_Newton_Steps.do_Newton_Step(x,jfz,piv,fz,dd_err);
    dd_res := Max_Norm(fz);
    declare
      nz : constant string
         := DoblDobl_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    err := to_double(dd_err);
    rco := to_double(dd_rco);
    res := to_double(dd_res);
  end DoblDobl_Newton_Step_on_Polynomial_Homotopy;

  procedure QuadDobl_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                err,rco,res : out double_float ) is

    h : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.DoblDobl_Homotopy_System;
    n : constant integer32 := h'last;
    x : QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Complex_Vector_Strings.Parse(z.all);
    qd_re : constant quad_double
          := Create(Standard_Complex_Numbers.REAL_PART(t));
    qd_im : constant quad_double
          := Create(Standard_Complex_Numbers.IMAG_PART(t));
    qd_t : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(qd_re,qd_im);
    fz : constant QuadDobl_Complex_Vectors.Vector(h'range)
       := Varbprec_Homotopy.Eval(x,qd_t);
    jfz : QuadDobl_Complex_Matrices.Matrix(h'range,x'range)
        := Varbprec_Homotopy.Diff(x,qd_t);
    piv : Standard_Integer_Vectors.Vector(x'range);
    qd_err,qd_rco,qd_res : quad_double;

  begin
    lufco(jfz,n,piv,qd_rco);
    Varbprec_Complex_Newton_Steps.do_Newton_Step(x,jfz,piv,fz,qd_err);
    qd_res := Max_Norm(fz);
    declare
      nz : constant string
         := QuadDobl_Complex_Vector_Strings.Write(x);
    begin
      Clear(z);
      z := new string'(nz);
    end;
    err := to_double(qd_err);
    rco := to_double(qd_rco);
    res := to_double(qd_res);
  end QuadDobl_Newton_Step_on_Polynomial_Homotopy;

  procedure Multprec_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                prcn : in natural32; err,rco,res : out double_float ) is

    size : constant natural32 := Decimal_to_Size(prcn);
    h : constant Multprec_Complex_Poly_Systems.Link_to_Poly_Sys
      := Varbprec_Homotopy.Multprec_Homotopy_System(prcn);
    n : constant integer32 := h'last;
    x : Multprec_Complex_Vectors.Vector(1..n)
      := Multprec_Complex_Vector_Strings.Parse(z.all);
    fz : Multprec_Complex_Vectors.Vector(h'range);
    jfz : Multprec_Complex_Matrices.Matrix(h'range,x'range);
    piv : Standard_Integer_Vectors.Vector(x'range);
    mp_re : Floating_Number := Create(Standard_Complex_Numbers.REAL_PART(t));
    mp_im : Floating_Number := Create(Standard_Complex_Numbers.IMAG_PART(t));
    mp_t : Multprec_Complex_Numbers.Complex_Number
         :=  Multprec_Complex_Numbers.Create(mp_re,mp_im);
    mp_rco,mp_err,mp_res : Floating_Number;

  begin
    Multprec_Complex_Vector_Tools.Set_Size(x,size);
    fz := Varbprec_Homotopy.Eval(x,mp_t,prcn);
    jfz := Varbprec_Homotopy.Diff(x,mp_t,prcn);
    lufco(jfz,n,piv,mp_rco);
    rco := Round(mp_rco);
    Varbprec_Complex_Newton_Steps.do_Newton_Step(x,jfz,piv,fz,mp_err);
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
    Multprec_Complex_Matrices.Clear(jfz);
    Multprec_Complex_Vectors.Clear(fz);
    Multprec_Complex_Vectors.Clear(x);
    Clear(mp_re); Clear(mp_im);
    Multprec_Complex_Numbers.Clear(mp_t);
    Clear(mp_rco); Clear(mp_err); Clear(mp_res);
  end Multprec_Newton_Step_on_Polynomial_Homotopy;

  procedure do_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                loss : in integer32; want : in natural32;
                err,rco,res : out double_float ) is

    precision : constant natural32 := natural32(-loss) + want;

  begin
    if precision <= 16 then
      Standard_Newton_Step_on_Polynomial_Homotopy(z,t,err,rco,res);
    elsif precision <= 32 then
      DoblDobl_Newton_Step_on_Polynomial_Homotopy(z,t,err,rco,res);
    elsif precision <= 64 then
      QuadDobl_Newton_Step_on_Polynomial_Homotopy(z,t,err,rco,res);
    else
      Multprec_Newton_Step_on_Polynomial_Homotopy(z,t,precision,err,rco,res);
    end if;
  end do_Newton_Step_on_Polynomial_Homotopy;

-- PART IV : sequence of Newton steps to the wanted accuracy

  procedure Newton_Steps_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                want,maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float ) is

    precision : natural32;
    err_accu,res_accu : integer32;

  begin
    for i in 1..maxitr loop
      loss := Estimate_Loss_for_Polynomial_Homotopy(z.all,t,maxprc);
      precision := natural32(-loss) + want;
      do_Newton_Step_on_Polynomial_Homotopy(z,t,loss,want,err,rco,res);
      if err = 0.0
       then err_accu := integer32(precision);
       else err_accu := abs(integer32(log10(err)));
      end if;
      if res = 0.0
       then res_accu := integer32(precision);
       else res_accu := abs(integer32(log10(res)));
      end if;
      exit when ((err_accu >= integer32(want))
             and (res_accu >= integer32(want)));
    end loop;
  end Newton_Steps_on_Polynomial_Homotopy;

  procedure Newton_Steps_on_Polynomial_Homotopy
              ( file : in file_type; z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                want,maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float ) is

    precision : natural32;
    err_accu,res_accu : integer32;

  begin
    for i in 1..maxitr loop
      put(file,"estimated loss in step "); put(file,i,1); put(file," : ");
      loss := Estimate_Loss_for_Polynomial_Homotopy(z.all,t,maxprc);
      put(file,loss,1);
      precision := natural32(-loss) + want;
      put(file,", precision : "); put(file,precision,1); new_line(file);
      do_Newton_Step_on_Polynomial_Homotopy(z,t,loss,want,err,rco,res);
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
      exit when ((err_accu >= integer32(want))
             and (res_accu >= integer32(want)));
    end loop;
  end Newton_Steps_on_Polynomial_Homotopy;

end Varbprec_Corrector_Steps;
