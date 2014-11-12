with Standard_Complex_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with DoblDobl_Complex_Numbers;
with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with DoblDobl_Complex_Poly_Functions;    use DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Numbers;
with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Poly_Functions;    use QuadDobl_Complex_Poly_Functions;
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
                
end Varbprec_Complex_Newton_Steps;
