with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;

package body QuadDobl_Durand_Kerner is

  function Horner ( p : QuadDobl_Complex_Vectors.Vector;
                    x : Complex_Number ) return Complex_Number is

    res : Complex_Number;

  begin
    if p'first = p'last then
      res := p(p'first);
    else
      res := p(p'last);
      for i in reverse p'first..(p'last-1) loop
        res := res*x + p(i);
      end loop;
    end if;
    return res;
  end Horner;

  function Derivative ( p : QuadDobl_Complex_Vectors.Vector )
                      return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(0..p'last-1);
    qd_i : quad_double;

  begin
    for i in 1..p'last loop
      qd_i := create(i);
      res(i-1) := qd_i*p(i);
    end loop;
    return res;
  end Derivative;

  procedure Newton ( p,dp : in QuadDobl_Complex_Vectors.Vector;
                     z : in out Complex_Number;
                     err,rco,res : out double_float ) is

    pz : constant Complex_Number := Horner(p,z);
    dpz : constant Complex_Number := Horner(dp,z);
    dz : constant Complex_Number := pz/dpz;
    qd_rco,qd_err,qd_res : quad_double;

  begin
    qd_rco := AbsVal(dpz); rco := to_double(qd_rco);
    qd_err := AbsVal(dz);  err := to_double(qd_err);
    qd_res := AbsVal(pz);  res := to_double(qd_res);
    z := z - dz;
  end Newton;

  procedure Newton ( p,dp : in QuadDobl_Complex_Vectors.Vector;
                     z : in out QuadDobl_Complex_Vectors.Vector;
                     err,rco,res : out Standard_Floating_Vectors.Vector ) is
  begin
    for i in z'range loop
      Newton(p,dp,z(i),err(i),rco(i),res(i));
    end loop;
  end Newton;

  procedure DK ( p : in QuadDobl_Complex_Vectors.Vector;
                 z,r : in out QuadDobl_Complex_Vectors.Vector ) is

    function Compute_q ( i : integer32;
                         a : QuadDobl_Complex_Vectors.Vector )
                       return Complex_Number is

    -- DESCRIPTION :
    --   Computes the quotient needed in the Durand-Kerner step.

      res : Complex_Number;

    begin
      res := create(integer(1));
      for j in a'range loop
        if j /= i
         then res := res*(a(i)-a(j));
        end if;
      end loop;
      return res;
    end Compute_q;

  begin
    for i in z'range loop
      z(i) := z(i) - Horner(p,z(i))/Compute_q(i,z);
      r(i) := Horner(p,z(i));
    end loop;
  end DK;

  procedure Reporting_Durand_Kerner
              ( p : in QuadDobl_Complex_Vectors.Vector;
                z,r : in out QuadDobl_Complex_Vectors.Vector;
                maxsteps : in natural32; eps : in double_float;
                nb : out natural32; fail : out boolean ) is

    pp : QuadDobl_Complex_Vectors.Vector(p'range);
    one : constant Complex_Number := create(integer(1));

  begin
    if p(p'last) /= one then
      for i in p'range loop
        pp(i) := p(i)/p(p'last);
      end loop;
    else
      for i in p'range loop
        pp(i) := p(i);
      end loop;
    end if;
    fail := true;
    for k in 1..maxsteps loop
      nb := k;
      DK(pp,z,r);  
      Write(k,z,r);
      if Max_Norm(r) <= eps
       then fail := false; exit;
      end if;
    end loop;
  exception
    when others => fail := true;
  end Reporting_Durand_Kerner;

  procedure Silent_Durand_Kerner
              ( p : in QuadDobl_Complex_Vectors.Vector;
                z,r : in out QuadDobl_Complex_Vectors.Vector;
                maxsteps : in natural32; eps : in double_float;
                nb : out natural32; fail : out boolean ) is

    pp : QuadDobl_Complex_Vectors.Vector(p'range);
    one : constant Complex_Number := create(integer(1));

  begin
    if p(p'last) /= one then
      for i in p'range loop
        pp(i) := p(i)/p(p'last);
      end loop;
    else
      for i in p'range loop
        pp(i) := p(i);
      end loop;
    end if;
    fail := true;
    for k in 1..maxsteps loop
      nb := k;
      DK(pp,z,r);  
      if Max_Norm(r) <= eps
       then fail := false; exit;
      end if;
    end loop;
  exception
    when others => fail := true;
  end Silent_Durand_Kerner;

end QuadDobl_Durand_Kerner;
