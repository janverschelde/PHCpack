with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;

package body Multprec_Durand_Kerner is

  function Horner ( p : Multprec_Complex_Vectors.Vector;
                    x : Complex_Number ) return Complex_Number is

    res : Complex_Number;

  begin
    if p'first = p'last then
      Copy(p(p'first),res);
    else
      Copy(p(p'last),res);
      for i in reverse p'first..(p'last-1) loop
        Mul(res,x);
        Add(res,p(i));
      end loop;
    end if;
    return res;
  end Horner;

  function Derivative ( p : Multprec_Complex_Vectors.Vector )
                      return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(0..p'last-1);
    mpi : Floating_Number;
    mpc : Complex_Number;

  begin
    for i in 1..p'last loop
      mpi := Create(double_float(i));
      mpc := Multprec_Complex_Numbers.Create(mpi);
      res(i-1) := mpc*p(i);
      Clear(mpi);
      Multprec_Complex_Numbers.Clear(mpc);
    end loop;
    return res;
  end Derivative;

  procedure Newton ( p,dp : in Multprec_Complex_Vectors.Vector;
                     z : in out Complex_Number;
                     err,rco,res : out Floating_Number ) is

    pz : Complex_Number := Horner(p,z);
    dpz : Complex_Number := Horner(dp,z);
    dz : Complex_Number := pz/dpz;

  begin
    rco := AbsVal(dpz);
    err := AbsVal(dz);
    res := AbsVal(pz);
    Sub(z,dz);
    Clear(pz); Clear(dpz); Clear(dz);
  end Newton;

  procedure Newton ( p,dp : in Multprec_Complex_Vectors.Vector;
                     z : in out Multprec_Complex_Vectors.Vector;
                     err,rco,res : out Multprec_Floating_Vectors.Vector ) is
  begin
    for i in z'range loop
      Newton(p,dp,z(i),err(i),rco(i),res(i));
    end loop;
  end Newton;

-- AUXILIARIES :

  procedure DK ( p : in Multprec_Complex_Vectors.Vector;
                 z,res : in out Multprec_Complex_Vectors.Vector;
                 start : in natural32 ) is

  -- DESCRIPTION :
  --   Computes one step in the Durand-Kerner iteration.

    eva,fac : Complex_Number;

    function Compute_q 
               ( i : integer32;
                 a : Multprec_Complex_Vectors.Vector )
               return Complex_Number is

    -- DESCRIPTION :
    --   Computes the quotient needed in the Durand-Kerner step.

      res,diff : Complex_Number;

    begin
      res := Create(integer(1));
      for j in a'range loop
        if j /= i then
          diff := a(i) - a(j);
          Mul(res,diff);
          Clear(diff);
        end if;
      end loop;
      return res;
    end Compute_q;

  begin
    for i in integer32(start)..z'last loop
      eva := Horner(p,z(i));
      fac := Compute_q(i,z);
      Div(eva,fac);
      Sub(z(i),eva);
      Clear(fac); Clear(eva); Clear(res(i));
      res(i) := Horner(p,z(i));
    end loop;
    for i in z'first..(integer32(start)-1) loop
      eva := Horner(p,z(i));
      fac := Compute_q(i,z);
      Div(eva,fac);
      Sub(z(i),eva);
      Clear(fac); Clear(eva); Clear(res(i));
      res(i) := Horner(p,z(i));
    end loop;
  end DK;

-- TARGET ROUTINES :

  procedure Reporting_Durand_Kerner
              ( p : in Multprec_Complex_Vectors.Vector;
                z,res : in out Multprec_Complex_Vectors.Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean ) is

    pp : Multprec_Complex_Vectors.Vector(p'range);
    nrm : Floating_Number;
    start : natural32;
    stop : boolean := false;

  begin
    for i in res'range loop
      res(i) := Create(integer(1));
    end loop;
    if p(p'last) /= Create(integer(1)) then
      for i in p'range loop
        pp(i) := p(i) / p(p'last);
      end loop;
      for k in 1..maxsteps loop
        nb := k;
        start := k mod natural32(z'last);
        if start = 0
         then start := 1;
        end if;
        DK(pp,z,res,start);  
        Write(k,z,res);
        nrm := Max_Norm(res);
        stop := (nrm < eps);
        Clear(nrm);
        exit when stop; 
      end loop;
      Clear(pp);
    else
      for k in 1..maxsteps loop
        nb := k;
        start := k mod natural32(z'last);
        if start = 0
         then start := 1;
        end if;
        DK(p,z,res,start);  
        Write(k,z,res);
        nrm := Max_Norm(res);
        stop := (nrm < eps);
        Clear(nrm);
        exit when stop;
      end loop;
    end if;
    fail := not stop;
  end Reporting_Durand_Kerner;

  procedure Silent_Durand_Kerner
              ( p : in Multprec_Complex_Vectors.Vector;
                z,res : in out Multprec_Complex_Vectors.Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean ) is

    pp : Multprec_Complex_Vectors.Vector(p'range);
    nrm : Floating_Number;
    start : natural32;
    stop : boolean := false;

  begin
    for i in res'range loop
      res(i) := Create(integer(1));
    end loop;
    if p(p'last) /= Create(integer(1)) then
      for i in p'range loop
        pp(i) := p(i) / p(p'last);
      end loop;
      for k in 1..maxsteps loop
        nb := k;
        start := k mod natural32(z'last);
        if start = 0
         then start := 1;
        end if;
        DK(pp,z,res,start);  
        nrm := Max_Norm(res);
        stop := (nrm < eps);
        Clear(nrm);
        exit when stop;
      end loop;
      Clear(pp);
    else
      for k in 1..maxsteps loop
        nb := k;
        start := k mod natural32(z'last);
        if start = 0
         then start := 1;
        end if;
        DK(p,z,res,start);  
        nrm := Max_Norm(res);
        stop := (nrm < eps);
        Clear(nrm);
        exit when stop;
      end loop;
    end if;
    fail := not stop;
  end Silent_Durand_Kerner;

end Multprec_Durand_Kerner;
