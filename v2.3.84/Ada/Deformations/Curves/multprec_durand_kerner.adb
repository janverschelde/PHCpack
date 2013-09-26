with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;

package body Multprec_Durand_Kerner is

-- AUXILIARIES :

  function Horner ( p : Vector; x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns (..((a[n]*x + a[n-1])*x + a[n-2])*x + .. + a[1])*x + a[0].

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

  procedure DK ( p : in Vector; z,res : in out Vector;
                 start : in natural32 ) is

  -- DESCRIPTION :
  --   Computes one step in the Durand-Kerner iteration.

    eva,fac : Complex_Number;

    function Compute_q ( i : integer32; a : Vector ) return Complex_Number is

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
              ( p : in Vector; z,res : in out Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean ) is

    pp : Vector(p'range);
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
              ( p : in Vector; z,res : in out Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean ) is

    pp : Vector(p'range);
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
