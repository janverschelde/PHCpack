with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Standard_Durand_Kerner;
--with Multprec_Durand_Kerner;

package body Hybrid_Durand_Kerner is

  function Eval ( p : Vector; x : Complex_Number ) return Complex_Number is

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
  end Eval;

  function Diff ( p : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns the coefficients of the derivative polynomial.

    res : Vector(0..p'last-1);

  begin
    for i in res'range loop
      res(i) := Create(i+1);
      Mul(res(i),p(i+1));
    end loop;
    return res;
  end Diff;

  generic

    with procedure Write ( step : in natural32; z,res : in Vector );

  procedure Reporting_Multprec_Newton_Refiner
              ( p : in Vector; z,res : in out Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Does a number of multi-precision refinements with Newton's method
  --   until the desired precision eps is reached for all roots or 
  --   until the number of steps is larger than maxsteps.

  procedure Reporting_Multprec_Newton_Refiner
              ( p : in Vector; z,res : in out Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean ) is

    nrm : Floating_Number;
    eva,der : Complex_Number;
    dp : Vector(0..p'last-1) := Diff(p);

  begin
    nb := maxsteps;
    for i in 1..maxsteps loop
      for j in z'range loop
        der := Eval(dp,z(j));     
        eva := res(j)/der;
        Sub(z(j),eva);
        Clear(eva); Clear(der); Clear(res(j));
        res(j) := Eval(p,z(j));
      end loop;
      Write(i,z,res);
      nrm := Max_Norm(res);
      fail := (nrm > eps);
      Clear(nrm);
      if not fail
       then nb := i; exit;
      end if;
    end loop;
    Clear(dp);
  end Reporting_Multprec_Newton_Refiner;

  procedure Silent_Multprec_Newton_Refiner
              ( p : in Vector; z,res : in out Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean ) is

    nrm : Floating_Number;
    eva,der : Complex_Number;
    dp : Vector(0..p'last-1) := Diff(p);

  begin
    nb := maxsteps;
    for i in 1..maxsteps loop
      for j in z'range loop
        der := Eval(dp,z(j));
        eva := res(j)/der;
        Sub(z(j),eva);
        Clear(eva); Clear(der); Clear(res(j));
        res(j) := Eval(p,z(j));
      end loop;
      nrm := Max_Norm(res);
      fail := (nrm > eps);
      Clear(nrm);
      if not fail
       then nb := i; exit;
      end if;
    end loop;
    Clear(dp);
  end Silent_Multprec_Newton_Refiner;

-- TARGET PROCEDURES :

  procedure Reporting_Durand_Kerner
              ( p : in Vector; z,res : in out Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean ) is

    stp : constant Standard_Complex_Vectors.Vector(p'range) := Round(p);
    stz : Standard_Complex_Vectors.Vector(z'range) := Round(z);
    stres : Standard_Complex_Vectors.Vector(res'range) := stz;
    st_eps : constant double_float := 1.0E-13;
    st_nb,mp_nb : natural32 := 0;
  
    procedure Standard_Write
                ( step : in natural32;
                  z,res : in Standard_Complex_Vectors.Vector ) is

      mpz : constant Multprec_Complex_Vectors.Vector(z'range) := Create(z);
      mpres : constant Multprec_Complex_Vectors.Vector(z'range) := Create(res);

    begin
      Write(step,mpz,mpres);
    end Standard_Write;

    procedure stdk is new
      Standard_Durand_Kerner.Reporting_Durand_Kerner(Standard_Write);

   -- procedure mpdk is new
   --   Multprec_Durand_Kerner.Reporting_Durand_Kerner(Write);

    procedure mpnr is new
      Reporting_Multprec_Newton_Refiner(Write);

  begin
    stdk(stp,stz,stres,maxsteps,st_eps,st_nb,fail);
    z := Create(stz);
    res := Create(stres);
   -- mpdk(p,z,res,maxsteps,eps,mp_nb,fail);
    mpnr(p,z,res,maxsteps,eps,mp_nb,fail);
    nb := st_nb + mp_nb;
  end Reporting_Durand_Kerner;

  procedure Silent_Durand_Kerner
              ( p : in Vector; z,res : in out Vector;
                maxsteps : in natural32; eps : in Floating_Number;
                nb : out natural32; fail : out boolean ) is

    stp : constant Standard_Complex_Vectors.Vector(p'range) := Round(p);
    stz : Standard_Complex_Vectors.Vector(z'range) := Round(z);
    stres : Standard_Complex_Vectors.Vector(res'range) := stz;
    st_eps : constant double_float := 1.0E-13;
    st_nb,mp_nb : natural32 := 0;

  begin
    Standard_Durand_Kerner.Silent_Durand_Kerner
      (stp,stz,stres,maxsteps,st_eps,st_nb,fail);
    z := Create(stz);
    res := Create(stres);
   -- Multprec_Durand_Kerner.Silent_Durand_Kerner
   --   (p,z,res,maxsteps,eps,mp_nb,fail);
    Silent_Multprec_Newton_Refiner
      (p,z,res,maxsteps,eps,mp_nb,fail);
    nb := st_nb + mp_nb;
  end Silent_Durand_Kerner;

end Hybrid_Durand_Kerner;
