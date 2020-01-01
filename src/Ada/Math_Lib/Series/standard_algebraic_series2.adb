with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Standard_Dense_Series2_io;          use Standard_Dense_Series2_io;

package body Standard_Algebraic_Series2 is

  function sqrt ( c : Series; i : natural32;
                  verbose : boolean := false ) return Series is

    res : Series(c.deg) := Create(0.0,c.deg);
    cpc : constant Series(c.deg) := c;
    wrk,dx : Series(c.deg);
    fac : constant Complex_Number := Create(0.5);
    tol : constant double_float := 1.0E-13;

  begin
    if AbsVal(cpc.cff(0)) < tol then
      res := Create(0.0,c.deg);
    else
      res := Create(0.0);
      res.cff(0) := Standard_Complex_Numbers_Polar.Root(cpc.cff(0),2,i);
      for i in 0..c.deg loop
        wrk := res*res - cpc;
        dx := fac*wrk/res;
        if verbose then
          put("update dx at degree = "); put(res.deg,1);
          put_line(" :"); put(dx);
        end if;
        res := res - dx;
      end loop;
    end if;
    return res;
  end sqrt;

  function Root ( c : Series; n,i : natural32;
                  verbose : boolean := false ) return Series is

    res : Series(c.deg) := Create(0.0,c.deg);
    cpc : constant Series(c.deg) := c;
    wrk,dx : Series(c.deg);
    denominator : constant double_float := 1.0/double_float(n);
    fac : constant Complex_Number := Create(denominator);

  begin
    res.cff(0) := Standard_Complex_Numbers_Polar.Root(cpc.cff(0),n,i);
    for i in 0..c.deg loop
      wrk := res**n - cpc;
      dx := fac*wrk/(res**(n-1));
      if verbose then
        put("update dx at degree = "); put(res.deg,1);
        put_line(" :"); put(dx);
      end if;
      res := res - dx;
    end loop;
    return res;
  end Root;

end Standard_Algebraic_Series2;
