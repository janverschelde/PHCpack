with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;
with DoblDobl_Dense_Series_io;           use DoblDobl_Dense_Series_io;

package body DoblDobl_Algebraic_Series is

  function sqrt ( c : Series; i : natural32;
                  verbose : boolean := false ) return Series is

    res : Series;
    cpc : Series := c;
    wrk,dx : Series;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(integer32(1));
    two : constant double_double := create(integer32(2));
    half : constant double_double := one/two;
    fac : constant Complex_Number := Create(half);
    tol : constant double_float := 1.0E-24;

  begin
    if AbsVal(cpc.cff(0)) < tol then
      res := Create(0.0);
      return res;
    else
      res.cff(0) := DoblDobl_Complex_Numbers_Polar.Root(cpc.cff(0),2,i);
      res.cff(1) := Create(zero);
      res.deg := 1;
      cpc.deg := 1;
      loop
        wrk := res*res - cpc;
        dx := fac*wrk/res;
        if verbose then
          put("evaluation at degree = "); put(res.deg,1);
          put_line(" :"); put(wrk);
          put("update dx at degree = "); put(res.deg,1);
          put_line(" :"); put(dx);
        end if;
        res := res - dx;
        exit when res.deg >= c.deg;
        cpc.deg := 2*cpc.deg;
        res.deg := 2*res.deg;
      end loop;
      if res.deg = c.deg
       then return res;
       else return Create(res,c.deg);
      end if;
    end if;
  end sqrt;

  function Root ( c : Series; n,i : natural32;
                  verbose : boolean := false ) return Series is

    res : Series;
    cpc : Series := c;
    wrk,dx : Series;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    ddn : constant double_double := create(n);
    denominator : constant double_double := one/ddn;
    fac : constant Complex_Number := Create(denominator);

  begin
    res.cff(0) := DoblDobl_Complex_Numbers_Polar.Root(cpc.cff(0),n,i);
    res.cff(1) := Create(zero);
    res.deg := 1;
    cpc.deg := 1;
    loop
      wrk := res**n - cpc;
      dx := fac*wrk/(res**(n-1));
      if verbose then
        put("update dx at degree = "); put(res.deg,1);
        put_line(" :"); put(dx);
      end if;
      res := res - dx;
      exit when res.deg >= c.deg;
      cpc.deg := 2*cpc.deg;
      res.deg := 2*res.deg;
    end loop;
    if res.deg = c.deg
     then return res;
     else return Create(res,c.deg);
    end if;
  end Root;

end DoblDobl_Algebraic_Series;
