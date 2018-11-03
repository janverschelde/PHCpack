with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;
with DoblDobl_Complex_Series_io;         use DoblDobl_Complex_Series_io;

package body DoblDobl_Algebraic_Series3 is

  function sqrt ( c : Series; i : natural32;
                  verbose : boolean := false ) return Series is

    res : Series(c.deg) := Create(0,c.deg);
    cpc : Series(c.deg) := c;
    wrk,dx : Series(c.deg);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(integer32(1));
    two : constant double_double := create(integer32(2));
    ctwo : constant Complex_Number := create(two);
    half : constant double_double := one/two;
    fac : constant Complex_Number := Create(half);
    tol : constant double_float := 1.0E-24;

  begin
    if AbsVal(cpc.cff(0)) < tol then
      res := Create(0,c.deg);
    else
      res := Create(0);
      res.cff(0) := DoblDobl_Complex_Numbers_Polar.Root(cpc.cff(0),2,i);
      for i in 0..c.deg loop
        wrk := res*res - cpc;
        dx := fac*wrk/res;
        if verbose then
          put("evaluation at degree = "); put(res.deg,1);
          put_line(" :"); put(wrk);
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

    res : Series(c.deg) := Create(0,c.deg);
    cpc : Series(c.deg) := c;
    wrk,dx : Series(c.deg);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    ddn : constant double_double := create(n);
    denominator : constant double_double := one/ddn;
    fac : constant Complex_Number := Create(denominator);

  begin
    res.cff(0) := DoblDobl_Complex_Numbers_Polar.Root(cpc.cff(0),n,i);
    for i in 0..c.deg loop
      wrk := res**integer(n) - cpc;
      dx := fac*wrk/(res**integer(n-1));
      if verbose then
        put("update dx at degree = "); put(res.deg,1);
        put_line(" :"); put(dx);
      end if;
      res := res - dx;
    end loop;
    return res;
  end Root;

end DoblDobl_Algebraic_Series3;
