with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with HexaDobl_Complex_Numbers;           use HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers_Polar;
with HexaDobl_Complex_Series_io;         use HexaDobl_Complex_Series_io;

package body HexaDobl_Complex_Algebraic_Series is

  function sqrt ( c : Series; i : natural32;
                  verbose : boolean := false ) return Series is

    res : Series(c.deg) := Create(0,c.deg);
    cpc : constant Series(c.deg) := c;
    wrk,dx : Series(c.deg);
    one : constant hexa_double := create(integer32(1));
    two : constant hexa_double := create(integer32(2));
    half : constant hexa_double := one/two;
    fac : constant Complex_Number := Create(half);
    tol : constant double_float := 1.0E-48;

  begin
    if AbsVal(cpc.cff(0)) < tol then
      res := Create(0,c.deg);
    else
      res := Create(0,c.deg);
      res.cff(0) := HexaDobl_Complex_Numbers_Polar.Root(cpc.cff(0),2,i);
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
    cpc : constant Series(c.deg) := c;
    wrk,dx : Series(c.deg);
    one : constant hexa_double := create(1.0);
    ddn : constant hexa_double := create(n);
    denominator : constant hexa_double := one/ddn;
    fac : constant Complex_Number := Create(denominator);

  begin
    res.cff(0) := HexaDobl_Complex_Numbers_Polar.Root(cpc.cff(0),n,i);
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

end HexaDobl_Complex_Algebraic_Series;
