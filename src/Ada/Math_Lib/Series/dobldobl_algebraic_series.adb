with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
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
    one : constant double_double := create(1.0);
    two : constant double_double := create(2.0);
    half : constant double_double := one/two;
    fac : constant Complex_Number := Create(half);

  begin
    res.cff(0) := DoblDobl_Complex_Numbers_Polar.Root(cpc.cff(0),2,i);
    res.cff(1) := Create(zero);
    res.order := 1;
    cpc.order := 1;
    loop
      wrk := res*res - cpc;
      dx := fac*wrk/res;
      if verbose then
        put("update dx at order = "); put(res.order,1);
        put_line(" :"); put(dx);
      end if;
      res := res - dx;
      exit when res.order >= c.order;
      cpc.order := 2*cpc.order;
      res.order := 2*res.order;
    end loop;
    if res.order = c.order
     then return res;
     else return Create(res,c.order);
    end if;
  end sqrt;

  function Root ( c : Series; n,i : natural32;
                  verbose : boolean := false ) return Series is

    res : Series;
    cpc : Series := c;
    wrk,dx : Series;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    denominator : constant double_double := one/create(n);
    fac : constant Complex_Number := Create(denominator);

  begin
    res.cff(0) := DoblDobl_Complex_Numbers_Polar.Root(cpc.cff(0),n,i);
    res.cff(1) := Create(zero);
    res.order := 1;
    cpc.order := 1;
    loop
      wrk := res**n - cpc;
      dx := fac*wrk/(res**(n-1));
      if verbose then
        put("update dx at order = "); put(res.order,1);
        put_line(" :"); put(dx);
      end if;
      res := res - dx;
      exit when res.order >= c.order;
      cpc.order := 2*cpc.order;
      res.order := 2*res.order;
    end loop;
    if res.order = c.order
     then return res;
     else return Create(res,c.order);
    end if;
  end Root;

end DoblDobl_Algebraic_Series;
