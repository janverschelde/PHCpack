with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers_Polar;

package body Convergence_Radius_Estimates is

  function Is_Zero ( z : Standard_Complex_Numbers.Complex_Number )
                   return boolean is

    use Standard_Complex_Numbers;

  begin
    if REAL_PART(z) + 1.0 /= 1.0
     then return false;
     else return (IMAG_PART(z) + 1.0 = 1.0);
    end if;
  end Is_Zero;

  function Is_Zero ( z : DoblDobl_Complex_Numbers.Complex_Number )
                   return boolean is

    use DoblDobl_Complex_Numbers;

    one : constant double_double := create(1.0);

  begin
    if REAL_PART(z) + one /= one
     then return false;
     else return (IMAG_PART(z) + one = one);
    end if;
  end Is_Zero;

  function Is_Zero ( z : QuadDobl_Complex_Numbers.Complex_Number )
                   return boolean is

    use QuadDobl_Complex_Numbers;

    one : constant quad_double := create(1.0);

  begin
    if REAL_PART(z) + one /= one
     then return false;
     else return (IMAG_PART(z) + one = one);
    end if;
  end Is_Zero;

  procedure Fabry ( c : in Standard_Complex_Vectors.Vector;
                    z : out Standard_Complex_Numbers.Complex_Number;
                    e : out double_float; fail : out boolean ) is

    use Standard_Complex_Numbers;

  begin
    fail := Is_Zero(c(c'last));
    if not fail then
      z := c(c'last-1)/c(c'last);
      if Is_Zero(c(c'last-1))
       then e := 1.0;
       else e := AbsVal(z - c(c'last-2)/c(c'last-1));
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in DoblDobl_Complex_Vectors.Vector;
                    z : out DoblDobl_Complex_Numbers.Complex_Number;
                    e : out double_double; fail : out boolean ) is

    use DoblDobl_Complex_Numbers;

  begin
    fail := Is_Zero(c(c'last));
    if not fail then
      z := c(c'last-1)/c(c'last);
      if Is_Zero(c(c'last-1))
       then e := create(1.0);
       else e := AbsVal(z - c(c'last-2)/c(c'last-1));
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in QuadDobl_Complex_Vectors.Vector;
                    z : out QuadDobl_Complex_Numbers.Complex_Number;
                    e : out quad_double; fail : out boolean ) is

    use QuadDobl_Complex_Numbers;

  begin
    fail := Is_Zero(c(c'last));
    if not fail then
      z := c(c'last-1)/c(c'last);
      if Is_Zero(c(c'last-1))
       then e := create(1.0);
       else e := AbsVal(z - c(c'last-2)/c(c'last-1));
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in Standard_Complex_VecVecs.VecVec;
                    z : out Standard_Complex_Numbers.Complex_Number;
                    r : out double_float;
                    e : out double_float; fail : out boolean;
                    verbose : in boolean := true ) is

    use Standard_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : double_float;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail);
      if verbose then
        if kfail
         then put_line("zero last coefficient");
         else put(zk); put("  error estimate :"); put(ek,3); new_line;
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := Standard_Complex_Numbers_Polar.Radius(z);
        else
          rad := Standard_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( c : in DoblDobl_Complex_VecVecs.VecVec;
                    z : out DoblDobl_Complex_Numbers.Complex_Number;
                    r : out double_double;
                    e : out double_double; fail : out boolean;
                    verbose : in boolean := true ) is

    use DoblDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : double_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail);
      if verbose then
        if kfail
         then put_line("zero last coefficient");
         else put(zk); put("  error estimate : "); put(ek,3); new_line;
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := DoblDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := DoblDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( c : in QuadDobl_Complex_VecVecs.VecVec;
                    z : out QuadDobl_Complex_Numbers.Complex_Number;
                    r : out quad_double;
                    e : out quad_double; fail : out boolean;
                    verbose : in boolean := true ) is

    use QuadDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : quad_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail);
      if verbose then
        if kfail
         then put_line("zero last coefficient");
         else put(zk); put("  error estimate : "); put(ek,3); new_line;
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := QuadDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := QuadDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

end Convergence_Radius_Estimates;
