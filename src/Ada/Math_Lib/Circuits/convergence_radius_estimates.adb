with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Penta_Double_Numbers_io;            use Penta_Double_Numbers_io;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with TripDobl_Complex_Numbers_io;        use TripDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with PentDobl_Complex_Numbers_io;        use PentDobl_Complex_Numbers_io;
with OctoDobl_Complex_Numbers_io;        use OctoDobl_Complex_Numbers_io;
with DecaDobl_Complex_Numbers_io;        use DecaDobl_Complex_Numbers_io;
with HexaDobl_Complex_Numbers_io;        use HexaDobl_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers_Polar;
with TripDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers_Polar;
with PentDobl_Complex_Numbers_Polar;
with OctoDobl_Complex_Numbers_Polar;
with DecaDobl_Complex_Numbers_Polar;
with HexaDobl_Complex_Numbers_Polar;

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

  function Is_Zero ( z : TripDobl_Complex_Numbers.Complex_Number )
                   return boolean is

    use TripDobl_Complex_Numbers;

    one : constant triple_double := create(1.0);

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

  function Is_Zero ( z : PentDobl_Complex_Numbers.Complex_Number )
                   return boolean is

    use PentDobl_Complex_Numbers;

    one : constant penta_double := create(1.0);

  begin
    if REAL_PART(z) + one /= one
     then return false;
     else return (IMAG_PART(z) + one = one);
    end if;
  end Is_Zero;

  function Is_Zero ( z : OctoDobl_Complex_Numbers.Complex_Number )
                   return boolean is

    use OctoDobl_Complex_Numbers;

    one : constant octo_double := create(1.0);

  begin
    if REAL_PART(z) + one /= one
     then return false;
     else return (IMAG_PART(z) + one = one);
    end if;
  end Is_Zero;

  function Is_Zero ( z : DecaDobl_Complex_Numbers.Complex_Number )
                   return boolean is

    use DecaDobl_Complex_Numbers;

    one : constant deca_double := create(1.0);

  begin
    if REAL_PART(z) + one /= one
     then return false;
     else return (IMAG_PART(z) + one = one);
    end if;
  end Is_Zero;

  function Is_Zero ( z : HexaDobl_Complex_Numbers.Complex_Number )
                   return boolean is

    use HexaDobl_Complex_Numbers;

    one : constant hexa_double := create(1.0);

  begin
    if REAL_PART(z) + one /= one
     then return false;
     else return (IMAG_PART(z) + one = one);
    end if;
  end Is_Zero;

  procedure Fabry ( c : in Standard_Complex_Vectors.Vector;
                    z : out Standard_Complex_Numbers.Complex_Number;
                    e : out double_float; fail : out boolean;
                    offset : in integer32 := 0 ) is

    use Standard_Complex_Numbers;
    nm1 : constant double_float := double_float(c'last - 1);

  begin
    fail := Is_Zero(c(c'last-offset));
    if not fail then
      if offset = 0 then
        z := c(c'last-1)/c(c'last);
        if Is_Zero(c(c'last-1))
         then e := 1.0;
         else e := AbsVal(z - c(c'last-2)/c(c'last-1))*nm1;
        end if;
      else -- offset should be > 0, typically 2
        z := c(c'last-1-offset)/c(c'last-offset);
        if Is_Zero(c(c'last))
         then e := 1.0;
         else e := AbsVal(z - c(c'last-1)/c(c'last))*nm1;
        end if;
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in DoblDobl_Complex_Vectors.Vector;
                    z : out DoblDobl_Complex_Numbers.Complex_Number;
                    e : out double_double; fail : out boolean;
                    offset : in integer32 := 0 ) is

    use DoblDobl_Complex_Numbers;
    nm1 : constant double_float := double_float(c'last - 1);

  begin
    fail := Is_Zero(c(c'last-offset));
    if not fail then
      if offset = 0 then
        z := c(c'last-1)/c(c'last);
        if Is_Zero(c(c'last-1))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-2)/c(c'last-1))*nm1;
        end if;
      else -- offset should be > 0, typically 2
        z := c(c'last-1-offset)/c(c'last-offset);
        if Is_Zero(c(c'last))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-1)/c(c'last))*nm1;
        end if;
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in TripDobl_Complex_Vectors.Vector;
                    z : out TripDobl_Complex_Numbers.Complex_Number;
                    e : out triple_double; fail : out boolean;
                    offset : in integer32 := 0 ) is

    use TripDobl_Complex_Numbers;
    nm1 : constant double_float := double_float(c'last - 1);

  begin
    fail := Is_Zero(c(c'last-offset));
    if not fail then
      if offset = 0 then
        z := c(c'last-1)/c(c'last);
        if Is_Zero(c(c'last-1))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-2)/c(c'last-1))*nm1;
        end if;
      else -- offset should be > 0, typically 2
        z := c(c'last-1-offset)/c(c'last-offset);
        if Is_Zero(c(c'last))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-1)/c(c'last))*nm1;
        end if;
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in QuadDobl_Complex_Vectors.Vector;
                    z : out QuadDobl_Complex_Numbers.Complex_Number;
                    e : out quad_double; fail : out boolean;
                    offset : in integer32 := 0 ) is

    use QuadDobl_Complex_Numbers;
    nm1 : constant double_float := double_float(c'last - 1);

  begin
    fail := Is_Zero(c(c'last-offset));
    if not fail then
      if offset = 0 then
        z := c(c'last-1)/c(c'last);
        if Is_Zero(c(c'last-1))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-2)/c(c'last-1))*nm1;
        end if;
      else -- offset should be > 0, typically 2
        z := c(c'last-1-offset)/c(c'last-offset);
        if Is_Zero(c(c'last))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-1)/c(c'last))*nm1;
        end if;
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in PentDobl_Complex_Vectors.Vector;
                    z : out PentDobl_Complex_Numbers.Complex_Number;
                    e : out penta_double; fail : out boolean;
                    offset : in integer32 := 0 ) is

    use PentDobl_Complex_Numbers;
    nm1 : constant double_float := double_float(c'last - 1);

  begin
    fail := Is_Zero(c(c'last-offset));
    if not fail then
      if offset = 0 then
        z := c(c'last-1)/c(c'last);
        if Is_Zero(c(c'last-1))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-2)/c(c'last-1))*nm1;
        end if;
      else -- offset should be > 0, typically 2
        z := c(c'last-1-offset)/c(c'last-offset);
        if Is_Zero(c(c'last))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-1)/c(c'last))*nm1;
        end if;
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in OctoDobl_Complex_Vectors.Vector;
                    z : out OctoDobl_Complex_Numbers.Complex_Number;
                    e : out octo_double; fail : out boolean;
                    offset : in integer32 := 0 ) is

    use OctoDobl_Complex_Numbers;
    nm1 : constant double_float := double_float(c'last - 1);

  begin
    fail := Is_Zero(c(c'last-offset));
    if not fail then
      if offset = 0 then
        z := c(c'last-1)/c(c'last);
        if Is_Zero(c(c'last-1))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-2)/c(c'last-1))*nm1;
        end if;
      else -- offset should be > 0, typically 2
        z := c(c'last-1-offset)/c(c'last-offset);
        if Is_Zero(c(c'last))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-1)/c(c'last))*nm1;
        end if;
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in DecaDobl_Complex_Vectors.Vector;
                    z : out DecaDobl_Complex_Numbers.Complex_Number;
                    e : out deca_double; fail : out boolean;
                    offset : in integer32 := 0 ) is

    use DecaDobl_Complex_Numbers;
    nm1 : constant double_float := double_float(c'last - 1);

  begin
    fail := Is_Zero(c(c'last-offset));
    if not fail then
      if offset = 0 then
        z := c(c'last-1)/c(c'last);
        if Is_Zero(c(c'last-1))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-2)/c(c'last-1))*nm1;
        end if;
      else -- offset should be > 0, typically 2
        z := c(c'last-1-offset)/c(c'last-offset);
        if Is_Zero(c(c'last))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-1)/c(c'last))*nm1;
        end if;
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in HexaDobl_Complex_Vectors.Vector;
                    z : out HexaDobl_Complex_Numbers.Complex_Number;
                    e : out hexa_double; fail : out boolean;
                    offset : in integer32 := 0 ) is

    use HexaDobl_Complex_Numbers;
    nm1 : constant double_float := double_float(c'last - 1);

  begin
    fail := Is_Zero(c(c'last-offset));
    if not fail then
      if offset = 0 then
        z := c(c'last-1)/c(c'last);
        if Is_Zero(c(c'last-1))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-2)/c(c'last-1))*nm1;
        end if;
      else -- offset should be > 0, typically 2
        z := c(c'last-1-offset)/c(c'last-offset);
        if Is_Zero(c(c'last))
         then e := create(1.0);
         else e := AbsVal(z - c(c'last-1)/c(c'last))*nm1;
        end if;
      end if;
    end if;
  end Fabry;

  procedure Fabry ( c : in Standard_Complex_VecVecs.VecVec;
                    z : out Standard_Complex_Numbers.Complex_Number;
                    r : out double_float;
                    e : out double_float; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use Standard_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : double_float;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
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
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use DoblDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : double_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
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

  procedure Fabry ( c : in TripDobl_Complex_VecVecs.VecVec;
                    z : out TripDobl_Complex_Numbers.Complex_Number;
                    r : out triple_double;
                    e : out triple_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use TripDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : triple_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail
         then put_line("zero last coefficient");
         else put(zk); put("  error estimate : "); put(ek,3); new_line;
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := TripDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := TripDobl_Complex_Numbers_Polar.Radius(zk);
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
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use QuadDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : quad_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
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

  procedure Fabry ( c : in PentDobl_Complex_VecVecs.VecVec;
                    z : out PentDobl_Complex_Numbers.Complex_Number;
                    r : out penta_double;
                    e : out penta_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use PentDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : penta_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail
         then put_line("zero last coefficient");
         else put(zk); put("  error estimate : "); put(ek,3); new_line;
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := PentDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := PentDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( c : in OctoDobl_Complex_VecVecs.VecVec;
                    z : out OctoDobl_Complex_Numbers.Complex_Number;
                    r : out octo_double;
                    e : out octo_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use OctoDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : octo_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail
         then put_line("zero last coefficient");
         else put(zk); put("  error estimate : "); put(ek,3); new_line;
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := OctoDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := OctoDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( c : in DecaDobl_Complex_VecVecs.VecVec;
                    z : out DecaDobl_Complex_Numbers.Complex_Number;
                    r : out deca_double;
                    e : out deca_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use DecaDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : deca_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail
         then put_line("zero last coefficient");
         else put(zk); put("  error estimate : "); put(ek,3); new_line;
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := DecaDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := DecaDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( c : in HexaDobl_Complex_VecVecs.VecVec;
                    z : out HexaDobl_Complex_Numbers.Complex_Number;
                    r : out hexa_double;
                    e : out hexa_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use HexaDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : hexa_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail
         then put_line("zero last coefficient");
         else put(zk); put("  error estimate : "); put(ek,3); new_line;
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := HexaDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := HexaDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( file : in file_type;
                    c : in Standard_Complex_VecVecs.VecVec;
                    z : out Standard_Complex_Numbers.Complex_Number;
                    r : out double_float;
                    e : out double_float; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use Standard_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : double_float;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail then
          put_line(file,"zero last coefficient");
        else
          put(file,zk); put(file,"  error estimate :");
          put(file,ek,3); new_line(file);
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

  procedure Fabry ( file : in file_type;
                    c : in DoblDobl_Complex_VecVecs.VecVec;
                    z : out DoblDobl_Complex_Numbers.Complex_Number;
                    r : out double_double;
                    e : out double_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use DoblDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : double_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail then
          put_line(file,"zero last coefficient");
        else
          put(file,zk); put(file,"  error estimate : ");
          put(file,ek,3); new_line(file);
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

  procedure Fabry ( file : in file_type;
                    c : in TripDobl_Complex_VecVecs.VecVec;
                    z : out TripDobl_Complex_Numbers.Complex_Number;
                    r : out triple_double;
                    e : out triple_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use TripDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : triple_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail then
          put_line(file,"zero last coefficient");
        else
          put(file,zk); put(file,"  error estimate : ");
          put(file,ek,3); new_line(file);
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := TripDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := TripDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( file : in file_type;
                    c : in QuadDobl_Complex_VecVecs.VecVec;
                    z : out QuadDobl_Complex_Numbers.Complex_Number;
                    r : out quad_double;
                    e : out quad_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use QuadDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : quad_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail then
          put_line(file,"zero last coefficient");
        else
          put(file,zk); put(file,"  error estimate : ");
          put(file,ek,3); new_line(file);
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

  procedure Fabry ( file : in file_type;
                    c : in PentDobl_Complex_VecVecs.VecVec;
                    z : out PentDobl_Complex_Numbers.Complex_Number;
                    r : out penta_double;
                    e : out penta_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use PentDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : penta_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail then
          put_line(file,"zero last coefficient");
        else
          put(file,zk); put(file,"  error estimate : ");
          put(file,ek,3); new_line(file);
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := PentDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := PentDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( file : in file_type;
                    c : in OctoDobl_Complex_VecVecs.VecVec;
                    z : out OctoDobl_Complex_Numbers.Complex_Number;
                    r : out octo_double;
                    e : out octo_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use OctoDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : octo_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail then
          put_line(file,"zero last coefficient");
        else
          put(file,zk); put(file,"  error estimate : ");
          put(file,ek,3); new_line(file);
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := OctoDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := OctoDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( file : in file_type;
                    c : in DecaDobl_Complex_VecVecs.VecVec;
                    z : out DecaDobl_Complex_Numbers.Complex_Number;
                    r : out deca_double;
                    e : out deca_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use DecaDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : deca_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail then
          put_line(file,"zero last coefficient");
        else
          put(file,zk); put(file,"  error estimate : ");
          put(file,ek,3); new_line(file);
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := DecaDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := DecaDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

  procedure Fabry ( file : in file_type;
                    c : in HexaDobl_Complex_VecVecs.VecVec;
                    z : out HexaDobl_Complex_Numbers.Complex_Number;
                    r : out hexa_double;
                    e : out hexa_double; fail : out boolean;
                    offset : in integer32 := 0;
                    verbose : in boolean := true ) is

    use HexaDobl_Complex_Numbers;

    zk : Complex_Number;
    ek,rad : hexa_double;
    kfail : boolean;

  begin
    fail := true;
    for k in c'range loop
      Fabry(c(k).all,zk,ek,kfail,offset);
      if verbose then
        if kfail then
          put_line(file,"zero last coefficient");
        else
          put(file,zk); put(file,"  error estimate : ");
          put(file,ek,3); new_line(file);
        end if;
      end if;
      if not kfail then
        if k = c'first then
          z := zk; e := ek;
          r := HexaDobl_Complex_Numbers_Polar.Radius(z);
        else
          rad := HexaDobl_Complex_Numbers_Polar.Radius(zk);
          if rad < r
           then z := zk; e := ek; r := rad;
          end if;
        end if;
        fail := false;
      end if;
    end loop;
  end Fabry;

-- WRAPPERS :

  procedure Apply_Fabry
              ( c : in Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,e : double_float;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( file : in file_type;
                c : in Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,e : double_float;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(file,c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,e : double_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( file : in file_type;
                c : in DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,e : double_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(file,c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in TripDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : TripDobl_Complex_Numbers.Complex_Number;
    r,e : triple_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( file : in file_type;
                c : in TripDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : TripDobl_Complex_Numbers.Complex_Number;
    r,e : triple_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(file,c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : QuadDobl_Complex_Numbers.Complex_Number;
    r,e : quad_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( file : in file_type;
                c : in QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : QuadDobl_Complex_Numbers.Complex_Number;
    r,e : quad_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(file,c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in PentDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : PentDobl_Complex_Numbers.Complex_Number;
    r,e : penta_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( file : in file_type;
                c : in PentDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : PentDobl_Complex_Numbers.Complex_Number;
    r,e : penta_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(file,c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in OctoDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : OctoDobl_Complex_Numbers.Complex_Number;
    r,e : octo_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( file : in file_type;
                c : in OctoDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : OctoDobl_Complex_Numbers.Complex_Number;
    r,e : octo_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(file,c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in DecaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : DecaDobl_Complex_Numbers.Complex_Number;
    r,e : deca_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( file : in file_type;
                c : in DecaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : DecaDobl_Complex_Numbers.Complex_Number;
    r,e : deca_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(file,c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in HexaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : HexaDobl_Complex_Numbers.Complex_Number;
    r,e : hexa_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( file : in file_type;
                c : in HexaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : HexaDobl_Complex_Numbers.Complex_Number;
    r,e : hexa_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(file,c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

end Convergence_Radius_Estimates;
