with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Random_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Random_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Double_Double_Vectors;
with Double_Double_Vectors_io;           use Double_Double_Vectors_io;
with Double_Double_Matrices;
with Quad_Double_Vectors;
with Quad_Double_Vectors_io;             use Quad_Double_Vectors_io;
with Quad_Double_Matrices;
with Multprec_Floating_Vectors;
with Multprec_Floating_Vectors_io;       use Multprec_Floating_Vectors_io;
with Multprec_Floating_Matrices;
with Standard_vLpRs_Tables;
with DoblDobl_vLpRs_Tables;
with QuadDobl_vLpRs_Tables;
with Standard_vLpRs_Algorithm;
with DoblDobl_vLpRs_Algorithm;
with QuadDobl_vLpRs_Algorithm;
with Multprec_vLpRs_Tables;
with Multprec_vLpRs_Algorithm;
with Standard_Winding_Numbers;
with DoblDobl_Winding_Numbers;
with QuadDobl_Winding_Numbers;
with Multprec_Winding_Numbers;

procedure ts_vlprs is

-- DESCRIPTION :
--   Test on the extrapolation algorithm on power series.

-- GENERATION AND EVALUATION OF A POWER SERIES :

  function Ask_Power_Series_Exponents
               ( n : integer32 ) return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Asks for a vector of n exponents.

    res : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0);

  begin
    put("Give the leading exponent ( -1, 0, or 1 ) : ");
    get(res(1));
    for i in 2..res'last loop
      res(i) := res(i-1) + 1;
    end loop;
    return res;
  end Ask_Power_Series_Exponents;

  function Standard_Random_Power_Series_Coefficients
               ( n : integer32 ) return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..n of randomly generated floating-point
  --   numbers, which represent the coefficients of a power series.

    res : Standard_Floating_Vectors.Vector(1..n);

  begin
    for i in res'range loop
      res(i) := Random;
    end loop;
    return res;
  end Standard_Random_Power_Series_Coefficients;

  function DoblDobl_Random_Power_Series_Coefficients
               ( n : integer32 ) return Double_Double_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..n of randomly generated floating-point
  --   numbers, which represent the coefficients of a power series.

    res : Double_Double_Vectors.Vector(1..n);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Random_Numbers.Random;
    end loop;
    return res;
  end DoblDobl_Random_Power_Series_Coefficients;

  function QuadDobl_Random_Power_Series_Coefficients
               ( n : integer32 ) return Quad_Double_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..n of randomly generated floating-point
  --   numbers, which represent the coefficients of a power series.

    res : Quad_Double_Vectors.Vector(1..n);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Random_Numbers.Random;
    end loop;
    return res;
  end QuadDobl_Random_Power_Series_Coefficients;

  function Multprec_Random_Power_Series_Coefficients
               ( n : integer32; s : natural32 )
               return Multprec_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..n of randomly generated floating-point
  --   numbers of size s, which represent the coefficients of a power series.

    res : Multprec_Floating_Vectors.Vector(1..n);

  begin
    for i in res'range loop
      res(i) := Random(s);
    end loop;
    return res;
  end Multprec_Random_Power_Series_Coefficients;

  function Standard_Eval_Power_Series
                ( s : double_float; a : Standard_Floating_Vectors.Vector;
                  m : Standard_Integer_Vectors.Vector; w : natural32 )
                return double_float is

  -- DESCRIPTION :
  --   Returns the sum of a_i*s**(m(i)/w).

    res : double_float := 0.0;

  begin
    for i in m'range loop
      res := res + a(i)*s**(double_float(m(i))/double_float(w));
    end loop;
    return res;
  end Standard_Eval_Power_Series;

  function DoblDobl_Eval_Power_Series
                ( s : double_double; a : Double_Double_Vectors.Vector;
                  m : Standard_Integer_Vectors.Vector; w : natural32 )
                return double_double is

  -- DESCRIPTION :
  --   Returns the sum of a_i*s**(m(i)/w).

    res : double_double := create(0.0);
    dd_m : double_double;
    dd_w : constant double_double := create(double_float(w));

  begin
    for i in m'range loop
      dd_m := Double_Double_Numbers.create(double_float(m(i)));
      res := res + a(i)*s**(dd_m/dd_w);
    end loop;
    return res;
  end DoblDobl_Eval_Power_Series;

  function QuadDobl_Eval_Power_Series
                ( s : quad_double; a : Quad_Double_Vectors.Vector;
                  m : Standard_Integer_Vectors.Vector; w : natural32 )
                return quad_double is

  -- DESCRIPTION :
  --   Returns the sum of a_i*s**(m(i)/w).

    res : quad_double := create(0.0);
    qd_m : quad_double;
    qd_w : constant quad_double := create(double_float(w));

  begin
    for i in m'range loop
      qd_m := Quad_Double_Numbers.create(double_float(m(i)));
      res := res + a(i)*s**(qd_m/qd_w);
    end loop;
    return res;
  end QuadDobl_Eval_Power_Series;

  function Multprec_Eval_Power_Series
                ( s : Floating_Number;
                  a : Multprec_Floating_Vectors.Vector;
                  m : Standard_Integer_Vectors.Vector; w : natural32 )
                return Floating_Number is

  -- DESCRIPTION :
  --   Returns the sum of a_i*s**(m(i)/w).

    res : Floating_Number := Create(integer(0));
    acc : Floating_Number;

  begin
    for i in m'range loop
      acc := s**(double_float(m(i))/double_float(w));
      Mul(acc,a(i));
      Add(res,acc);
      Clear(acc);
    end loop;
    return res;
  end Multprec_Eval_Power_Series;

 -- procedure Standard_Eval_Log_Power_Series
 --              ( s : in double_float; a : in Standard_Floating_Vectors.Vector;
 --                m : in Standard_Integer_Vectors.Vector;
 --                w : in natural32; logs,logx : out double_float ) is

  -- DESCRIPTION :
  --   On return are the logarithms of s and the result of the power series.

 -- begin
 --   logs := LOG10(s);
 --   logx := LOG10(abs(Standard_Eval_Power_Series(s,a,m,w)));
 -- end Standard_Eval_Log_Power_Series;

 -- procedure DoblDobl_Eval_Log_Power_Series
 --              ( s : in double_double; a : in Double_Double_Vectors.Vector;
 --                m : in Standard_Integer_Vectors.Vector;
 --                w : in natural32; logs,logx : out double_double ) is

  -- DESCRIPTION :
  --   On return are the logarithms of s and the result of the power series.

 -- begin
 --   logs := LOG10(s);
 --   logx := LOG10(abs(DoblDobl_Eval_Power_Series(s,a,m,w)));
 -- end DoblDobl_Eval_Log_Power_Series;

 -- procedure QuadDobl_Eval_Log_Power_Series
 --              ( s : in quad_double; a : in Quad_Double_Vectors.Vector;
 --                m : in Standard_Integer_Vectors.Vector;
 --                w : in natural32; logs,logx : out quad_double ) is

  -- DESCRIPTION :
  --   On return are the logarithms of s and the result of the power series.

 -- begin
 --   logs := LOG10(s);
 --   logx := LOG10(abs(QuadDobl_Eval_Power_Series(s,a,m,w)));
 -- end QuadDobl_Eval_Log_Power_Series;

 -- procedure Multprec_Eval_Log_Power_Series
 --              ( s : in Floating_Number;
 --                a : in Multprec_Floating_Vectors.Vector;
 --                m : in Standard_Integer_Vectors.Vector;
 --                w : in natural32; logs,logx : out Floating_Number ) is

  -- DESCRIPTION :
  --   On return are the logarithms of s and the result of the power series.

 --   acc : Floating_Number := Multprec_Eval_Power_Series(s,a,m,w);
 --   abs_acc : Floating_Number := AbsVal(acc);

 -- begin
 --   logs := LOG10(s);
 --   logx := LOG10(abs_acc);
 --   Clear(acc); Clear(abs_acc);
 -- end Multprec_Eval_Log_Power_Series;

  procedure Standard_Eval_Log_Power_Series
               ( s,a : in Standard_Floating_Vectors.Vector;
                 m : in Standard_Integer_Vectors.Vector; w : in natural32;
                 logs,logx : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   On return are the logarithms of s and the result of the power series.

  begin
    for i in s'range loop
      logs(i) := LOG10(s(i));
      logx(i) := LOG10(abs(Standard_Eval_Power_Series(s(i),a,m,w)));
    end loop;
  end Standard_Eval_Log_Power_Series;

  procedure DoblDobl_Eval_Log_Power_Series
               ( s,a : in Double_Double_Vectors.Vector;
                 m : in Standard_Integer_Vectors.Vector; w : in natural32;
                 logs,logx : out Double_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   On return are the logarithms of s and the result of the power series.

  begin
    for i in s'range loop
      logs(i) := LOG10(s(i));
      logx(i) := LOG10(abs(DoblDobl_Eval_Power_Series(s(i),a,m,w)));
    end loop;
  end DoblDobl_Eval_Log_Power_Series;

  procedure QuadDobl_Eval_Log_Power_Series
               ( s,a : in Quad_Double_Vectors.Vector;
                 m : in Standard_Integer_Vectors.Vector; w : in natural32;
                 logs,logx : out Quad_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   On return are the logarithms of s and the result of the power series.

  begin
    for i in s'range loop
      logs(i) := LOG10(s(i));
      logx(i) := LOG10(abs(QuadDobl_Eval_Power_Series(s(i),a,m,w)));
    end loop;
  end QuadDobl_Eval_Log_Power_Series;

  procedure Multprec_Eval_Log_Power_Series
               ( s,a : in Multprec_Floating_Vectors.Vector;
                 m : in Standard_Integer_Vectors.Vector; w : in natural32;
                 logs,logx : out Multprec_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   On return are the logarithms of s and the result of the power series.

    acc,abs_acc : Floating_Number;

  begin
    for i in s'range loop
      logs(i) := LOG10(s(i));
      acc := Multprec_Eval_Power_Series(s(i),a,m,w);
      abs_acc := AbsVal(acc);
      logx(i) := LOG10(abs_acc);
      Clear(acc); Clear(abs_acc);
    end loop;
  end Multprec_Eval_Log_Power_Series;

-- OUTPUT OF RESULTS :

  procedure Standard_Conclude
              ( estle,esterr : in double_float; 
                m : in integer32; winding : in natural32 ) is

  -- DESCRIPTION :
  --   Writes a concluding statement on the extrapolation for
  --   the winding number and leading exponent.

  -- ON ENTRY :
  --   estle     leading exponent computed by extrapolation;
  --   esterr    estimated error on the computed winding number;
  --   m         leading exponent of the power series;
  --   winding   exact winding number.

    tol_zero : constant double_float := 1.0E-3;
    r : constant double_float := double_float(m)/double_float(winding);

  begin
    put("exact leading exponent : "); put(r); new_line;
    put("the extrapolated value : "); put(estle); new_line;
    put("  with estimated error :");  put(esterr,3);
    if esterr < tol_zero
     then put_line("  okay.");
     else put_line("  at risk!!!");
    end if;
    if estle < -tol_zero then
      put_line("Series diverges to infinity.");
    elsif estle < tol_zero then
      put_line("Series converges to a nonzero constant.");
    else
      put_line("Series converges to zero.");
    end if;
  end Standard_Conclude;

  procedure DoblDobl_Conclude
              ( estle,esterr : in double_double; 
                m : in integer32; winding : in natural32 ) is

  -- DESCRIPTION :
  --   Writes a concluding statement on the extrapolation for
  --   the winding number and leading exponent.

  -- ON ENTRY :
  --   estle     leading exponent computed by extrapolation;
  --   esterr    estimated error on the computed winding number;
  --   m         leading exponent of the power series;
  --   winding   exact winding number.

    tol_zero : constant double_float := 1.0E-3;
    d_m : constant double_float := double_float(m);
    d_w : constant double_float := double_float(winding);
    r : constant double_double := create(d_m)/create(d_w);
 
  begin
    put("exact leading exponent : "); put(r); new_line;
    put("the extrapolated value : "); put(estle); new_line;
    put("  with estimated error : ");  put(esterr,3);
    if esterr < tol_zero
     then put_line("  okay.");
     else put_line("  at risk!!!");
    end if;
    if estle < -tol_zero then
      put_line("Series diverges to infinity.");
    elsif estle < tol_zero then
      put_line("Series converges to a nonzero constant.");
    else
      put_line("Series converges to zero.");
    end if;
  end DoblDobl_Conclude;

  procedure QuadDobl_Conclude
              ( estle,esterr : in quad_double; 
                m : in integer32; winding : in natural32 ) is

  -- DESCRIPTION :
  --   Writes a concluding statement on the extrapolation for
  --   the winding number and leading exponent.

  -- ON ENTRY :
  --   estle     leading exponent computed by extrapolation;
  --   esterr    estimated error on the computed winding number;
  --   m         leading exponent of the power series;
  --   winding   exact winding number.

    tol_zero : constant double_float := 1.0E-3;
    d_m : constant double_float := double_float(m);
    d_w : constant double_float := double_float(winding);
    qd_m : constant quad_double := create(d_m);
    qd_w : constant quad_double := create(d_w);
    r : constant quad_double := qd_m/qd_w;
 
  begin
    put("exact leading exponent : "); put(r); new_line;
    put("the extrapolated value : "); put(estle); new_line;
    put("  with estimated error : ");  put(esterr,3);
    if esterr < tol_zero
     then put_line("  okay.");
     else put_line("  at risk!!!");
    end if;
    if estle < -tol_zero then
      put_line("Series diverges to infinity.");
    elsif estle < tol_zero then
      put_line("Series converges to a nonzero constant.");
    else
      put_line("Series converges to zero.");
    end if;
  end QuadDobl_Conclude;

  procedure Multprec_Conclude
              ( estle,esterr : in Floating_Number; 
                m : in integer32; winding : in natural32 ) is

  -- DESCRIPTION :
  --   Writes a concluding statement on the extrapolation for
  --   the winding number and leading exponent.

  -- ON ENTRY :
  --   estle     leading exponent computed by extrapolation;
  --   esterr    estimated error on the computed winding number;
  --   m         leading exponent of the power series;
  --   winding   exact winding number.

    tol_zero : constant double_float := 1.0E-3;
    f_m : constant double_float := double_float(m);
    f_w : constant double_float := double_float(winding);
    r : constant Floating_Number := Create(f_m)/create(f_w);
 
  begin
    put("exact leading exponent : "); put(r); new_line;
    put("the extrapolated value : "); put(estle); new_line;
    put("  with estimated error :");  put(esterr,3);
    if esterr < tol_zero
     then put_line("  okay.");
     else put_line("  at risk!!!");
    end if;
    if estle < -tol_zero then
      put_line("Series diverges to infinity.");
    elsif estle < tol_zero then
      put_line("Series converges to a nonzero constant.");
    else
      put_line("Series converges to zero.");
    end if;
  end Multprec_Conclude;

  procedure Standard_Write
              ( l,v : in Standard_Floating_Vectors.Vector;
                m : in integer32; winding : in natural32 ) is

    w,err,prev_w,esterr : double_float;
   -- tol_zero : constant double_float := 1.0E-3;

  begin
    prev_w := 0.0;
    for i in v'range loop
      put(v(i),3); put(l(i),3); 
      w := v(i)/l(i); put(" "); put(w);
      esterr := abs(w-prev_w);
      err := abs(w - double_float(m)/double_float(winding));
      put(err,3); new_line;
      prev_w := w;
    end loop;
    Standard_Conclude(w,esterr,m,winding);
  end Standard_Write;

  procedure DoblDobl_Write
              ( l,v : in Double_Double_Vectors.Vector;
                m : in integer32; winding : in natural32 ) is

    w,err,prev_w,esterr : double_double;
   -- tol_zero : constant double_float := 1.0E-3;
    dd_m,dd_w,r : double_double;

  begin
    prev_w := create(0.0);
    for i in v'range loop
      put(v(i),3); put(" "); put(l(i),3); 
      w := v(i)/l(i); put(" "); put(w);
      esterr := abs(w-prev_w);
      dd_m := create(double_float(m));
      dd_w := create(double_float(winding));
      r := dd_m/dd_w;
      err := abs(w - r);
      put(" "); put(err,3); new_line;
      prev_w := w;
    end loop;
    DoblDobl_Conclude(w,esterr,m,winding);
  end DoblDobl_Write;

  procedure QuadDobl_Write
              ( l,v : in Quad_Double_Vectors.Vector;
                m : in integer32; winding : in natural32 ) is

    w,err,prev_w,esterr : quad_double;
   -- tol_zero : constant double_float := 1.0E-3;
    qd_m,qd_w,r : quad_double;

  begin
    prev_w := create(0.0);
    for i in v'range loop
      put(v(i),3); put(" "); put(l(i),3); 
      w := v(i)/l(i); put(" "); put(w);
      esterr := abs(w-prev_w);
      qd_m := create(double_float(m));
      qd_w := create(double_float(winding));
      r := qd_m/qd_w;
      err := abs(w - r);
      put(" "); put(err,3); new_line;
      prev_w := w;
    end loop;
    QuadDobl_Conclude(w,esterr,m,winding);
  end QuadDobl_Write;

  procedure Multprec_Write
              ( l,v : in Multprec_Floating_Vectors.Vector;
                m : in integer32; winding : in natural32 ) is

    w,err,esterr,prev_w : Floating_Number;
    fw : constant double_float := double_float(m)/double_float(winding);
    frw : Floating_Number := Create(fw);
   -- tol_zero : constant double_float := 1.0E-3;

  begin
    prev_w := Create(integer(0));
    for i in v'range loop
      put(v(i),3,3,3); put(l(i),3,3,3); 
      w := v(i)/l(i); put(" "); put(w);
      Clear(err); Clear(esterr);
      err := w - prev_w;
      esterr := AbsVal(err);
      Clear(err);
      err := AbsVal(w - frw);
      put(err,3,3,3); new_line;
      Copy(w,prev_w);
    end loop;
    Multprec_Conclude(w,esterr,m,winding);
    Clear(w); Clear(frw); Clear(prev_w); Clear(esterr);
  end Multprec_Write;

-- CALLING the EXTRAPOLATORS :

  procedure Standard_Winding_Number
              ( w : in natural32; ratio : in double_float;
                output : in boolean;
                s,logs : in out Standard_Floating_Vectors.Vector;
                logx : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies extrapolation to estimate the winding number.

  -- ON ENTRY :
  --   w        exact winding number;
  --   ratio    ratio used in decreasing sequence of s-values;
  --   output   flag to indicate whether extra output is wanted;
  --   s        geometrically decreasing sequence of s-values;
  --   logs     logarithms of the s-values.

  -- ON RETURN :
  --   s        s-values adjusted with correct winding number;
  --   logs     corresponding logarithms of adjusted s-values.

    use Standard_Floating_Vectors;
    use Standard_Winding_Numbers;

    log10h : constant double_float := LOG10(ratio);
    csc_dif : Vector(logx'first..logx'last-1);
    csc_err : Vector(logx'first..logx'last-2);
    ext_err : Vector(logx'first..logx'last-3);
    error : double_float;
    estwin : integer32;

  begin
    if output
     then Extrapolate_on_Errors_full
            (Standard_Output,logx,ratio,log10h,
             csc_dif,csc_err,ext_err,estwin,error);
     else Extrapolate_on_Errors_full
            (logx,ratio,log10h,csc_dif,csc_err,ext_err,estwin,error);
    end if;
    put("Estimate : "); put(estwin,1);
    put(" with error : "); put(error,3);
    if estwin = integer32(w)
     then put_line("  correct.");
     else put_line("  wrong!!!");
    end if;
    for i in s'range loop
      s(i) := s(i)**(1.0/double_float(w));
      logs(i) := logs(i)/double_float(w); 
    end loop;
  end Standard_Winding_Number;

  procedure DoblDobl_Winding_Number
              ( w : in natural32; ratio : in double_double;
                output : in boolean;
                s,logs : in out Double_Double_Vectors.Vector;
                logx : in Double_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies extrapolation to estimate the winding number.

  -- ON ENTRY :
  --   w        exact winding number;
  --   ratio    ratio used in decreasing sequence of s-values;
  --   output   flag to indicate whether extra output is wanted;
  --   s        geometrically decreasing sequence of s-values;
  --   logs     logarithms of the s-values.

  -- ON RETURN :
  --   s        s-values adjusted with correct winding number;
  --   logs     corresponding logarithms of adjusted s-values.

    use Double_Double_Vectors;
    use DoblDobl_Winding_Numbers;

    log10h : constant double_double := LOG10(ratio);
    csc_dif : Vector(logx'first..logx'last-1);
    csc_err : Vector(logx'first..logx'last-2);
    ext_err : Vector(logx'first..logx'last-3);
    error,dd_w : double_double;
    estwin : integer32;

  begin
    if output
     then Extrapolate_on_Errors_full
            (Standard_Output,logx,ratio,log10h,
             csc_dif,csc_err,ext_err,estwin,error);
     else Extrapolate_on_Errors_full
            (logx,ratio,log10h,csc_dif,csc_err,ext_err,estwin,error);
    end if;
    put("Estimate : "); put(estwin,1);
    put(" with error : "); put(error,3);
    if estwin = integer32(w)
     then put_line("  correct.");
     else put_line("  wrong!!!");
    end if;
    for i in s'range loop
      dd_w := create(double_float(w));
      s(i) := s(i)**(create(1.0)/dd_w);
      logs(i) := logs(i)/dd_w; 
    end loop;
  end DoblDobl_Winding_Number;

  procedure QuadDobl_Winding_Number
              ( w : in natural32; ratio : in quad_double;
                output : in boolean;
                s,logs : in out Quad_Double_Vectors.Vector;
                logx : in Quad_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies extrapolation to estimate the winding number.

  -- ON ENTRY :
  --   w        exact winding number;
  --   ratio    ratio used in decreasing sequence of s-values;
  --   output   flag to indicate whether extra output is wanted;
  --   s        geometrically decreasing sequence of s-values;
  --   logs     logarithms of the s-values.

  -- ON RETURN :
  --   s        s-values adjusted with correct winding number;
  --   logs     corresponding logarithms of adjusted s-values.

    use Quad_Double_Vectors;
    use QuadDobl_Winding_Numbers;

    log10h : constant quad_double := LOG10(ratio);
    csc_dif : Vector(logx'first..logx'last-1);
    csc_err : Vector(logx'first..logx'last-2);
    ext_err : Vector(logx'first..logx'last-3);
    error,qd_w : quad_double;
    estwin : integer32;

  begin
    if output
     then Extrapolate_on_Errors_full
            (Standard_Output,logx,ratio,log10h,
             csc_dif,csc_err,ext_err,estwin,error);
     else Extrapolate_on_Errors_full
            (logx,ratio,log10h,csc_dif,csc_err,ext_err,estwin,error);
    end if;
    put("Estimate : "); put(estwin,1);
    put(" with error : "); put(error,3);
    if estwin = integer32(w)
     then put_line("  correct.");
     else put_line("  wrong!!!");
    end if;
    for i in s'range loop
      qd_w := create(double_float(w));
      s(i) := s(i)**(Quad_Double_Numbers.Create(1.0)/qd_w);
      logs(i) := logs(i)/qd_w; 
    end loop;
  end QuadDobl_Winding_Number;
   
  procedure Multprec_Winding_Number 
              ( w : in natural32; ratio : in Floating_Number;
                output : in boolean;
                s,logs : in out Multprec_Floating_Vectors.Vector;
                logx : in Multprec_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies extrapolation to estimate the winding number.

  -- ON ENTRY :
  --   w        exact winding number;
  --   ratio    ratio used in decreasing sequence of s-values;
  --   output   flag to indicate whether extra output is wanted;
  --   s        geometrically decreasing sequence of s-values;
  --   logs     logarithms of the s-values.

  -- ON RETURN :
  --   s        s-values adjusted with correct winding number;
  --   logs     corresponding logarithms of adjusted s-values.

    error,acc : Floating_Number;
    estwin : integer32;

    use Multprec_Winding_Numbers;

  begin
    if output
     then Extrapolate_on_Errors(Standard_Output,logx,ratio,estwin,error);
     else Extrapolate_on_Errors(logx,ratio,estwin,error);
    end if;
    put("Estimate : "); put(estwin,1);
    put(" with error : "); put(error,3);
    if estwin = integer32(w)
     then put_line("  correct.");
     else put_line("  wrong!!!");
    end if;
    for i in s'range loop
      acc := s(i)**(1.0/double_float(w));
      Copy(acc,s(i));
      Clear(acc);
      Clear(logs(i));
      logs(i) := LOG10(s(i));
    end loop;
  end Multprec_Winding_Number;

  procedure Standard_Leading_Exponent
              ( m,r,le : in integer32;
                s,logs,logx : in Standard_Floating_Vectors.Vector ) is
                
  -- DESCRIPTION :
  --   Applies the vLpRs-algorithm to extrapolate the leading exponent
  --   after the correct winding number has adjusted the s-values.

  -- ON ENTRY :
  --   m        number of sample points;
  --   r        order of the extrapolator;
  --   le       leading exponent of the power series;
  --   s        s-values, geometrically decreasing sequence;
  --   logs     logarithms of the s-values;
  --   logx     logarithms of the absolute values of the samples.

    use Standard_Floating_Vectors;
    use Standard_Floating_Matrices;
    use Standard_vLpRs_Tables;
    use Standard_vLpRs_Algorithm;

    srp,dsp : Vector(1..r-1) := (1..r-1 => 0.0);
    p : Vector(0..r-1) := (0..r-1 => 0.0);
    l,v : Vector(0..r) := (0..r => 0.0);
    rt1,rt2 : Matrix(1..r-1,1..r-1);
    estle,err,prev_estle,esterr : double_float;

  begin
    put_line("Extrapolating using the correct winding number ... ");
    rt1(1,1) := 0.0; rt2(1,1) := 0.0;
   -- vlprs_full takes all sample points 
   -- vlprs_full(r,s,logs,logx,srp,dsp,p,l,v,rt1,rt2);
    vL_full(s(0..r),logs(0..r),logx(0..r),srp,dsp,p,l,v,rt1,rt2);
    put_line("The extrapolated values : ");
    Standard_Write(l,v,le,1);
    if m > r+1 then
      put_line("Continuing the extrapolation incrementally ... ");
      rt1 := rt2;
      prev_estle := v(r)/l(r);
      for i in r+1..m-1 loop
        vLpRs_pipe(s(i),logs(i),logx(i),srp,dsp,p,l,v,rt1,rt2);
        put(v(r),3,3,3); put(l(r),3,3,3);
        estle := v(r)/l(r); put(" "); put(estle);
        esterr := abs(estle - prev_estle);
        err := abs(estle - double_float(le));
        put(err,3,3,3); new_line;
        prev_estle := estle;
      end loop;
      Standard_Conclude(estle,esterr,le,1);
    end if;
  end Standard_Leading_Exponent;

  procedure DoblDobl_Leading_Exponent
              ( m,r,le : in integer32;
                s,logs,logx : in Double_Double_Vectors.Vector ) is
                
  -- DESCRIPTION :
  --   Applies the vLpRs-algorithm to extrapolate the leading exponent
  --   after the correct winding number has adjusted the s-values.

  -- ON ENTRY :
  --   m        number of sample points;
  --   r        order of the extrapolator;
  --   w        exact winding number;
  --   le       leading exponent of the power series;
  --   s        s-values, geometrically decreasing sequence;
  --   logs     logarithms of the s-values;
  --   logx     logarithms of the absolute values of the samples.

    use Double_Double_Vectors;
    use Double_Double_Matrices;
    use DoblDobl_vLpRs_Tables;
    use DoblDobl_vLpRs_Algorithm;

    srp,dsp : Vector(1..r-1) := (1..r-1 => create(0.0));
    p : Vector(0..r-1) := (0..r-1 => create(0.0));
    l,v : Vector(0..r) := (0..r => create(0.0));
    rt1,rt2 : Matrix(1..r-1,1..r-1);
    estle,err,prev_estle,esterr : double_double;

  begin
    put_line("Extrapolating using the correct winding number ... ");
    rt1(1,1) := create(0.0); rt2(1,1) := create(0.0);
   -- vlprs_full takes all sample points 
   -- vlprs_full(r,s,logs,logx,srp,dsp,p,l,v,rt1,rt2);
    vL_full(s(0..r),logs(0..r),logx(0..r),srp,dsp,p,l,v,rt1,rt2);
    put_line("The extrapolated values : ");
    DoblDobl_Write(l,v,le,1);
    if m > r+1 then
      put_line("Continuing the extrapolation incrementally ... ");
      rt1 := rt2;
      prev_estle := v(r)/l(r);
      for i in r+1..m-1 loop
        vLpRs_pipe(s(i),logs(i),logx(i),srp,dsp,p,l,v,rt1,rt2);
        put(v(r),3); put(" "); put(l(r),3);
        estle := v(r)/l(r); put(" "); put(estle);
        esterr := abs(estle - prev_estle);
        err := abs(estle - Double_Double_Numbers.create(double_float(le)));
        put(" "); put(err,3); new_line;
        prev_estle := estle;
      end loop;
      DoblDobl_Conclude(estle,esterr,le,1);
    end if;
  end DoblDobl_Leading_Exponent;

  procedure QuadDobl_Leading_Exponent
              ( m,r,le : in integer32;
                s,logs,logx : in Quad_Double_Vectors.Vector ) is
                
  -- DESCRIPTION :
  --   Applies the vLpRs-algorithm to extrapolate the leading exponent
  --   after the correct winding number has adjusted the s-values.

  -- ON ENTRY :
  --   m        number of sample points;
  --   r        order of the extrapolator;
  --   w        exact winding number;
  --   le       leading exponent of the power series;
  --   s        s-values, geometrically decreasing sequence;
  --   logs     logarithms of the s-values;
  --   logx     logarithms of the absolute values of the samples.

    use Quad_Double_Vectors;
    use Quad_Double_Matrices;
    use QuadDobl_vLpRs_Tables;
    use QuadDobl_vLpRs_Algorithm;

    srp,dsp : Vector(1..r-1) := (1..r-1 => create(0.0));
    p : Vector(0..r-1) := (0..r-1 => create(0.0));
    l,v : Vector(0..r) := (0..r => create(0.0));
    rt1,rt2 : Matrix(1..r-1,1..r-1);
    estle,err,prev_estle,esterr : quad_double;

  begin
    put_line("Extrapolating using the correct winding number ... ");
    rt1(1,1) := create(0.0); rt2(1,1) := create(0.0);
   -- vlprs_full takes all sample points 
   -- vlprs_full(r,s,logs,logx,srp,dsp,p,l,v,rt1,rt2);
    vL_full(s(0..r),logs(0..r),logx(0..r),srp,dsp,p,l,v,rt1,rt2);
    put_line("The extrapolated values : ");
    QuadDobl_Write(l,v,le,1);
    if m > r+1 then
      put_line("Continuing the extrapolation incrementally ... ");
      rt1 := rt2;
      prev_estle := v(r)/l(r);
      for i in r+1..m-1 loop
        vLpRs_pipe(s(i),logs(i),logx(i),srp,dsp,p,l,v,rt1,rt2);
        put(v(r),3); put(" "); put(l(r),3);
        estle := v(r)/l(r); put(" "); put(estle);
        esterr := abs(estle - prev_estle);
        err := abs(estle - Quad_Double_Numbers.create(double_float(le)));
        put(" "); put(err,3); new_line;
        prev_estle := estle;
      end loop;
      QuadDobl_Conclude(estle,esterr,le,1);
    end if;
  end QuadDobl_Leading_Exponent;

  procedure Multprec_Leading_Exponent
              ( m,r,le : in integer32;
                s,logs,logx : in Multprec_Floating_Vectors.Vector ) is
                
  -- DESCRIPTION :
  --   Applies the vLpRs-algorithm to extrapolate the leading exponent
  --   after the correct winding number has adjusted the s-values.

  -- ON ENTRY :
  --   m        number of sample points;
  --   r        order of the extrapolator;
  --   w        exact winding number;
  --   le       leading exponent of the power series;
  --   s        s-values, geometrically decreasing sequence;
  --   logs     logarithms of the s-values;
  --   logx     logarithms of the absolute values of the samples.

    use Multprec_Floating_Vectors;
    use Multprec_Floating_Matrices;
    use Multprec_vLpRs_Tables;
    use Multprec_vLpRs_Algorithm;

    srp,dsp : Vector(1..r-1) := (1..r-1 => Create(integer(0)));
    p : Vector(0..r-1) := (0..r-1 => Create(integer(0)));
    l,v : Vector(0..r) := (0..r => Create(integer(0)));
    rt1,rt2 : Matrix(1..r-1,1..r-1);
    prev_estle,estle,acc,err,esterr : Floating_Number;

  begin
    put_line("Extrapolating using the correct winding number ... ");
   -- vlprs_pipe(Standard_Output,r,s,logs,logx,srp,dsp,p,l,v,rt1,rt2);
    rt1(1,1) := Create(integer(0)); rt2(1,1) := Create(integer(0));
   -- vlprs_full takes all sample points
   -- vlprs_full(r,s,logs,logx,srp,dsp,p,l,v,rt1,rt2);
    vL_full(s(0..r),logs(0..r),logx(0..r),srp,dsp,p,l,v,rt1,rt2);
    put_line("The extrapolated values : ");
    Multprec_Write(l,v,le,1);
    if m > r+1 then
      put_line("Continuing the extrapolation incrementally ... ");
      Copy(rt2,rt1);
      prev_estle := v(r)/l(r);
      for i in r+1..m-1 loop
        vLpRs_pipe(s(i),logs(i),logx(i),srp,dsp,p,l,v,rt1,rt2);
        put(v(r),3,3,3); put(l(r),3,3,3);
        Clear(estle);
        estle := v(r)/l(r); put(" "); put(estle);
        acc := estle - prev_estle;
        esterr := AbsVal(acc); Clear(acc);
        acc := Create(le);
        Sub(acc,estle);
        err := AbsVal(acc); Clear(acc);
        put(err,3,3,3); new_line;
        Copy(estle,prev_estle);
      end loop;
      Multprec_Conclude(estle,esterr,le,1);
    end if;
  end Multprec_Leading_Exponent;

  procedure Standard_Extrapolate
              ( m,r : in integer32; w : in natural32;
                a : in Standard_Floating_Vectors.Vector;
                exp : in Standard_Integer_Vectors.Vector ) is

    s,logs,logx : Standard_Floating_Vectors.Vector(0..m-1);
    smax,ratio : double_float := 0.0;
    ans : character;
    output : boolean;

  begin
    put("Give largest s-value : "); get(smax);
    put("Give decreasing ratio (< 1.0) : "); get(ratio);
    s(0) := smax;
    for i in 1..s'last loop
      s(i) := ratio*s(i-1);
    end loop;
    put("Do you want extra output ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if output then
      put("The sequence of "); put(m,1);
      put_line(" s-values : "); put_line(s);
    end if;
    Standard_Eval_Log_Power_Series(s,a,exp,w,logs,logx);
    Standard_Winding_Number(w,ratio,output,s,logs,logx);
    Standard_Leading_Exponent(m,r,exp(exp'first),s,logs,logx);
  end Standard_Extrapolate;

  procedure DoblDobl_Extrapolate
              ( m,r : in integer32; w : in natural32;
                a : in Double_Double_Vectors.Vector;
                exp : in Standard_Integer_Vectors.Vector ) is

    s,logs,logx : Double_Double_Vectors.Vector(0..m-1);
    smax,ratio : double_double := create(0.0);
    ans : character;
    output : boolean;

  begin
    put("Give largest s-value : "); skip_line; get(smax);
    put("Give decreasing ratio (< 1.0) : "); get(ratio);
    s(0) := smax;
    for i in 1..s'last loop
      s(i) := ratio*s(i-1);
    end loop;
    put("Do you want extra output ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if output then
      put("The sequence of "); put(m,1);
      put_line(" s-values : "); put_line(s);
    end if;
    DoblDobl_Eval_Log_Power_Series(s,a,exp,w,logs,logx);
    DoblDobl_Winding_Number(w,ratio,output,s,logs,logx);
    DoblDobl_Leading_Exponent(m,r,exp(exp'first),s,logs,logx);
  end DoblDobl_Extrapolate;

  procedure QuadDobl_Extrapolate
              ( m,r : in integer32; w : in natural32;
                a : in Quad_Double_Vectors.Vector;
                exp : in Standard_Integer_Vectors.Vector ) is

    s,logs,logx : Quad_Double_Vectors.Vector(0..m-1);
    smax,ratio : quad_double := create(0.0);
    ans : character;
    output : boolean;

  begin
    put("Give largest s-value : "); skip_line; get(smax);
    put("Give decreasing ratio (< 1.0) : "); get(ratio);
    s(0) := smax;
    for i in 1..s'last loop
      s(i) := ratio*s(i-1);
    end loop;
    put("Do you want extra output ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if output then
      put("The sequence of "); put(m,1);
      put_line(" s-values : "); put_line(s);
    end if;
    QuadDobl_Eval_Log_Power_Series(s,a,exp,w,logs,logx);
    QuadDobl_Winding_Number(w,ratio,output,s,logs,logx);
    QuadDobl_Leading_Exponent(m,r,exp(exp'first),s,logs,logx);
  end QuadDobl_Extrapolate;

  procedure Multprec_Extrapolate
              ( m,r : in integer32; w,sz : in natural32;
                a : in Multprec_Floating_Vectors.Vector;
                exp : in Standard_Integer_Vectors.Vector ) is

    use Multprec_Floating_Vectors;

    s,logs,logx : Vector(0..m-1);
    smax,ratio : Floating_Number;
    ans : character;
    output : boolean;

  begin
    put("Give largest s-value : "); get(smax);
    put("Give decreasing ratio (< 1.0) : "); get(ratio);
    Set_Size(smax,sz);
    Set_Size(ratio,sz);
    s(0) := smax;
    for i in 1..s'last loop
      s(i) := ratio*s(i-1);
      Set_Size(s(i),sz);
    end loop;
    put("Do you want extra output ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if output then
      put("The sequence of "); put(m,1);
      put_line(" s-values : "); put_line(s);
    end if;
    Multprec_Eval_Log_Power_Series(s,a,exp,w,logs,logx);
    Multprec_Winding_Number(w,ratio,output,s,logs,logx);
    Multprec_Leading_Exponent(m,r,exp(exp'first),s,logs,logx);
  end Multprec_Extrapolate;

  procedure Standard_Test_vlprs ( n,m,r : in integer32; w : in natural32 ) is

  -- DESCRIPTION :
  --   Test on the extrapolation algorithm in standard arithmetic.

    exp : constant Standard_Integer_Vectors.Vector(1..n)
        := Ask_Power_Series_Exponents(n);
    pcf : constant Standard_Floating_Vectors.Vector(1..n)
        := Standard_Random_Power_Series_Coefficients(n);

  begin
    put("the exponents : "); put(exp); new_line;
    put_line("the coefficients : "); put_line(pcf);
    Standard_Extrapolate(m,r,w,pcf,exp);
  end Standard_Test_vlprs;

  procedure DoblDobl_Test_vlprs ( n,m,r : in integer32; w : in natural32 ) is

  -- DESCRIPTION :
  --   Test on the extrapolation algorithm in standard arithmetic.

    exp : constant Standard_Integer_Vectors.Vector(1..n)
        := Ask_Power_Series_Exponents(n);
    pcf : constant Double_Double_Vectors.Vector(1..n)
        := DoblDobl_Random_Power_Series_Coefficients(n);

  begin
    put("the exponents : "); put(exp); new_line;
    put_line("the coefficients : "); put_line(pcf);
    DoblDobl_Extrapolate(m,r,w,pcf,exp);
  end DoblDobl_Test_vlprs;

  procedure QuadDobl_Test_vlprs ( n,m,r : in integer32; w : in natural32 ) is

  -- DESCRIPTION :
  --   Test on the extrapolation algorithm in standard arithmetic.

    exp : constant Standard_Integer_Vectors.Vector(1..n)
        := Ask_Power_Series_Exponents(n);
    pcf : constant Quad_Double_Vectors.Vector(1..n)
        := QuadDobl_Random_Power_Series_Coefficients(n);

  begin
    put("the exponents : "); put(exp); new_line;
    put_line("the coefficients : "); put_line(pcf);
    QuadDobl_Extrapolate(m,r,w,pcf,exp);
  end QuadDobl_Test_vlprs;

  procedure Multprec_Test_vlprs ( n,m,r : in integer32; w,s : in natural32 ) is

  -- DESCRIPTION :
  --   Test on the extrapolation algorithm in multiprecision arithmetic.

    exp : constant Standard_Integer_Vectors.Vector(1..n)
        := Ask_Power_Series_Exponents(n);
    pcf : constant Multprec_Floating_Vectors.Vector(1..n)
        := Multprec_Random_Power_Series_Coefficients(n,s);

  begin
    put("the exponents : "); put(exp); new_line;
    put_line("the coefficients : "); put_line(pcf);
    Multprec_Extrapolate(m,r,w,s,pcf,exp);
  end Multprec_Test_vlprs;

  procedure Test_vlprs is

    n,m,r : integer32 := 0;
    w,d,s : natural32 := 0;
    ans : character;
 
  begin
    put("Give the order of the extrapolator : "); get(r);
    loop
      put("Give the number of points : "); get(m);
      exit when m > r;
      put("should be larger than "); put(r,1); put_line(".  Please retry.");
    end loop;
    put("Give the length of the power series : "); get(n);
    put("Give the winding number : "); get(w);
    new_line;
    put_line("MENU for aritmetic to extrapolate : ");
    put_line("  1. extrapolation in standard double arithmetic;");
    put_line("  2. extrapolation in double double arithmetic;");
    put_line("  3. extrapolation in quad double arithmetic;");
    put_line("  4. extrapolation in multiprecision arithmetic.");
    put("Type 1, 2, 3, or 4 to select arithmetic : ");
    Ask_Alternative(ans,"1234");
    new_line;
    case ans is
      when '1' => Standard_Test_vlprs(n,m,r,w);
      when '2' => DoblDobl_Test_vlprs(n,m,r,w);
      when '3' => QuadDobl_Test_vlprs(n,m,r,w);
      when others
        => put("Give the number of decimal places (<15 is standard) : ");
           get(d);
           if d < 15
            then Standard_Test_vlprs(n,m,r,w);
            else s := Decimal_to_Size(d);
                 Multprec_Test_vlprs(n,m,r,w,s);
          end if;
    end case;
  end Test_vlprs;

begin
  Test_vlprs;
end ts_vlprs;
