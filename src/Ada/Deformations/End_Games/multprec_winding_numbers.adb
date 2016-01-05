with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;

package body Multprec_Winding_Numbers is

  function Consecutive_Errors ( logx : in Vector ) return Vector is

    dif : Multprec_Floating_Vectors.Vector(logx'first..logx'last-1);
    res : Multprec_Floating_Vectors.Vector(logx'first..logx'last-2);
    acc : Floating_Number;

  begin
    for i in dif'range loop             -- consecutive differences
      dif(i) := logx(i) - logx(i+1);
    end loop;
    for i in res'range loop
      acc := dif(i) - dif(i+1);
      res(i) := AbsVal(acc);
      Clear(acc);
    end loop;
    return res;
  end Consecutive_Errors;

  procedure Write_Extrapolation_Errors
              ( file : in file_type; e : in Vector;
                log10h : in Floating_Number; ew : in integer32 ) is

    w,ewf,error : Floating_Number;

  begin
    ewf := Create(ew);
    for i in e'range loop
      w := log10h/e(i);
      error := ewf - w;
      if error < 0.0 then Min(error); end if;
      put(file,"estimate : "); put(file,w);
      put(file,"  error : "); put(file,error,3);
      new_line(file);
    end loop;
  end Write_Extrapolation_Errors;

-- DRIVER PROCEDURES :

  procedure Extrapolate_on_Errors
                 ( logx : in Vector; h : in Floating_Number;
                   ew : out integer32; error : out Floating_Number ) is

    err : Vector(logx'first..logx'last-2) := Consecutive_Errors(logx);
    e : Vector(err'first..err'last-1);
    log10h : Floating_Number := LOG10(h);
    acc,d,w,q,mw : Floating_Number;
    fw : double_float;

  begin
    for i in e'range loop
      e(i) := LOG10(err(i+1));
      acc := LOG10(err(i));
      Sub(e(i),acc);
      Clear(acc);
      Clear(w);
      w := log10h/e(i);
    end loop;
    for k in 1..e'last loop
      fw := Round(w);
      ew := integer32(fw);
      mw := Create(ew);
      d := Create(integer(1));
      q := Create(k);
      Div(q,mw);
      Clear(mw);
      acc := h**q;
      Sub(acc,d);
      Clear(d);
      d := acc;
      Clear(q);
      for i in e'first..e'last-k loop
        acc := e(i)-e(i+1);
        Div(acc,d);
        Add(e(i),acc);
        Clear(acc);
        Clear(w);
        w := log10h/e(i);
      end loop;
      Clear(d);
    end loop;
    fw := Round(w);
    ew := integer32(fw);
    error := Create(ew); 
    Sub(error,w);
    if error < 0.0 then Min(error); end if;
    Clear(err); Clear(e); Clear(log10h); Clear(w);
  end Extrapolate_on_Errors;

  procedure Extrapolate_on_Errors
                 ( file : in file_type;
                   logx : in Vector; h : in Floating_Number;
                   ew : out integer32; error : out Floating_Number ) is

    err : Vector(logx'first..logx'last-2) := Consecutive_Errors(logx);
    e : Vector(err'first..err'last-1);
    log10h : Floating_Number := LOG10(h);
    acc,d,w,q : Floating_Number;
    fw : double_float;

  begin
    for i in e'range loop
      e(i) := LOG10(err(i+1));
      acc := LOG10(err(i));
      Sub(e(i),acc);
      Clear(acc);
      Clear(w);
      w := log10h/e(i);
      put(file,"estimate ("); put(file,i,1); put(file,",0) : ");
      put(file,w); new_line(file);
    end loop;
    for k in 1..e'last loop
      d := Create(integer(1));
      q := Create(k);
      Div(q,w);
      acc := h**q;
      Sub(acc,d);
      Clear(d);
      d := acc;
      Clear(q);
      for i in e'first..e'last-k loop
        acc := e(i)-e(i+1);
        Div(acc,d);
        Add(e(i),acc);
        Clear(acc);
        Clear(w);
        w := log10h/e(i);
        put(file,"estimate ("); put(file,i,1); put(file,",");
        put(file,k,1); put(file,") : ");
        put(file,w); new_line(file);
      end loop;
      Clear(d);
    end loop;
    fw := Round(w);
    ew := integer32(fw);
    error := Create(ew); 
    Sub(error,w);
    if error < 0.0 then Min(error); end if;
    Write_Extrapolation_Errors(file,e,log10h,ew);
    Clear(err); Clear(e); Clear(log10h); Clear(w);
  end Extrapolate_on_Errors;

end Multprec_Winding_Numbers;
