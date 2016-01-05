with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
--with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;

package body DoblDobl_Winding_Numbers is

  function Consecutive_Differences ( logx : in Vector ) return Vector is

    res : Vector(logx'first..logx'last-1);

  begin
    for i in res'range loop
      res(i) := logx(i) - logx(i+1);
    end loop;
    return res;
  end Consecutive_Differences;

  function Consecutive_Errors ( difs : in Vector ) return Vector is

    res : Vector(difs'first..difs'last-1);

  begin
    for i in res'range loop
      res(i) := abs(difs(i) - difs(i+1));
    end loop;
    return res;
  end Consecutive_Errors;

  procedure Write_Extrapolation_Errors
              ( file : in file_type; e : in Vector;
                log10h : in double_double; ew : in integer32 ) is

    w,error : double_double;
    d_ew : constant double_float := create(ew);

  begin
    for i in e'range loop
      w := log10h/e(i);
      error := abs(w - d_ew);
      put(file,"estimate : "); put(file,w);
      put(file,"  error : "); put(file,error,3);
      new_line(file);
    end loop;
  end Write_Extrapolation_Errors;

-- DRIVER PROCEDURES :

  procedure Extrapolate_on_Errors_full
              ( logx : in Vector; h,log10h : in double_double;
                csc_dif,csc_err,ext_err : out Vector;
                ew : out integer32; error : out double_double ) is

    d,w,dd_k : double_double;
    d_ew,d_k : double_float;

  begin
    csc_dif := Consecutive_Errors(logx);
    csc_err := Consecutive_Errors(csc_dif);
    for i in ext_err'range loop
      ext_err(i) := LOG10(csc_err(i+1)) - LOG10(csc_err(i));
      w := log10h/ext_err(i);
    end loop;
    for k in 1..ext_err'last loop
      d_k := create(k); 
      dd_k := create(d_k)/w;
      d := h**(dd_k) - 1.0;
      for i in ext_err'first..ext_err'last-k loop
        ext_err(i) := ext_err(i) + (ext_err(i)-ext_err(i+1))/d;
        w := log10h/ext_err(i);
      end loop;
    end loop;
    ew := integer32(hi_part(w));
    d_ew := double_float(ew);
    error := abs(w-d_ew);
  end Extrapolate_on_Errors_full;

  procedure Extrapolate_on_Errors_full
              ( file : in file_type;
                logx : in Vector; h,log10h : in double_double;
                csc_dif,csc_err,ext_err : out Vector;
                ew : out integer32; error : out double_double ) is

    d,w,dd_k : double_double;
    d_ew,d_k : double_float;

  begin
    csc_dif := Consecutive_Errors(logx);
    csc_err := Consecutive_Errors(csc_dif);
    for i in ext_err'range loop
      ext_err(i) := LOG10(csc_err(i+1)) - LOG10(csc_err(i));
      w := log10h/ext_err(i);
      put(file,"estimate ("); put(file,i,1); put(file,",0) : ");
      put(file,w); new_line(file);
    end loop;
    ew := integer32(hi_part(w));
    for k in 1..ext_err'last loop
      d_k := create(k); d_ew := create(ew);
      dd_k := create(d_k)/create(d_ew);
      d := h**(dd_k) - 1.0;
      for i in ext_err'first..ext_err'last-k loop
        ext_err(i) := ext_err(i) + (ext_err(i)-ext_err(i+1))/d;
        w := log10h/ext_err(i);
        put(file,"estimate ("); put(file,i,1); put(file,",");
        put(file,k,1); put(file,") : ");
        put(file,w); new_line(file);
      end loop;
      ew := integer32(hi_part(w));
    end loop;
    d_ew := double_float(ew);
    error := abs(w-d_ew);
    Write_Extrapolation_Errors(file,ext_err,log10h,ew);
  end Extrapolate_on_Errors_full;

--  procedure Extrapolate_on_Errors_pipe
--              ( prev_logx,new_logx,h,log10h : in double_double;
--                csc_dif,csc_err,ext_err : in out Vector;
--                ew : out integer; error : out double_double ) is
--  begin
--    for i in csc_dif'first+1..csc_dif'last loop
--      csc_dif(i-1) := csc_dif(i);
--    end loop;
--    csc_dif(csc_dif'last) := prev_logx - new_logx;
--    for i in csc_err'first+1..csc_err'last loop
--      csc_err(i-1) := csc_err(i);
--    end loop;
--    csc_err(csc_err'last) := abs(csc_dif(csc_dif'last-1)
--                               - csc_dif(csc_dif'last));
--    for i in ext_err'first+1..ext_err'last loop
--      ext_err(i-1) := ext_err(i);
--    end loop;
--    ext_err(ext_err'last) := LOG10(csc_err(csc_err'last))
--                           - LOG10(csc_err(csc_err'last-1));
--  end Extrapolate_on_Errors_pipe;

--  procedure Extrapolate_on_Errors_pipe
--              ( file : in file_type;
--                logx : in Vector; h,log10h : in double_double;
--                csc_dif,csc_err,ext_err : in out Vector;
--                ew : out integer; error : out double_double ) is
--  begin
--    null;
--  end Extrapolate_on_Errors_pipe;

end DoblDobl_Winding_Numbers;
