with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;
with Quad_Double_Vectors_io;             use Quad_Double_Vectors_io;
with Quad_Double_Vector_Norms;           use Quad_Double_Vector_Norms;
with Quad_Double_Matrices;               use Quad_Double_Matrices;
with Quad_Double_Matrices_io;            use Quad_Double_Matrices_io;
--with Standard_Integer_Vectors;
--with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with QuadDobl_vLpRs_Algorithm;           use QuadDobl_vLpRs_Algorithm;

package body Directions_of_QuadDobl_Paths is

-- FIRST ORDER EXTRAPOLATION  :

  procedure Affine_Update_Direction
                ( t,prev_t,target : in Complex_Number;
                  x,prevx : in QuadDobl_Complex_Vectors.Vector;
                  prevdls,prevstep : in out quad_double;
                  prevdiff,v : in out Quad_Double_Vectors.Vector ) is

    newdls : quad_double;
    newstep : constant quad_double := AbsVal(t-prev_t);
    newdiff : Quad_Double_Vectors.Vector(prevdiff'range);
    ratio,factor : quad_double;

  begin
    for i in v'range loop
      newdiff(i) := LOG10(AbsVal(x(i))) - LOG10(AbsVal(prevx(i)));
    end loop;
    newdls := LOG10(AbsVal(target-prev_t)) - LOG10(AbsVal(target-t));
    if not is_zero(prevdls) then
      ratio := prevstep/newstep;
      factor := prevdls - ratio*newdls;
      for i in v'range loop
        v(i) := (prevdiff(i) - ratio*newdiff(i))/factor;
      end loop;
    end if;
    prevdiff := newdiff;
    prevstep := newstep;
    prevdls := newdls;
  end Affine_Update_Direction;

  procedure Projective_Update_Direction
                ( t,prev_t,target : in Complex_Number;
                  x,prevx : in QuadDobl_Complex_Vectors.Vector;
                  prevdls,prevstep : in out quad_double;
                  prevdiff,v : in out Quad_Double_Vectors.Vector ) is

    newdls : quad_double;
    newstep : constant quad_double := AbsVal(t-prev_t);
    newdiff : Quad_Double_Vectors.Vector(prevdiff'range);
    ratio,factor : quad_double;

  begin
    for i in v'range loop
      newdiff(i) := ( LOG10(AbsVal(x(i))) - LOG10(AbsVal(x(x'last))) )
      - ( LOG10(AbsVal(prevx(i))) - LOG10(AbsVal(prevx(prevx'last))) );
    end loop;
    newdls := LOG10(AbsVal(target-prev_t)) - LOG10(AbsVal(target-t));
    if not is_zero(prevdls) then -- if prevdls /= 0.0 then
      ratio := prevstep/newstep;
      factor := prevdls - ratio*newdls;
      for i in v'range loop
        v(i) := (prevdiff(i) - ratio*newdiff(i))/factor;
      end loop;
    end if;
    prevdiff := newdiff;
    prevstep := newstep;
    prevdls := newdls;
  end Projective_Update_Direction;

-- DATA MANAGEMENT FOR HIGHER-ORDER EXTRAPOLATION :

  procedure Shift_Up ( v : in out Quad_Double_Vectors.Vector;
                       x : in quad_double ) is

  begin
    for i in reverse v'first..(v'last-1) loop
      v(i+1) := v(i);
    end loop;
    v(v'first) := x;
  end Shift_Up;

  procedure Extrapolation_Window
                ( r,m : in integer32; t,target : in Complex_Number;
                  x : in QuadDobl_Complex_Vectors.Vector;
                  dt,s,logs : in out Quad_Double_Vectors.Vector;
                  logx : in out Quad_Double_VecVecs.VecVec ) is

    use Quad_Double_Vectors;

  begin
    if (r = s'last) and (logx(r) /= null) then
      for i in s'first+1..s'last loop   -- shift the data
        s(i-1) := s(i);
        dt(i-1) := dt(i);
        logs(i-1) := logs(i);
        logx(i-1).all := logx(i).all;
      end loop;
    end if;
    dt(r) := (AbsVal(target-t));
    s(r) := (dt(r))**(1.0/Quad_Double_Numbers.create(m));
    logs(r) := LOG10(s(r));
  end Extrapolation_Window;

  procedure Refresh_Window
               ( r,m : in integer32;
                 dt : in Quad_Double_Vectors.Vector;
                 s,logs : in out Quad_Double_Vectors.Vector ) is
  begin
    for i in s'first..r loop
      s(i) := (dt(i))**(1.0/Quad_Double_Numbers.Create(m));
      logs(i) := LOG10(s(i));
    end loop;
  end Refresh_Window;

  procedure Write_Update_Information
                ( file : in file_type;
                  r,m : in integer32; s,logs : in quad_double;
                  logx : in Quad_Double_Vectors.Vector ) is
  begin
    put(file,"extrapolation with r = "); put(file,r,1);
    put(file," and m = "); put(file,m,1); --put_line(file," : ");
   -- put(file,"m : "); put(file,s); 
   -- put(file,"  log(s) : "); put(file,logs);
   -- put(file,"   log|x(s)| : ");
   -- put(file,logx);
   -- for i in logx'range loop
   --   put(file,logx(i));
   -- end loop;
    new_line(file);
  end Write_Update_Information;

  procedure Affine_Update_Extrapolation_Data
                ( r,m : in integer32; t,target : in Complex_Number;
                  x : in QuadDobl_Complex_Vectors.Vector;
                  dt,s,logs : in out Quad_Double_Vectors.Vector;
                  logx,wvl0,wvl1,wvltmp : in out Quad_Double_VecVecs.VecVec ) is

    use Quad_Double_Vectors;

  begin
    Extrapolation_Window(r,m,t,target,x,dt,s,logs,logx);
    if logx(r) = null
     then logx(r) := new Quad_Double_Vectors.Vector(x'range);
    end if;
    if r > 0 then
      if wvl0(r) = null then
         wvl0(r) := new Quad_Double_Vectors.Vector'(x'range => create(0.0));
         wvl1(r) := new Quad_Double_Vectors.Vector'(x'range => create(0.0));
         wvltmp(r) := new Quad_Double_Vectors.Vector'(x'range => create(0.0));
      end if;
    end if;
    for i in x'range loop
      logx(r)(i) := LOG10(AbsVal(x(i)));
    end loop;
  end Affine_Update_Extrapolation_Data;

  procedure Projective_Update_Extrapolation_Data
                ( r,m : in integer32; t,target : in Complex_Number;
                  x : in QuadDobl_Complex_Vectors.Vector;
                  dt,s,logs : in out Quad_Double_Vectors.Vector;
                  logx : in out Quad_Double_VecVecs.VecVec ) is

    use Quad_Double_Vectors;

  begin
    Extrapolation_Window(r,m,t,target,x,dt,s,logs,logx);
    if logx(r) = null
     then logx(r) := new Quad_Double_Vectors.Vector(x'first..x'last-1);
    end if;
    for i in x'first..x'last-1 loop
      logx(r)(i) := LOG10(AbsVal(x(i))) - LOG10(AbsVal(x(x'last)));
    end loop;
  end Projective_Update_Extrapolation_Data;

  procedure Update_Errors
              ( r : in integer32;
                errorv : in out Quad_Double_Vectors.Vector;
                error : out quad_double;
                wvl0,wvl1,wvltmp : in out Quad_Double_VecVecs.VecVec ) is
  begin
    if r = 1 then
      for i in errorv'range loop
        errorv(i) := abs(wvltmp(r)(i) - wvl1(r)(i));
      end loop;
    else
      for i in errorv'range loop
        errorv(i) := abs(wvltmp(r-1)(i) - wvl1(r-1)(i));
      end loop;
    end if;
    error := Sum_Norm(errorv);
    for i in 1..r loop
      wvl0(i).all := wvl1(i).all;
      wvl1(i).all := wvltmp(i).all;
    end loop;
  end Update_Errors;

-- ESTIMATING THE WINDING NUMBER :

  procedure Frequency_of_Estimate
               ( newest : in integer32; max : in natural32;
                 m,estm : in out integer32; cnt : in out natural32;
                 newm : out boolean ) is
  begin
    if cnt = 0 then                                       -- initial estimate
      estm := newest;
      cnt := 1;
    elsif newest = estm then         -- update frequency for current estimate
      cnt := cnt + 1;
    else 
      cnt := 1;                                         -- new estimate found
      estm := newest;
    end if;
    if estm /= m then        -- look for modification of current cycle number
      if (cnt >= max) --and (eps <= 0.1)
       then m := estm; newm := true;
       else newm := false;
      end if;
     else newm := false;
    end if;
  end Frequency_of_Estimate;

  procedure Extrapolate_on_Errors
               ( r : in integer32; h : in quad_double;
                 err : in Quad_Double_Vectors.Vector;
                 estm : out Quad_Double_Vectors.Vector ) is

    hm,exterr : Quad_Double_Vectors.Vector(1..r+1);
    dlog : constant quad_double := log10(h);
    f : quad_double;

  begin
    for j in exterr'range loop
      exterr(j) := log10(err(j-1)) - log10(err(j));
    end loop;
    estm(1) := dlog/exterr(1);                         -- 0th order estimate
   -- if (m(1) < 0.0001) or (m(1) > 1000.0)              -- avoid divergence
   --  then m(r+1) := m(1);
   --  else 
          hm(1) := h**(1.0/estm(1));
          for k in 1..r loop
            f := hm(k) - 1.0;
            for j in 1..r-k+1 loop
              exterr(j) := exterr(j+1) + (exterr(j+1) - exterr(j))/f;
            end loop;
            estm(k+1) := dlog/exterr(1);
   --  exit when ((m(k+1) < 0.0001) or (m(k+1) > 1000.0));
            hm(k+1) := h**(Quad_Double_Numbers.create(k+1)/estm(k+1));
          end loop;
   -- end if;
  exception
    when others => null;
  end Extrapolate_on_Errors;

  procedure Accuracy_of_Estimates
               ( estm : in Quad_Double_Vectors.Vector;
                 success : out boolean; k : out integer32;
                 estwin : out integer32; eps : out quad_double ) is

    res,wrk : integer32;
    best_eps,prev_eps : quad_double;

  begin
    k := estm'first-1;                -- take zero order as the best
    res := integer32(hihi_part(estm(estm'first)));
    eps := abs(estm(estm'first) - Quad_Double_Numbers.create(res));
    best_eps := eps;
    success := true;                  -- assume extrapolation worked
    for i in estm'first+1..estm'last loop
      wrk := integer32(hihi_part(estm(i)));
      eps := abs(estm(i) - Quad_Double_Numbers.create(wrk));
      for j in estm'first..i-1 loop   -- compare with previous estimates
        prev_eps := abs(estm(j) - Quad_Double_Numbers.create(wrk));
        if prev_eps > eps
         then success := false;       -- extrapolation failed
        end if;
        exit when (not success);
      end loop;
      exit when (not success);
      if eps < best_eps then
        res := wrk;                   -- update the result
        k := i-1;                     -- its order
        best_eps := eps;              -- and its accuracy
      end if;
    end loop;
    estwin := res;
    eps := best_eps;
  end Accuracy_of_Estimates;

  procedure Estimate0
               ( r : in integer32; max : in natural32;
                 m,estm : in out integer32; cnt : in out natural32;
                 dt : in Quad_Double_Vectors.Vector;
                 err,newerr : in quad_double; rat,eps : out quad_double;
                 newm : out boolean ) is

    dferr : constant quad_double := log10(newerr) - log10(err);
    h : constant quad_double := dt(r)/dt(r-1);
    dfsr1 : constant quad_double := log10(h);
    ratio : constant quad_double := abs(dfsr1/dferr);
    res : integer32 := integer32(hihi_part(ratio));
    accuracy : constant quad_double
             := abs(Quad_Double_Numbers.create(res) - ratio);

  begin
    if res = 0
     then res := 1;
    end if;
    Frequency_of_Estimate(res,max,m,estm,cnt,newm);
    rat := ratio; eps := accuracy;
  end Estimate0;

  procedure Estimate_Winding_Number
               ( r : in integer32;
                 max : in natural32; m,estm : in out integer32;
                 cnt : in out natural32; h : in quad_double;
                 diferr : in Quad_Double_Vectors.Vector;
                 rat,eps : out quad_double; newm : out boolean ) is

    res,order : integer32;
    accuracy : quad_double;
    esm : Quad_Double_Vectors.Vector(1..r+1);
    estgood : boolean;

  begin
    Extrapolate_on_Errors(r,h,diferr(0..r+1),esm);
    Accuracy_of_Estimates(esm,estgood,order,res,accuracy);
    if res <= 0
     then res := m; -- keep the current value for the winding number
     else Frequency_of_Estimate(res,max,m,estm,cnt,newm);
    end if;
    rat := esm(order+1);
    eps := accuracy;
  end Estimate_Winding_Number;

  procedure Estimate_Winding_Number
               ( file : in file_type; r : in integer32;
                 max : in natural32; m,estm : in out integer32;
                 cnt : in out natural32; h : in quad_double;
                 diferr : in Quad_Double_Vectors.Vector;
                 rat,eps : out quad_double; newm : out boolean ) is

    res,order : integer32;
    accuracy : quad_double;
    esm : Quad_Double_Vectors.Vector(1..r+1);
    estgood : boolean;

  begin
    Extrapolate_on_Errors(r,h,diferr(0..r+1),esm);
    put(file,"estm(0.."); put(file,r,1); put(file,") : ");
    for i in esm'range loop
      put(file," "); put(file,esm(i),3);
    end loop;
    Accuracy_of_Estimates(esm,estgood,order,res,accuracy);
    if res <= 0 then
      put_line(file,"  wrong result.");
      res := m; -- keep the current value for the winding number
    else
      if estgood then
        put_line(file,"  extrapolation succeeded.");
      else
        put(file,"  extrapolation failed, order = ");
        put(file,order,1); new_line(file);
      end if;
      Frequency_of_Estimate(res,max,m,estm,cnt,newm);
    end if;
    rat := esm(order+1);
    eps := accuracy;
  end Estimate_Winding_Number;

-- APPLYING THE vLpRs-Algorithm :

  function vLpRs_Extrapolate
                ( r : integer32; s,logs : Quad_Double_Vectors.Vector;
                  logx : Quad_Double_VecVecs.VecVec )
                return Quad_Double_Vectors.Vector is

    res : Quad_Double_Vectors.Vector(logx(r).all'range);
    logx1 : Quad_Double_Vectors.Vector(0..r);
    rt1,rt2 : Matrix(1..r-1,1..r-1);
    srp,dsp : Quad_Double_Vectors.Vector(1..r-1) := (0..r-1 => create(0.0));
    p : Quad_Double_Vectors.Vector(0..r-1) := (0..r-1 => create(0.0));
    l,v : Quad_Double_Vectors.Vector(0..r) := (0..r => create(1.0));

  begin
    rt1(1,1) := create(0.0); rt2(1,1) := create(0.0);
    for i in res'range loop
      for j in logx1'range loop
        logx1(j) := logx(j)(i);
      end loop;
      vLpRs_full(r,s(0..r),logs(0..r),logx1(0..r),srp,dsp,p,l,v,rt1,rt2);
      res(i) := v(r)/l(r);
    end loop;
    return res;
  end vLpRs_Extrapolate;

  function vLpRs_Extrapolate
                ( file : file_type; r : integer32;
                  s,logs : Quad_Double_Vectors.Vector;
                  logx : Quad_Double_VecVecs.VecVec )
                return Quad_Double_Vectors.Vector is

    res : Quad_Double_Vectors.Vector(logx(r).all'range);
    logx1 : Quad_Double_Vectors.Vector(0..r);
    rt1,rt2 : Matrix(1..r-1,1..r-1);
    srp,dsp : Quad_Double_Vectors.Vector(1..r-1) := (0..r-1 => create(0.0));
    p : Quad_Double_Vectors.Vector(0..r-1) := (0..r-1 => create(0.0));
    l,v : Quad_Double_Vectors.Vector(0..r) := (0..r => create(0.0));

  begin
    for i in rt1'range(1) loop
      for j in rt1'range(2) loop
        rt1(i,j) := create(0.0); rt2(i,j) := create(0.0);
      end loop;
    end loop;
    for i in res'range loop
      for j in logx1'range loop
        logx1(j) := logx(j)(i);
      end loop;
      vLpRs_pipe(file,r,s(0..r),logs(0..r),logx1(0..r),srp,dsp,p,l,v,rt1,rt2);
      res(i) := v(r)/l(r);
    end loop;
    return res;
  end vLpRs_Extrapolate;

  procedure vLpRs_Extrapolate
                ( r : in integer32;
                  s,logs : in Quad_Double_Vectors.Vector;
                  logx,wvl0 : in Quad_Double_VecVecs.VecVec;
                  wvl1 : in out Quad_Double_VecVecs.VecVec;
                  w,wv,wl : out Quad_Double_Vectors.Vector ) is

    n : constant integer32 := logx(r)'length;
    logx1 : Quad_Double_Vectors.Vector(0..r);
    rt1,rt2 : Matrix(1..r-1,1..r-1);
    srp,dsp : Quad_Double_Vectors.Vector(1..r-1) := (1..r-1 => create(0.0));
    p : Quad_Double_Vectors.Vector(0..r-1) := (0..r-1 => create(0.0));
    l,v : Quad_Double_Vectors.Vector(0..r) := (0..r => create(0.0));
    error : Quad_Double_Vectors.Vector(1..r) := (1..r => create(0.0));

  begin
    for i in rt1'range(1) loop
      for j in rt1'range(2) loop
        rt1(i,j) := create(0.0); rt2(i,j) := create(0.0);
      end loop;
    end loop;
    for i in logx(r)'range loop
      for j in logx1'range loop
        logx1(j) := logx(j)(i);
      end loop;
      vLpRs_pipe(r,s(0..r),logs(0..r),logx1(0..r),srp,dsp,p,l,v,rt1,rt2);
      w(i) := v(r)/l(r);
      for j in 1..r loop
        wvl1(j)(i) := v(j)/l(j);
        error(j) := abs(wvl1(j)(i)-wvl0(j)(i));
      end loop;
    end loop;
    wv := v; wl := l;
  end vLpRs_Extrapolate;

  procedure vLpRs_Extrapolate
                ( file : in file_type; r : in integer32;
                  s,logs : in Quad_Double_Vectors.Vector;
                  logx,wvl0 : in Quad_Double_VecVecs.VecVec;
                  wvl1 : in out Quad_Double_VecVecs.VecVec;
                  w,wv,wl : out Quad_Double_Vectors.Vector ) is

    n : constant integer32 := logx(r)'length;
    logx1 : Quad_Double_Vectors.Vector(0..r);
    rt1,rt2 : Matrix(1..r-1,1..r-1);
    srp,dsp : Quad_Double_Vectors.Vector(1..r-1) := (1..r-1 => create(0.0));
    p : Quad_Double_Vectors.Vector(0..r-1) := (0..r-1 => create(0.0));
    l,v : Quad_Double_Vectors.Vector(0..r) := (0..r => create(0.0));
    error : Quad_Double_Vectors.Vector(1..r) := (1..r => create(0.0));

  begin
    for i in rt1'range(1) loop
      for j in rt1'range(2) loop
        rt1(i,j) := create(0.0); rt2(i,j) := create(0.0);
      end loop;
    end loop;
    for i in logx(r)'range loop
      for j in logx1'range loop
        logx1(j) := logx(j)(i);
      end loop;
      vLpRs_pipe(file,r,s(0..r),logs(0..r),logx1(0..r),srp,dsp,p,l,v,rt1,rt2);
      w(i) := v(r)/l(r);
      put(file,s(r),3);
      for j in 1..r loop
        wvl1(j)(i) := v(j)/l(j);
        error(j) := abs(wvl1(j)(i)-wvl0(j)(i));
        put(file," "); put(file,error(j),3);
      end loop;
      new_line(file);
    end loop;
    wv := v; wl := l;
  end vLpRs_Extrapolate;

-- HIGHER-ORDER EXTRAPOLATION (target routines) :

  procedure Affine_Update_Direction
                ( r,m,estm : in out integer32;
                  cntm : in out natural32; thresm : in natural32;
                  er : in out integer32; t,target : in Complex_Number;
                  x : in QuadDobl_Complex_Vectors.Vector;
                  dt,s,logs : in out Quad_Double_Vectors.Vector;
                  logx,wvl0,wvl1,wvltmp : in out Quad_Double_VecVecs.VecVec;
                  v,diferr : in out Quad_Double_Vectors.Vector;
                  error : in out quad_double ) is

    use Quad_Double_Vectors;
    res,errorv : Quad_Double_Vectors.Vector(v'range);
    wv,wl : Quad_Double_Vectors.Vector(0..r);
    ind : integer32 := 1;      -- errors based on first-order extrapolation
    rat,eps : quad_double;
    newm : boolean := false;

  begin
    Affine_Update_Extrapolation_Data
      (r,m,t,target,x,dt,s,logs,logx,wvl0,wvl1,wvltmp);
    if r >= 1 then
      vLpRs_Extrapolate(r,s,logs,logx,wvl1,wvltmp,res,wv,wl);
      if r = 1 then diferr(0) := create(1.0); end if;
      for i in errorv'range loop
        errorv(i) := abs(wvltmp(ind)(i) - wvl1(ind)(i));
      end loop;
      Shift_Up(diferr,Sum_Norm(errorv));
      if er < diferr'last-1 then er := er+1; end if;
    end if;
    if er >= 1 and (diferr(0) < diferr(1)) then
      Estimate_Winding_Number
        (er,thresm,m,estm,cntm,dt(r)/dt(r-1),diferr,rat,eps,newm);
      if newm then
        Refresh_Window(r,m,dt,s,logs);
        vLpRs_Extrapolate(r,s,logs,logx,wvl1,wvltmp,res,wv,wl);
        for i in errorv'range loop
          errorv(i) := abs(wvltmp(ind)(i) - wvl1(ind)(i));
        end loop;
        Shift_Up(diferr,Sum_Norm(errorv)); er := -2;
      end if;
    end if;
    if r >= 1 then
      v := res;
      Update_Errors(r,errorv,error,wvl0,wvl1,wvltmp);
    end if;
    if r < s'last then r := r+1; end if;
  end Affine_Update_Direction;

  procedure Affine_Update_Direction
                ( file : in file_type; 
                  r,m,estm : in out integer32; 
                  cntm : in out natural32; thresm : in natural32;
                  er : in out integer32; t,target : in Complex_Number;
                  x : in QuadDobl_Complex_Vectors.Vector;
                  dt,s,logs : in out Quad_Double_Vectors.Vector;
                  logx,wvl0,wvl1,wvltmp : in out Quad_Double_VecVecs.VecVec;
                  v,diferr : in out Quad_Double_Vectors.Vector;
                  error : in out quad_double ) is

    use Quad_Double_Vectors;
   -- res,errorv,newerrv : Quad_Double_Vectors.Vector(v'range);
    res,errorv : Quad_Double_Vectors.Vector(v'range);
    wv,wl : Quad_Double_Vectors.Vector(0..r);
    ind : integer32 := 1;      -- errors based on first-order extrapolation
    rat,eps : quad_double;
   -- newerr : quad_double;
   -- mv,cntmv,estmv : Standard_Integer_Vectors.Vector(v'range);
    newm : boolean := false;

  begin
    Affine_Update_Extrapolation_Data
      (r,m,t,target,x,dt,s,logs,logx,wvl0,wvl1,wvltmp);
    Write_Update_Information(file,r,m,s(r),logs(r),logx(r).all);
    if r >= 1 then
      vLpRs_Extrapolate(file,r,s,logs,logx,wvl1,wvltmp,res,wv,wl);
      if r = 1 then diferr(0) := create(1.0); end if;
      for i in errorv'range loop
        errorv(i) := abs(wvltmp(ind)(i) - wvl1(ind)(i));
      end loop;
      Shift_Up(diferr,Sum_Norm(errorv));
      if er < diferr'last-1 then er := er+1; end if;
    end if;
    if er >= 1 and (diferr(0) < diferr(1)) then
     -- Estimate0(r,thresm,m,estm,cntm,diferr(1),diferr(0),rat,eps,newm);
      Estimate_Winding_Number
        (file,er,thresm,m,estm,cntm,dt(r)/dt(r-1),diferr,rat,eps,newm);
      put(file,"Ratio for m : "); put(file,rat,3);
      put(file," and accuracy : "); put(file,eps,3); new_line(file);
      if newm then
        Refresh_Window(r,m,dt,s,logs);
        put_line(file,"Extrapolation after adjusting the m-value :");
        vLpRs_Extrapolate(file,r,s,logs,logx,wvl1,wvltmp,res,wv,wl);
        for i in errorv'range loop
          errorv(i) := abs(wvltmp(ind)(i) - wvl1(ind)(i));
        end loop;
        Shift_Up(diferr,Sum_Norm(errorv)); er := -2;
        put(file,"old direction : "); put(file,v); new_line(file);
        put(file,"new direction : "); put(file,res); new_line(file);
      end if;
    end if;
    if r >= 1 then
      v := res;
      Update_Errors(r,errorv,error,wvl0,wvl1,wvltmp);
    end if;
    if r < s'last then r := r+1; end if;
  end Affine_Update_Direction;

  procedure Projective_Update_Direction
                ( r,m,estm : in out integer32;
                  cntm : in out natural32; thresm : in natural32;
                  er : in out integer32; t,target : in Complex_Number;
                  x : in QuadDobl_Complex_Vectors.Vector;
                  dt,s,logs : in out Quad_Double_Vectors.Vector;
                  logx : in out Quad_Double_VecVecs.VecVec;
                  prevv,v : in out Quad_Double_Vectors.Vector;
                  error : in out quad_double ) is

    use Quad_Double_Vectors;
    res : Quad_Double_Vectors.Vector(v'range);
   -- eps : quad_double;
    newerr : quad_double;
    newm : boolean := false;

  begin
    Projective_Update_Extrapolation_Data(r,m,t,target,x,dt,s,logs,logx);
    if r >= 1 then
      res := vLpRs_Extrapolate(r,s,logs,logx);
      error := Sum_Norm(res-prevv); newerr := Sum_Norm(res-v);
      prevv := v;
      v := res;
    end if;
   -- if r >= 2 and (newerr < error)
   --  then Estimate(r,r,thresm,m,estm,cntm,dt,s,logs,error,newerr,eps,newm);
   -- end if;
    if r < s'last
     then r := r+1;
    end if;
    error := newerr;
  end Projective_Update_Direction;

  procedure Projective_Update_Direction
                ( file : in file_type;
                  r,m,estm : in out integer32;
                  cntm : in out natural32; thresm : in natural32;
                  er : in out integer32; t,target : in Complex_Number;
                  x : in QuadDobl_Complex_Vectors.Vector;
                  dt,s,logs : in out Quad_Double_Vectors.Vector;
                  logx : in out Quad_Double_VecVecs.VecVec;
                  prevv,v : in out Quad_Double_Vectors.Vector;
                  error : in out quad_double ) is

    use Quad_Double_Vectors;
    res : Quad_Double_Vectors.Vector(v'range);
   -- eps : quad_double;
    newerr : quad_double;
    newm : boolean := false;

  begin
    Projective_Update_Extrapolation_Data(r,m,t,target,x,dt,s,logs,logx);
    Write_Update_Information(file,r,m,s(r),logs(r),logx(r).all);
    if r >= 1 then
      res := vLpRs_Extrapolate(file,r,s,logs,logx);
      error := Sum_Norm(res-prevv);
      newerr := Sum_Norm(res-v);
      prevv := v;
      v := res;
    end if;
    if r < s'last
     then r := r+1;
    end if;
   -- if r >= 2 and (newerr < error)
   --  then Estimate(r,r,thresm,m,estm,cntm,dt,s,logs,error,newerr,eps,newm);
   -- end if;
    error := newerr;
  end Projective_Update_Direction;

end Directions_of_QuadDobl_Paths;
