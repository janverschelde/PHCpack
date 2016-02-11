with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;   
--with Double_Double_Matrices_io;          use Double_Double_Matrices_io;
with DoblDobl_vLpRs_Tables;              use DoblDobl_vLpRs_Tables;

package body DoblDobl_vLpRs_Algorithm is

-- OUTPUT ROUTINES OF ERROR :

  procedure Write_Init ( file : in file_type; s,L,v : in Vector ) is

  -- DESCRIPTION :
  --   Writes the beginning of the error table.

    tmp,err : double_double;
    w : integer32;
    d_w : double_float;

  begin
    tmp := v(1)/L(1);
    w := integer32(hi_part(tmp));
    d_w := create(w);
    err := abs(tmp - create(d_w));
    put(file,s(0),3); new_line(file);
    put(file,s(1),3); put(file," "); put(file,err,3); new_line(file);
  end Write_Init;

  procedure Write ( file : in file_type; k : in integer32; 
                    s : in double_double; L,v : in Vector ) is

  -- DESCRIPTION :
  --   Writes an additional row of k columns of the error table.
  --   Note that when m is wrong, the outcome is no longer integer.

    tmp,err : double_double;
    w : integer32;
    d_w : double_float;

  begin
    put(file,s,3);
   -- w := v(k)/L(k); -- 
    w := integer32(hi_part(v(k)/L(k)));
    for i in 1..k loop
      tmp := v(i)/L(i);
     -- err := abs(tmp-w); -- 
      d_w := create(w);
      err := abs(tmp - create(d_w));
      put(file," "); put(file,err,3);
    end loop;
    new_line(file);
  end Write;

-- TARGET ROUTINES :

  procedure vLpRs_full
              ( r : in integer32; s,logs,logx : in Vector;
                srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix ) is
  begin
    vL_full(s(0..r),logs(0..r),logx(0..r),srp,dsp,p,L,v,rt1,rt2);
    rt1 := rt2;
    for k in r+1..s'last loop
      vlprs_pipe(s(k),logs(k),logx(k),srp,dsp,p,L,v,rt1,rt2);
    end loop;
  end vLpRs_full;

  procedure vLpRs_pipe
              ( r : in integer32; s,logs,logx : in Vector;
                srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix ) is
  begin
    p(0) := create(1.0);                                   -- initialization
    v(0..1) := logx(0..1);
    L(0..1) := logs(0..1);
    L_pipe(l(0..1),p(0..0),logs(1));
    v_pipe(v(0..1),p(0..0),logx(1));
    for k in 2..r loop
      p_full(s(0..k),srp(1..k-1),dsp(1..k-1),p(0..k-1),rt1,rt2);
      L_pipe(L(0..k),p(0..k-1),logs(k));                        -- extrapolate
      v_pipe(v(0..k),p(0..k-1),logx(k));
    end loop;
    rt1 := rt2;
    for k in r+1..s'last loop
      vlprs_pipe(s(k),logs(k),logx(k),srp,dsp,p,L,v,rt1,rt2);
    end loop;
  end vLpRs_pipe;

  procedure vLpRs_pipe
              ( file : in file_type;
                r : in integer32; s,logs,logx : in Vector;
                srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix ) is
  begin
    p(0) := create(1.0);                                   -- initialization
    v(0..1) := logx(0..1);
    L(0..1) := logs(0..1);
    L_pipe(l(0..1),p(0..0),logs(1));
    v_pipe(v(0..1),p(0..0),logx(1));
   -- Write_Init(file,s,l,v);                           -- write the error table
    for k in 2..r loop
      p_full(s(0..k),srp(1..k-1),dsp(1..k-1),p(0..k-1),rt1,rt2);
      L_pipe(L(0..k),p(0..k-1),logs(k));                        -- extrapolate
      v_pipe(v(0..k),p(0..k-1),logx(k));
     -- Write(file,k,s(k),l,v);                         -- write the error table
     -- put_line(file,"rt2 :"); put(file,rt2,3,3,3);
    end loop;
    rt1 := rt2;
    for k in r+1..s'last loop
      vlprs_pipe(file,s(k),logs(k),logx(k),srp,dsp,p,L,v,rt1,rt2);
     -- put_line(file,"rt2 : "); put(file,rt2,3,3,3);
    end loop;
  end vLpRs_pipe;

  procedure vLpRs_pipe
              ( s,logs,logx : in double_double;
                srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix ) is
  begin
    s_pipe(srp,s,dsp);
    RR_pipe(rt1,dsp,p,rt2);
    p_pipe(rt1,rt2,p);  rt1 := rt2;
    L_pipe(L,p,logs);
    v_pipe(v,p,logx);
  end vLpRs_pipe;

  procedure vLpRs_pipe
              ( file : in file_type; s,logs,logx : in double_double;
                srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix ) is
  begin
    vLpRs_pipe(s,logs,logx,srp,dsp,p,L,v,rt1,rt2);
   -- Write(file,L'last,s,L,v);
  end vLpRs_pipe;

end DoblDobl_vLpRs_Algorithm;
