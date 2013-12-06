with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;   
-- with Multprec_Floating_Matrices_io;      use Multprec_Floating_Matrices_io;
with Multprec_vLpRs_Tables;              use Multprec_vLpRs_Tables;

package body Multprec_vLpRs_Algorithm is

-- OUTPUT ROUTINES OF ERROR :

  procedure Write_Init ( file : in file_type; s,L,v : in Vector ) is

  -- DESCRIPTION :
  --   Writes the beginning of the error table.

    tmp,err,mpw,acc : Floating_Number;
    frw : double_float;
    w : integer32;

  begin
    tmp := v(1)/L(1);
    frw := Round(tmp);
    w := integer32(frw);
    mpw := Create(w);
    acc := tmp - mpw;
    err := AbsVal(acc);
    put(file,s(0),2,3,3); new_line(file);
    put(file,s(1),2,3,3); put(file,err,2,3,3); new_line(file);
    Clear(tmp); Clear(err); Clear(mpw); Clear(acc);
  end Write_Init;

  procedure Write ( file : in file_type; k : in integer32; 
                    s : in Floating_Number; L,v : in Vector ) is

  -- DESCRIPTION :
  --   Writes an additional row of k columns of the error table.
  --   Note that when m is wrong, the outcome is no longer integer.

    tmp,err,mpw,acc : Floating_Number;
    frw : double_float;
    w : integer32;

  begin
    put(file,s,2,3,3);
    tmp := v(k)/L(k);
    frw := Round(tmp); Clear(tmp);
    w := integer32(frw);
    mpw := Create(w);
    for i in 1..k loop
      tmp := v(i)/L(i);
      acc := tmp - mpw;
      err := AbsVal(acc);
      put(file,err,2,3,3);
      Clear(tmp); Clear(acc); Clear(err);
    end loop;
    new_line(file);
    Clear(mpw);
  end Write;

-- TARGET ROUTINES :

  procedure vLpRs_full
                ( r : in integer32; s,logs,logx : in Vector;
                  srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix ) is
  begin
    vL_full(s(0..r),logs(0..r),logx(0..r),srp,dsp,p,L,v,rt1,rt2);
    Copy(rt2,rt1);
    for k in r+1..s'last loop
      vlprs_pipe(s(k),logs(k),logx(k),srp,dsp,p,L,v,rt1,rt2);
    end loop;
  end vLpRs_full;

  procedure vLpRs_pipe
              ( file : in file_type; r : in integer32;
                s,logs,logx : in Vector; srp,dsp,p,L,v : in out Vector;
                rt1,rt2 : in out Matrix ) is
  begin
    Clear(p(0));
    p(0) := Create(integer(1));                         -- initialization
    Copy(logx(0),v(0)); Copy(logx(1),v(1));
    Copy(logs(0),L(0)); Copy(logs(1),L(1));
    L_pipe(L(0..1),p(0..0),logs(1));
    v_pipe(v(0..1),p(0..0),logx(1));
   -- Write_Init(file,s,L,v);                           -- write the error table
    for k in 2..r loop
      p_full(s(0..k),srp(1..k-1),dsp(1..k-1),p(0..k-1),rt1,rt2);
      L_pipe(L(0..k),p(0..k-1),logs(k));                        -- extrapolate
      v_pipe(v(0..k),p(0..k-1),logx(k));
     -- Write(file,k,s(k),L,v);                         -- write the error table
     -- put_line(file,"rt2 :"); put(file,rt2,3,3,3);
    end loop;
    Copy(rt2,rt1);
    for k in r+1..s'last loop
      vlprs_pipe(file,s(k),logs(k),logx(k),srp,dsp,p,L,v,rt1,rt2);
     -- put_line(file,"rt2 : "); put(file,rt2,3,3,3);
    end loop;
  end vLpRs_pipe;

  procedure vLpRs_pipe
              ( s,logs,logx : in Floating_Number;
                srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix ) is
  begin
    s_pipe(srp,s,dsp);
    RR_pipe(rt1,dsp,p,rt2);
    p_pipe(rt1,rt2,p);
    Copy(rt2,rt1);
    L_pipe(L,p,logs);
    v_pipe(v,p,logx);
  end vLpRs_pipe;

  procedure vLpRs_pipe
              ( file : in file_type; s,logs,logx : in Floating_Number;
                srp,dsp,p,L,v : in out Vector; rt1,rt2 : in out Matrix ) is
  begin
    vLpRs_pipe(s,logs,logx,srp,dsp,p,L,v,rt1,rt2);
   -- Write(file,l'last,s,L,v);
  end vLpRs_pipe;

end Multprec_vLpRs_Algorithm;
