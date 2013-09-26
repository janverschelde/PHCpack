with Standard_Integer_Numbers;           use Standard_integer_Numbers;

package body Multprec_vLpRs_Tables is

-- I. The v-table :

  procedure v_pipe ( v : in out Vector; p : in Vector;
                     vr : in Floating_Number ) is

    tmp0,tmp1,tmp2 : Floating_Number;

  begin
    Copy(v(0),tmp0);
    Copy(vr,v(0));
    for i in 1..v'last loop
      Copy(v(i),tmp1);
      Copy(tmp0,v(i));
      tmp2 := p(i-1)*v(i-1);
      Sub(v(i),tmp2);
      Clear(tmp2);
      Copy(tmp1,tmp0);
    end loop;
    Clear(tmp0);
    Clear(tmp1);
  end v_pipe;

  procedure v_pipe ( v,p : in Vector; vr : in Floating_Number;
                     vrp : in out Vector ) is

    tmp : Floating_Number;

  begin
    Copy(vr,vrp(0));
    for i in 1..v'last loop
      tmp := p(i-1)*vrp(i-1);
      Clear(vrp(i));
      vrp(i) := v(i-1) - tmp;
      Clear(tmp);
    end loop;
  end v_pipe;

-- II. The L-table :

  procedure L_pipe ( l : in out Vector; p : in Vector; 
                     lr : in Floating_Number ) is

    tmp0,tmp1,tmp2 : Floating_Number;

  begin
    Copy(l(0),tmp0);
    Copy(lr,l(0));
    for i in 1..l'last loop
      Copy(l(i),tmp1);
      tmp2 := p(i-1)*l(i-1);
      Clear(l(i));
      l(i) := tmp0 - tmp2;
      Clear(tmp2);
      Copy(tmp1,tmp0);
    end loop;
    Clear(tmp0);
    Clear(tmp1);
  end L_pipe;

  procedure L_pipe ( l,p : in Vector; lr : in Floating_Number;
                     lrp : in out Vector ) is

    tmp : Floating_Number;

  begin
    Copy(lr,lrp(0));
    for i in 1..l'last loop
      tmp := p(i-1)*lrp(i-1);
      Clear(lrp(i));
      lrp(i) := l(i-1) - tmp;
      Clear(tmp);
    end loop;
  end L_pipe;

-- The combined full version for v- and L-table :

  procedure vL_full ( s,l,v : in Vector; srp,dsp,p,lrp,vrp : out Vector;
                      rt1,rt2 : in out Matrix ) is

    srm1,tmpsrp,tmpdsp : Vector(1..s'last-1);
    tmpp : Vector(0..s'last-1);
    prev_lrp,new_lrp,prev_vrp,new_vrp : Vector(s'range);
    k_last : integer32;
    acc : Floating_Number;

  begin
    Copy(s(0),srm1(srm1'first));
    for i in srm1'first+1..srm1'last loop
      Clear(srm1(i));
      srm1(i) := s(0)*srm1(i-1);
    end loop;
    s_pipe(srm1,s(1),tmpsrp,tmpdsp);
    Copy(tmpsrp,srm1);
    for j in rt1'range(2) loop
      Copy(tmpdsp(j),rt1(1,j));
    end loop;
    Clear(tmpp(0));
    tmpp(0) := Create(integer(1));
    Copy(l(0),prev_lrp(0));  Copy(l(1),new_lrp(0));
    Copy(v(0),prev_vrp(0));  Copy(v(1),new_vrp(0));
    Clear(new_vrp(1)); new_vrp(1) := prev_vrp(0) - new_vrp(0);
    Clear(new_lrp(1)); new_lrp(1) := prev_lrp(0) - new_lrp(0);
    Copy(new_lrp(0),prev_lrp(0)); Copy(new_lrp(1),prev_lrp(1));
    Copy(new_vrp(0),prev_vrp(0)); Copy(new_vrp(1),prev_vrp(1));
    for i in 2..s'last loop
      s_pipe(srm1,s(i),tmpsrp,tmpdsp);
      Copy(tmpsrp,srm1);
      for j in rt2'range(2) loop
        Copy(tmpdsp(j),rt2(1,j));
      end loop;
      if i < s'last
       then k_last := i;    -- compute one additional row
       else k_last := i-1;
      end if;
      for k in 1..k_last loop
        Clear(tmpp(k));
        tmpp(k) := rt1(k,k)/rt2(k,k);
        for j in k+1..rt2'last(2) loop
          acc := tmpp(k)*rt2(k,j);
          Clear(rt2(k+1,j));
          rt2(k+1,j) := rt1(k,j) - acc;
          Clear(acc);
        end loop;
      end loop;
      Copy(rt2,rt1);
      Copy(v(i),new_vrp(0));
      Copy(l(i),new_lrp(0));
      for k in 1..i loop
        acc := tmpp(k-1)*new_vrp(k-1);
        Clear(new_vrp(k));
        new_vrp(k) := prev_vrp(k-1) - acc;
        Clear(acc);
        acc := tmpp(k-1)*new_lrp(k-1);
        Clear(new_lrp(k));
        new_lrp(k) := prev_lrp(k-1) - acc;
        Clear(acc);
      end loop;
      for j in 0..i loop
        Copy(new_vrp(j),prev_vrp(j));
        Copy(new_lrp(j),prev_lrp(j));
      end loop;
    end loop;
    Clear(srp); srp := tmpsrp;
    Clear(dsp); dsp := tmpdsp;
    Clear(p);   p := tmpp;
    Clear(lrp); lrp := new_lrp;
    Clear(vrp); vrp := new_vrp;
    Clear(prev_lrp); Clear(prev_vrp);
    Clear(srm1);
  end vL_full;

-- III. The p-table :

  procedure p_full ( s : in Vector; srp,dsp,p : out Vector;
                     rt1,rt2 : in out Matrix ) is
  begin
    RR_full(s,srp,dsp,p,rt1,rt2);
  end p_full;

  procedure p_pipe ( rt1,rt2 : in Matrix; p : out Vector ) is
  begin
    Clear(p(0));
    p(0) := Create(integer(1));
    for i in p'first+1..p'last loop
      Clear(p(i));
      p(i) := rt1(i,i)/rt2(i,i);
    end loop;
  end p_pipe;

-- IV. The R-table :

  procedure R_full ( s : in Vector; srp,dsp,p : out Vector;
                     rt1,rt2 : in out Matrix ) is

    srm1,tmpsrp,tmpdsp : Vector(1..s'last-1);
    tmpp : Vector(0..s'last-1);
    acc : Floating_Number;

  begin
    Copy(s(0),srm1(srm1'first));
    for i in srm1'first+1..srm1'last loop
      Clear(srm1(i));
      srm1(i) := s(0)*srm1(i-1);
    end loop;
    s_pipe(srm1,s(1),tmpsrp,tmpdsp);
    Copy(tmpsrp,srm1);
    for j in tmpdsp'range loop
      Copy(tmpdsp(j),rt1(1,j));
    end loop;
    Clear(tmpp(0));
    tmpp(0) := Create(integer(1));
    for i in 2..s'last loop 
      s_pipe(srm1,s(i),tmpsrp,tmpdsp);
      Copy(tmpsrp,srm1);
      for j in tmpdsp'range loop
        Copy(tmpdsp(j),rt2(1,j));
      end loop;
      for k in 1..i-1 loop
        Clear(tmpp(k));
        tmpp(k) := rt1(k,k)/rt2(k,k);
        for j in k+1..rt2'last(2) loop
          acc := tmpp(k)*rt2(k,j);
          Clear(rt2(k+1,j));
          rt2(k+1,j) := rt1(k,j) - acc;
          Clear(acc);
        end loop;
      end loop;
      Copy(rt2,rt1);
    end loop;
    Clear(srp); srp := tmpsrp;
    Clear(dsp); dsp := tmpdsp;
    Clear(p);   p := tmpp;
    Clear(srm1);
  end R_full;

  procedure RR_full ( s : in Vector; srp,dsp,p : out Vector;
                      rt1,rt2 : in out Matrix ) is

    srm1,tmpsrp,tmpdsp : Vector(1..s'last-1);
    tmpp : Vector(0..s'last-1);
    k_last,j_first : integer32;
    acc : Floating_Number;

  begin
    Copy(s(0),srm1(srm1'first));
    for i in srm1'first+1..srm1'last loop
      Clear(srm1(i));
      srm1(i) := s(0)*srm1(i-1);
    end loop;
    s_pipe(srm1,s(1),tmpsrp,tmpdsp);
    Copy(tmpsrp,srm1);
    for j in tmpdsp'range loop
      Clear(rt1(1,j));
      Copy(tmpdsp(j),rt1(1,j));
    end loop;
    Clear(tmpp(0));
    tmpp(0) := Create(integer(1));
    for i in 2..s'last loop
      s_pipe(srm1,s(i),tmpsrp,tmpdsp);
      Copy(tmpsrp,srm1);
      for j in tmpdsp'range loop
        Clear(rt2(1,j));
        Copy(tmpdsp(j),rt2(1,j));
      end loop;
      if i < s'last
       then k_last := i;    -- compute one additional row
       else k_last := i-1;
      end if;
      for k in 1..k_last loop
        Clear(tmpp(k));
        tmpp(k) := rt1(k,k)/rt2(k,k);
        if k < rt2'last(2)
         then j_first := k;
         else j_first := k+1;
        end if;
        for j in j_first..rt2'last(2) loop   -- start earlier with columns
          acc := tmpp(k)*rt2(k,j);
          Clear(rt2(k+1,j));
          rt2(k+1,j) := rt1(k,j) - acc;
          Clear(acc);
        end loop;
      end loop;
      Copy(rt2,rt1);
    end loop;
    Clear(srp); srp := tmpsrp;
    Clear(dsp); dsp := tmpdsp;
    Clear(p);   p := tmpp;
    Clear(srm1);
  end RR_full;

  procedure R_pipe ( rt1 : in Matrix; s,p : in Vector;
                     rt2 : in out Matrix ) is

    acc : Floating_Number;

  begin
    Copy(s(1),rt2(1,1));
    for j in 2..s'last loop
      Clear(rt2(1,j));
      Copy(s(j),rt2(1,j));
      for i in 2..j loop
        acc := p(i-1)*rt2(i-1,j);
        Clear(rt2(i,j));
        rt2(i,j) := rt1(i-1,j) - acc;
        Clear(acc);
      end loop;
    end loop;
  end R_pipe;

  procedure RR_pipe ( rt1 : in Matrix; s,p : in Vector;
                      rt2 : in out Matrix ) is

    i_last : integer32;
    acc : Floating_Number;

  begin
    Copy(s(1),rt2(1,1));
    for j in 2..s'last loop
      Copy(s(j),rt2(1,j));
      if j < rt2'last(1)
       then i_last := j+1;   -- compute one additional row
       else i_last := j;
      end if;
      for i in 2..i_last loop
        acc := p(i-1)*rt2(i-1,j);
        Clear(rt2(i,j));
        rt2(i,j) := rt1(i-1,j) - acc;
        Clear(acc);
      end loop;
    end loop;
  end RR_pipe;

-- V. The s-table :

  procedure s_full ( s : in Vector; srp,dsp : out Vector ) is

    srm1,tmpsrp : Vector(1..s'last-1);

  begin
    Copy(s(0),srm1(srm1'first));
    for i in srm1'first+1..srm1'last loop
      Clear(srm1(i));
      srm1(i) := s(0)*srm1(i-1);
    end loop;
    for i in 1..s'last loop
      s_pipe(srm1,s(i),tmpsrp,dsp);
      Copy(tmpsrp,srm1);
    end loop;
    Copy(tmpsrp,srp);
    Clear(srm1);
    Clear(tmpsrp);
  end s_full;

  procedure s_pipe ( srp : in out Vector; sr : in Floating_Number;
                     dsp : out Vector ) is

    srpow,tmp : Floating_Number;

  begin
    Copy(srp(1),tmp);
    Copy(sr,srpow);
    Copy(srpow,srp(1));
    Clear(dsp(1));
    dsp(1) := srpow - tmp;
    for i in 2..srp'last loop
      Mul(srpow,sr);
      Copy(srp(i),tmp);
      Copy(srpow,srp(i));
      Clear(dsp(i));
      dsp(i) := srpow - tmp;
    end loop;
    Clear(tmp);
    Clear(srpow);
  end s_pipe;

  procedure s_pipe ( sr1 : in Vector; sr : in Floating_Number;
                     srp,dsp : out Vector ) is

    srpow : Floating_Number;

  begin
    Copy(sr,srpow);
    Copy(srpow,srp(1));
    Clear(dsp(1));
    dsp(1) := srpow - sr1(1);
    for i in sr1'first+1..sr1'last loop
      Mul(srpow,sr);
      Copy(srpow,srp(i));
      Clear(dsp(i));
      dsp(i) := srpow - sr1(i);
    end loop;
    Clear(srpow);
  end s_pipe;

end Multprec_vLpRs_Tables;
