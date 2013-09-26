with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body DoblDobl_vLpRs_Tables is

-- I. The v-table :

  procedure v_pipe ( v : in out Vector; p : in Vector;
                     vr : in double_double ) is

    tmp0,tmp1 : double_double;

  begin
    tmp0 := v(0);
    v(0) := vr;
    for i in 1..v'last loop
      tmp1 := v(i);
      v(i) := tmp0 - p(i-1)*v(i-1);
      tmp0 := tmp1;
    end loop;
  end v_pipe;

  procedure v_pipe ( v,p : in Vector; vr : in double_double;
                     vrp : in out Vector ) is
  begin
    vrp(0) := vr;
    for i in 1..v'last loop
      vrp(i) := v(i-1) - p(i-1)*vrp(i-1);
    end loop;
  end v_pipe;

-- II. The L-table :

  procedure L_pipe ( l : in out Vector; p : in Vector; 
                     lr : in double_double ) is

    tmp0,tmp1 : double_double;

  begin
    tmp0 := l(0);
    l(0) := lr;
    for i in 1..l'last loop
      tmp1 := l(i);
      l(i) := tmp0 - p(i-1)*l(i-1);
      tmp0 := tmp1;
    end loop;
  end L_pipe;

  procedure L_pipe ( l,p : in Vector; lr : in double_double;
                     lrp : in out Vector ) is
  begin
    lrp(0) := lr;
    for i in 1..l'last loop
      lrp(i) := l(i-1) - p(i-1)*lrp(i-1);
    end loop;
  end L_pipe;

-- The combined full version for v- and L-table :

  procedure vL_full ( s,l,v : in Vector; srp,dsp,p,lrp,vrp : out Vector;
                      rt1,rt2 : in out Matrix ) is

    srm1,tmpsrp,tmpdsp : Vector(1..s'last-1);
    tmpp : Vector(0..s'last-1);
    prev_lrp,new_lrp,prev_vrp,new_vrp : Vector(s'range);
    k_last : integer32;

  begin
    srm1(srm1'first) := s(0);
    for i in srm1'first+1..srm1'last loop
      srm1(i) := s(0)*srm1(i-1);
    end loop;
    s_pipe(srm1,s(1),tmpsrp,tmpdsp); srm1 := tmpsrp;
    for j in rt1'range(2) loop
      rt1(1,j) := tmpdsp(j);
    end loop;
    tmpp(0) := create(1.0);
    prev_lrp(0) := l(0);  new_lrp(0) := l(1);
    prev_vrp(0) := v(0);  new_vrp(0) := v(1);
    new_vrp(1) := prev_vrp(0) - new_vrp(0);
    new_lrp(1) := prev_lrp(0) - new_lrp(0);
    prev_lrp(0..1) := new_lrp(0..1);
    prev_vrp(0..1) := new_vrp(0..1);
    for i in 2..s'last loop
      s_pipe(srm1,s(i),tmpsrp,tmpdsp); srm1 := tmpsrp;
      for j in rt2'range(2) loop
        rt2(1,j) := tmpdsp(j);
      end loop;
      if i < s'last
       then k_last := i;    -- compute one additional row
       else k_last := i-1;
      end if;
      for k in 1..k_last loop
        tmpp(k) := rt1(k,k)/rt2(k,k);
        for j in k+1..rt2'last(2) loop
          rt2(k+1,j) := rt1(k,j) - tmpp(k)*rt2(k,j);
        end loop;
      end loop;
      rt1 := rt2;
      new_vrp(0) := v(i);
      new_lrp(0) := l(i);
      for k in 1..i loop
        new_vrp(k) := prev_vrp(k-1) - tmpp(k-1)*new_vrp(k-1);
        new_lrp(k) := prev_lrp(k-1) - tmpp(k-1)*new_lrp(k-1);
      end loop;
      prev_vrp(0..i) := new_vrp(0..i);
      prev_lrp(0..i) := new_lrp(0..i);
    end loop;
    srp := tmpsrp;
    dsp := tmpdsp;
    p := tmpp;
    lrp := new_lrp;
    vrp := new_vrp;
  end vL_full;

-- III. The p-table :

  procedure p_full ( s : in Vector; srp,dsp,p : out Vector;
                     rt1,rt2 : in out Matrix ) is
  begin
    RR_full(s,srp,dsp,p,rt1,rt2);
  end p_full;

  procedure p_pipe ( rt1,rt2 : in Matrix; p : out Vector ) is
  begin
    p(0) := create(1.0);
    for i in p'first+1..p'last loop
      p(i) := rt1(i,i)/rt2(i,i);
    end loop;
  end p_pipe;

-- IV. The R-table :

  procedure R_full ( s : in Vector; srp,dsp,p : out Vector;
                     rt1,rt2 : in out Matrix ) is

    srm1,tmpsrp,tmpdsp : Vector(1..s'last-1);
    tmpp : Vector(0..s'last-1);

  begin
    srm1(srm1'first) := s(0);
    for i in srm1'first+1..srm1'last loop
      srm1(i) := s(0)*srm1(i-1);
    end loop;
    s_pipe(srm1,s(1),tmpsrp,tmpdsp); srm1 := tmpsrp;
    for j in tmpdsp'range loop
      rt1(1,j) := tmpdsp(j);
    end loop;
    tmpp(0) := create(1.0);
    for i in 2..s'last loop 
      s_pipe(srm1,s(i),tmpsrp,tmpdsp); srm1 := tmpsrp;
      for j in tmpdsp'range loop
        rt2(1,j) := tmpdsp(j);
      end loop;
      for k in 1..i-1 loop
        tmpp(k) := rt1(k,k)/rt2(k,k);
        for j in k+1..rt2'last(2) loop
          rt2(k+1,j) := rt1(k,j) - tmpp(k)*rt2(k,j);
        end loop;
      end loop;
      rt1 := rt2;
    end loop;
    srp := tmpsrp;
    dsp := tmpdsp;
    p := tmpp;
  end R_full;

  procedure RR_full ( s : in Vector; srp,dsp,p : out Vector;
                      rt1,rt2 : in out Matrix ) is

    srm1,tmpsrp,tmpdsp : Vector(1..s'last-1);
    tmpp : Vector(0..s'last-1);
    k_last,j_first : integer32;

  begin
    srm1(srm1'first) := s(0);
    for i in srm1'first+1..srm1'last loop
      srm1(i) := s(0)*srm1(i-1);
    end loop;
    s_pipe(srm1,s(1),tmpsrp,tmpdsp); srm1 := tmpsrp;
    for j in tmpdsp'range loop
      rt1(1,j) := tmpdsp(j);
    end loop;
    tmpp(0) := create(1.0);
    for i in 2..s'last loop
      s_pipe(srm1,s(i),tmpsrp,tmpdsp); srm1 := tmpsrp;
      for j in tmpdsp'range loop
        rt2(1,j) := tmpdsp(j);
      end loop;
      if i < s'last
       then k_last := i;    -- compute one additional row
       else k_last := i-1;
      end if;
      for k in 1..k_last loop
        tmpp(k) := rt1(k,k)/rt2(k,k);
        if k < rt2'last(2)
         then j_first := k;
         else j_first := k+1;
        end if;
        for j in j_first..rt2'last(2) loop   -- start earlier with columns
          rt2(k+1,j) := rt1(k,j) - tmpp(k)*rt2(k,j);
        end loop;
      end loop;
      rt1 := rt2;
    end loop;
    srp := tmpsrp;
    dsp := tmpdsp;
    p := tmpp;
  end RR_full;

  procedure R_pipe ( rt1 : in Matrix; s,p : in Vector;
                     rt2 : in out Matrix ) is
  begin
    rt2(1,1) := s(1);
    for j in 2..s'last loop
      rt2(1,j) := s(j);
      for i in 2..j loop
        rt2(i,j) := rt1(i-1,j) - p(i-1)*rt2(i-1,j);
      end loop;
    end loop;
  end R_pipe;

  procedure RR_pipe ( rt1 : in Matrix; s,p : in Vector;
                      rt2 : in out Matrix ) is

    i_last : integer32;

  begin
    rt2(1,1) := s(1);
    for j in 2..s'last loop
      rt2(1,j) := s(j);
      if j < rt2'last(1)
       then i_last := j+1;   -- compute one additional row
       else i_last := j;
      end if;
      for i in 2..i_last loop
        rt2(i,j) := rt1(i-1,j) - p(i-1)*rt2(i-1,j);
      end loop;
    end loop;
  end RR_pipe;

-- V. The s-table :

  procedure s_full ( s : in Vector; srp,dsp : out Vector ) is

    srm1,tmpsrp : Vector(1..s'last-1);

  begin
    srm1(srm1'first) := s(0);
    for i in srm1'first+1..srm1'last loop
      srm1(i) := s(0)*srm1(i-1);
    end loop;
    for i in 1..s'last loop
      s_pipe(srm1,s(i),tmpsrp,dsp);
      srm1 := tmpsrp;
    end loop;
    srp := tmpsrp;
  end s_full;

  procedure s_pipe ( srp : in out Vector; sr : in double_double;
                     dsp : out Vector ) is

    tmp : double_double := srp(1);
    srpow : double_double := sr;

  begin
    srp(1) := srpow;
    dsp(1) := srpow - tmp;
    for i in 2..srp'last loop
      srpow := srpow*sr;
      tmp := srp(i);
      srp(i) := srpow;
      dsp(i) := srpow - tmp;
    end loop;
  end s_pipe;

  procedure s_pipe ( sr1 : in Vector; sr : in double_double;
                     srp,dsp : out Vector ) is

    srpow : double_double := sr;

  begin
    srp(1) := srpow;
    dsp(1) := srpow - sr1(1);
    for i in sr1'first+1..sr1'last loop
      srpow := srpow*sr;
      srp(i) := srpow;
      dsp(i) := srpow - sr1(i);
    end loop;
  end s_pipe;

end DoblDobl_vLpRs_Tables;
