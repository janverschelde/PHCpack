with Double_Double_Basics;               use Double_Double_Basics;

package body Quad_Double_Renormalizations is

  procedure quick_renorm ( c0,c1,c2,c3,c4 : in out double_float ) is

    s,t0,t1,t2,t3 : double_float;

  begin
    quick_two_sum(c3,c4,s,t3);
    quick_two_sum(c2,s,s,t2);
    quick_two_sum(c1,s,s,t1);
    quick_two_sum(c0,s,c0,t0);
    quick_two_sum(t2,t3,s,t2);
    quick_two_sum(t1,s,s,t1);
    quick_two_sum(t0,s,c1,t0);
    quick_two_sum(t1,t2,s,t1);
    quick_two_sum(t0,s,c2,t0);
    c3 := t0 + t1;
  end quick_renorm;

  procedure renorm4 ( c0,c1,c2,c3 : in out double_float ) is

    s0,s1 : double_float;
    s2 : double_float := 0.0;
    s3 : double_float := 0.0;

  begin
    quick_two_sum(c2,c3,s0,c3);
    quick_two_sum(c1,s0,s0,c2);
    quick_two_sum(c0,s0,c0,c1);
    s0 := c0; s1 := c1;
    if s1 /= 0.0 then
      quick_two_sum(s1,c2,s1,s2);
      if s2 /= 0.0
       then quick_two_sum(s2,c3,s2,s3);
       else quick_two_sum(s1,c3,s1,s2);
      end if;
    else
      quick_two_sum(s0,c2,s0,s1);
      if s1 /= 0.0
       then quick_two_sum(s1,c3,s1,s2);
       else quick_two_sum(s0,c3,s0,s1);
      end if;
    end if;
    c0 := s0; c1 := s1; c2 := s2; c3 := s3;
  end renorm4;

  procedure renorm5 ( c0,c1,c2,c3,c4 : in out double_float ) is

    s0,s1 : double_float;
    s2 : double_float := 0.0;
    s3 : double_float := 0.0;

  begin
    quick_two_sum(c3,c4,s0,c4);
    quick_two_sum(c2,s0,s0,c3);
    quick_two_sum(c1,s0,s0,c2);
    quick_two_sum(c0,s0,c0,c1);
    s0 := c0; s1 := c1;
    quick_two_sum(c0,c1,s0,s1);
    if s1 /= 0.0 then
      quick_two_sum(s1,c2,s1,s2);
      if s2 /= 0.0 then
        quick_two_sum(s2,c3,s2,s3);
        if s3 /= 0.0
         then s3 := s3 + c4;
         else s2 := s2 + c4;
        end if;
      else
        quick_two_sum(s1,c3,s1,s2);
        if s2 /= 0.0
         then quick_two_sum(s2,c4,s2,s3);
         else quick_two_sum(s1,c4,s1,s2);
        end if;
      end if;
    else
      quick_two_sum(s0,c2,s0,s1);
      if s1 /= 0.0 then
        quick_two_sum(s1,c3,s1,s2);
        if s2 /= 0.0
         then quick_two_sum(s2,c4,s2,s3);
         else quick_two_sum(s1,c4,s1,s2);
        end if;
      else
        quick_two_sum(s0,c3,s0,s1);
        if s1 /= 0.0 
         then quick_two_sum(s1,c4,s1,s2);
         else quick_two_sum(s0,c4,s0,s1);
        end if;
      end if;
    end if;
    c0 := s0; c1 := s1; c2 := s2; c3 := s3;
  end renorm5;

  procedure three_sum ( a,b,c : in out double_float ) is

    t1,t2,t3 : double_float;

  begin
    two_sum(a,b,t1,t2);
    two_sum(c,t1,a,t3);
    two_sum(t2,t3,b,c);
  end three_sum;

  procedure three_sum2 ( a,b,c : in out double_float ) is

    t1,t2,t3 : double_float;

  begin
    two_sum(a,b,t1,t2);
    two_sum(c,t1,a,t3);
    b := t2 + t3;
  end three_sum2;

  procedure quick_three_accum
              ( a,b,s : in out double_float; c : in double_float ) is

    za,zb : boolean;

  begin
    two_sum(b,c,s,b);
    two_sum(a,s,s,a);
    za := (a /= 0.0);
    zb := (b /= 0.0);
    if za and zb then
      return;
    else
      if not zb
       then b := a; a := s;
       else a := s;
      end if;
      s := 0.0;
    end if;
  end quick_three_accum;

end Quad_Double_Renormalizations;
