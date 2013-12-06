package body Multprec_Extrapolators is

  function Extrapolate ( t,t0,t1,x0,x1 : Complex_Number )
                       return Complex_Number is

  --  This is plain linear extrapolation, via evaluation
  --  of the formula x0 + (x1-x0)*(t-t0)/(t1-t0).

    t10 : Complex_Number := t1 - t0;
    x10 : Complex_Number := x1 - x0;
    q10 : Complex_Number := x10/t10;
    dt0 : Complex_Number := t - t0;
    res : Complex_Number := q10*dt0;

  begin
    Add(res,x0);
    Clear(t10); Clear(x10); Clear(q10); Clear(dt0);
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,x0,x1,x2 : Complex_Number )
                       return Complex_Number is

  -- Quadratic extrapolation with divided differences.

    t10 : Complex_Number := t1 - t0;
    t20 : Complex_Number := t2 - t0;
    t21 : Complex_Number := t2 - t1;
    x10 : Complex_Number := x1 - x0;
    x20 : Complex_Number := x2 - x0;
    dt0 : Complex_Number := t - t0;
    dt1 : Complex_Number := t - t1;
    q01 : Complex_Number := x10/t10;
    q02 : Complex_Number := x20/t20;
    dq2 : Complex_Number := q02 - q01;
    q012 : Complex_Number := dq2/t21;
    res : Complex_Number; -- evaluate x0 + dt0*(q01 + q012*dt1)

  begin
    res := q012*dt1;
    Add(res,q01); Mul(res,dt0);
    Add(res,x0);
    Clear(t10); Clear(t20); Clear(t21); Clear(x10); Clear(x20);
    Clear(dt0); Clear(dt1); Clear(q01); Clear(q02); Clear(dq2);
    Clear(q012);
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3,x0,x1,x2,x3 : Complex_Number )
                       return Complex_Number is

  -- Cubic extrapolation with divided differences.

    t10 : Complex_Number := t1 - t0;
    t20 : Complex_Number := t2 - t0;
    t30 : Complex_Number := t3 - t0;
    t21 : Complex_Number := t2 - t1;
    t31 : Complex_Number := t3 - t1;
    t32 : Complex_Number := t3 - t2;
    x10 : Complex_Number := x1 - x0;
    x20 : Complex_Number := x2 - x0;
    x30 : Complex_Number := x3 - x0;
    dt0 : Complex_Number := t - t0;
    dt1 : Complex_Number := t - t1;
    dt2 : Complex_Number := t - t2;
    q01 : Complex_Number := x10/t10;
    q02 : Complex_Number := x20/t20;
    q03 : Complex_Number := x30/t30;
    dq2 : Complex_Number := q02 - q01;
    dq3 : Complex_Number := q03 - q01;
    q012 : Complex_Number := dq2/t21;
    q013 : Complex_Number := dq3/t31;
    dq23 : Complex_Number := q013 - q012;
    q0123 : Complex_Number := dq23/t32;
    res : Complex_Number;  -- evaluate x0 + dt0*(q01 + dt1*(q012 + dt2*q0123))

  begin
    res := dt2*q0123;
    Add(res,q012); Mul(res,dt1);
    Add(res,q01);  Mul(res,dt0);
    Add(res,x0);
    Clear(t10); Clear(t20); Clear(t30); Clear(t21); Clear(t31); Clear(t32);
    Clear(x10); Clear(x20); Clear(x30); Clear(dt0); Clear(dt1); Clear(dt2);
    Clear(q01); Clear(q02); Clear(q03); Clear(dq2); Clear(dq3);
    Clear(q012); Clear(q013); Clear(dq23); Clear(q0123);
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3,t4,x0,x1,x2,x3,x4 : Complex_Number )
                       return Complex_Number is

  -- Quartic extrapolation with divided differences.

    t10 : Complex_Number := t1 - t0;
    t20 : Complex_Number := t2 - t0;
    t30 : Complex_Number := t3 - t0;
    t40 : Complex_Number := t4 - t0;
    t21 : Complex_Number := t2 - t1;
    t31 : Complex_Number := t3 - t1;
    t41 : Complex_Number := t4 - t1;
    t32 : Complex_Number := t3 - t2;
    t42 : Complex_Number := t4 - t2;
    t43 : Complex_Number := t4 - t3;
    x10 : Complex_Number := x1 - x0;
    x20 : Complex_Number := x2 - x0;
    x30 : Complex_Number := x3 - x0;
    x40 : Complex_Number := x4 - x0;
    dt0 : Complex_Number := t - t0;
    dt1 : Complex_Number := t - t1;
    dt2 : Complex_Number := t - t2;
    dt3 : Complex_Number := t - t3;
    q01 : Complex_Number := x10/t10;
    q02 : Complex_Number := x20/t20;
    q03 : Complex_Number := x30/t30;
    q04 : Complex_Number := x40/t40;
    dq2 : Complex_Number := q02 - q01;
    dq3 : Complex_Number := q03 - q01;
    dq4 : Complex_Number := q04 - q01;
    q012 : Complex_Number := dq2/t21;
    q013 : Complex_Number := dq3/t31;
    q014 : Complex_Number := dq4/t41;
    dq23 : Complex_Number := q013 - q012;
    dq24 : Complex_Number := q014 - q012;
    q0123 : Complex_Number := dq23/t32;
    q0124 : Complex_Number := dq24/t42;
    dq34 : Complex_Number := q0124 - q0123;
    q01234 : Complex_Number := dq34/t43;
    res : Complex_Number;
    -- evaluate x0 + dt0*(q01 + dt1*(q012 + dt2*(q0123 + dt3*q01234)))

  begin
    res := dt3*q01234;
    Add(res,q0123); Mul(res,dt2);
    Add(res,q012);  Mul(res,dt1);
    Add(res,q01);   Mul(res,dt0);
    Add(res,x0);
    Clear(t10); Clear(t20); Clear(t30); Clear(t40); 
    Clear(t21); Clear(t31); Clear(t41); Clear(t32);
    Clear(t42); Clear(t43);
    Clear(x10); Clear(x20); Clear(x30); Clear(x40);
    Clear(dt0); Clear(dt1); Clear(dt2); Clear(dt3);
    Clear(q01); Clear(q02); Clear(q03); Clear(q04);
    Clear(dq2); Clear(dq3); Clear(dq4);
    Clear(q012); Clear(q013); Clear(q014);
    Clear(dq23); Clear(dq24);
    Clear(q0123); Clear(q0124); Clear(dq34); Clear(q01234);
    return res;
  end Extrapolate;

  function Extrapolate
             ( t,t0,t1,t2,t3,t4,t5,x0,x1,x2,x3,x4,x5 : Complex_Number )
             return Complex_Number is

  -- Quintic extrapolation with divided differences.

    t10 : Complex_Number := t1 - t0;
    t20 : Complex_Number := t2 - t0;
    t30 : Complex_Number := t3 - t0;
    t40 : Complex_Number := t4 - t0;
    t50 : Complex_Number := t5 - t0;
    t21 : Complex_Number := t2 - t1;
    t31 : Complex_Number := t3 - t1;
    t41 : Complex_Number := t4 - t1;
    t51 : Complex_Number := t5 - t1;
    t32 : Complex_Number := t3 - t2;
    t42 : Complex_Number := t4 - t2;
    t52 : Complex_Number := t5 - t2;
    t43 : Complex_Number := t4 - t3;
    t53 : Complex_Number := t5 - t3;
    t54 : Complex_Number := t5 - t4;
    x10 : Complex_Number := x1 - x0;
    x20 : Complex_Number := x2 - x0;
    x30 : Complex_Number := x3 - x0;
    x40 : Complex_Number := x4 - x0;
    x50 : Complex_Number := x5 - x0;
    dt0 : Complex_Number := t - t0;
    dt1 : Complex_Number := t - t1;
    dt2 : Complex_Number := t - t2;
    dt3 : Complex_Number := t - t3;
    dt4 : Complex_Number := t - t4;
    q01 : Complex_Number := x10/t10;
    q02 : Complex_Number := x20/t20;
    q03 : Complex_Number := x30/t30;
    q04 : Complex_Number := x40/t40;
    q05 : Complex_Number := x50/t50;
    dq2 : Complex_Number := q02 - q01;
    dq3 : Complex_Number := q03 - q01;
    dq4 : Complex_Number := q04 - q01;
    dq5 : Complex_Number := q05 - q01;
    q012 : Complex_Number := dq2/t21;
    q013 : Complex_Number := dq3/t31;
    q014 : Complex_Number := dq4/t41;
    q015 : Complex_Number := dq5/t51;
    dq23 : Complex_Number := q013 - q012;
    dq24 : Complex_Number := q014 - q012;
    dq25 : Complex_Number := q015 - q012;
    q0123 : Complex_Number := dq23/t32;
    q0124 : Complex_Number := dq24/t42;
    q0125 : Complex_Number := dq25/t52;
    dq34 : Complex_Number := q0124 - q0123;
    dq35 : Complex_Number := q0125 - q0123;
    q01234 : Complex_Number := dq34/t43;
    q01235 : Complex_Number := dq35/t53;
    dq45 : Complex_Number := q01235 - q01234;
    q012345 : Complex_Number := dq45/t54;
    res : Complex_Number;
    -- evaluate x0 + dt0*(q01
    --             + dt1*(q012 + dt2*(q0123 + dt3*(q01234 + dt4*q012345))))

  begin
    res := dt4*q012345;
    Add(res,q01234); Mul(res,dt3);
    Add(res,q0123);  Mul(res,dt2);
    Add(res,q012);   Mul(res,dt1);
    Add(res,q01);    Mul(res,dt0);
    Add(res,x0);
    Clear(t10); Clear(t20); Clear(t30); Clear(t40); Clear(t50);
    Clear(t21); Clear(t31); Clear(t41); Clear(t51);
    Clear(t32); Clear(t42); Clear(t52);
    Clear(t43); Clear(t53);
    Clear(t54);
    Clear(x10); Clear(x20); Clear(x30); Clear(x40); Clear(x50);
    Clear(dt0); Clear(dt1); Clear(dt2); Clear(dt3); Clear(dt4);
    Clear(q01); Clear(q02); Clear(q03); Clear(q04); Clear(q05);
    Clear(dq2); Clear(dq3); Clear(dq4); Clear(dq5);
    Clear(q012); Clear(q013); Clear(q014); Clear(q015);
    Clear(dq23); Clear(dq24); Clear(dq25);
    Clear(q0123); Clear(q0124); Clear(q0125);
    Clear(dq34); Clear(dq35);
    Clear(q01234); Clear(q01235);
    Clear(dq45); Clear(q012345);
    return res;
  end Extrapolate;

-- VECTOR VERSIONS :

  function Extrapolate ( t,t0,t1 : Complex_Number; x0,x1 : Vector )
                       return Vector is

  -- This is plain linear extrapolation, via evaluation
  -- of the formula x0 + (x1-x0)*(t-t0)/(t1-t0).

    res : Vector(x0'range);
    t10 : Complex_Number := t1 - t0;
    dt0 : Complex_Number := t - t0;
    x10,q10 : Complex_Number;

  begin
    for i in res'range loop
      x10 := x1(i) - x0(i);
      q10 := x10/t10;
      res(i) := q10*dt0;     -- res(i) := x0(i) + q10*dt0;
      Add(res(i),x0(i));
      Clear(x10); Clear(q10);
    end loop;
    Clear(t10); Clear(dt0);
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2 : Complex_Number; x0,x1,x2 : Vector )
                       return Vector is

  -- Quadratic extrapolation with divided differences.

    res : Vector(x0'range);
    t10 : Complex_Number := t1 - t0;
    t20 : Complex_Number := t2 - t0;
    t21 : Complex_Number := t2 - t1;
    dt0 : Complex_Number := t - t0;
    dt1 : Complex_Number := t - t1;
    x10,x20,q01,q02,dq2,q012 : Complex_Number;

  begin
    for i in x0'range loop
      x10 := x1(i) - x0(i); q01 := x10/t10;
      x20 := x2(i) - x0(i); q02 := x20/t20;
      dq2 := q02 - q01; q012 := dq2/t21;
      res(i) := q012*dt1;   -- res(i) := x0(i) + dt0*(q01 + q012*dt1);
      Add(res(i),q01);
      Mul(res(i),dt0);
      Add(res(i),x0(i));
      Clear(x10); Clear(x20); Clear(q01); Clear(q02);
      Clear(dq2); Clear(q012);
    end loop;
    Clear(t10); Clear(t20); Clear(t21); Clear(dt0); Clear(dt1);
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3 : Complex_Number;
                         x0,x1,x2,x3 : Vector ) return Vector is

  -- Cubic extrapolation with divided differences.

    res : Vector(x0'range);
    t10 : Complex_Number := t1 - t0;
    t20 : Complex_Number := t2 - t0;
    t30 : Complex_Number := t3 - t0;
    t21 : Complex_Number := t2 - t1;
    t31 : Complex_Number := t3 - t1;
    t32 : Complex_Number := t3 - t2;
    dt0 : Complex_Number := t - t0;
    dt1 : Complex_Number := t - t1;
    dt2 : Complex_Number := t - t2;
    x10,x20,x30,q01,q02,q03 : Complex_Number;
    dq2,dq3,q012,q013,dq23,q0123 : Complex_Number;

  begin
    for i in res'range loop
      x10 := x1(i) - x0(i); q01 := x10/t10;
      x20 := x2(i) - x0(i); q02 := x20/t20;
      x30 := x3(i) - x0(i); q03 := x30/t30;
      dq2 := q02 - q01; q012 := dq2/t21;
      dq3 := q03 - q01; q013 := dq3/t31;
      dq23 := q013 - q012; q0123 := dq23/t32;
     -- res(i) := x0(i) + dt0*(q01 + dt1*(q012 + dt2*q0123));
      res(i) := dt2*q0123;
      Add(res(i),q012); Mul(res(i),dt1);
      Add(res(i),q01);  Mul(res(i),dt0);
      Add(res(i),x0(i));
      Clear(x10); Clear(x20); Clear(x30);
      Clear(q01); Clear(q02); Clear(q03); 
      Clear(dq2); Clear(dq3); Clear(dq23);
      Clear(q012); Clear(q013); Clear(q0123); 
    end loop;
    Clear(t10); Clear(t20); Clear(t30);
    Clear(t21); Clear(t31); Clear(t32);
    Clear(dt0); Clear(dt1); Clear(dt2);
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3,t4 : Complex_Number;
                         x0,x1,x2,x3,x4 : Vector ) return Vector is

  -- Quartic extrapolation with divided differences.

    res : Vector(x0'range);
    t10 : Complex_Number := t1 - t0;
    t20 : Complex_Number := t2 - t0;
    t30 : Complex_Number := t3 - t0;
    t40 : Complex_Number := t4 - t0;
    t21 : Complex_Number := t2 - t1;
    t31 : Complex_Number := t3 - t1;
    t41 : Complex_Number := t4 - t1;
    t32 : Complex_Number := t3 - t2;
    t42 : Complex_Number := t4 - t2;
    t43 : Complex_Number := t4 - t3;
    dt0 : Complex_Number := t - t0;
    dt1 : Complex_Number := t - t1;
    dt2 : Complex_Number := t - t2;
    dt3 : Complex_Number := t - t3;
    x10,x20,x30,x40,q01,q02,q03,q04 : Complex_Number;
    dq2,dq3,dq4,q012,q013,q014 : Complex_Number;
    dq23,dq24,q0123,q0124,dq34,q01234 : Complex_Number;

  begin
    for i in res'range loop
      x10 := x1(i) - x0(i); q01 := x10/t10;
      x20 := x2(i) - x0(i); q02 := x20/t20;
      x30 := x3(i) - x0(i); q03 := x30/t30;
      x40 := x4(i) - x0(i); q04 := x40/t40;
      dq2 := q02 - q01; q012 := dq2/t21;
      dq3 := q03 - q01; q013 := dq3/t31;
      dq4 := q04 - q01; q014 := dq4/t41;
      dq23 := q013 - q012; q0123 := dq23/t32;
      dq24 := q014 - q012; q0124 := dq24/t42;
      dq34 := q0124 - q0123; q01234 := dq34/t43;
     -- res(i) := x0(i) + dt0*(q01 + dt1*(q012 + dt2*(q0123 + dt3*q01234)));
      res(i) := dt3*q01234;
      Add(res(i),q0123); Mul(res(i),dt2);
      Add(res(i),q012);  Mul(res(i),dt1);
      Add(res(i),q01);   Mul(res(i),dt0);
      Add(res(i),x0(i));
      Clear(x10); Clear(x20); Clear(x30); Clear(x40);
      Clear(q01); Clear(q02); Clear(q03); Clear(q04);
      Clear(dq2); Clear(dq3); Clear(dq4);
      Clear(q012); Clear(q013); Clear(q014);
      Clear(dq23); Clear(dq24); Clear(dq34); 
      Clear(q0123); Clear(q0124); Clear(q01234);
    end loop;
    Clear(t10); Clear(t20); Clear(t30); Clear(t40);
    Clear(t21); Clear(t31); Clear(t41);
    Clear(t32); Clear(t42);
    Clear(t43);
    Clear(dt0); Clear(dt1); Clear(dt2); Clear(dt3);
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3,t4,t5 : Complex_Number;
                         x0,x1,x2,x3,x4,x5 : Vector ) return Vector is

  -- Quintic extrapolation with divided differences.

    res : Vector(x0'range);
    t10 : Complex_Number := t1 - t0;
    t20 : Complex_Number := t2 - t0;
    t30 : Complex_Number := t3 - t0;
    t40 : Complex_Number := t4 - t0;
    t50 : Complex_Number := t5 - t0;
    t21 : Complex_Number := t2 - t1;
    t31 : Complex_Number := t3 - t1;
    t41 : Complex_Number := t4 - t1;
    t51 : Complex_Number := t5 - t1;
    t32 : Complex_Number := t3 - t2;
    t42 : Complex_Number := t4 - t2;
    t52 : Complex_Number := t5 - t2;
    t43 : Complex_Number := t4 - t3;
    t53 : Complex_Number := t5 - t3;
    t54 : Complex_Number := t5 - t4;
    dt0 : Complex_Number := t - t0;
    dt1 : Complex_Number := t - t1;
    dt2 : Complex_Number := t - t2;
    dt3 : Complex_Number := t - t3;
    dt4 : Complex_Number := t - t4;
    x10,x20,x30,x40,x50,q01,q02,q03,q04,q05 : Complex_Number;
    dq2,dq3,dq4,dq5,q012,q013,q014,q015 : Complex_Number;
    dq23,dq24,dq25,q0123,q0124,q0125,dq34,dq35 : Complex_Number;
    q01234,q01235,dq45,q012345 : Complex_Number;

  begin
    for i in res'range loop
      x10 := x1(i) - x0(i); q01 := x10/t10;
      x20 := x2(i) - x0(i); q02 := x20/t20;
      x30 := x3(i) - x0(i); q03 := x30/t30;
      x40 := x4(i) - x0(i); q04 := x40/t40;
      x50 := x5(i) - x0(i); q05 := x50/t50;
      dq2 := q02 - q01; q012 := dq2/t21;
      dq3 := q03 - q01; q013 := dq3/t31;
      dq4 := q04 - q01; q014 := dq4/t41;
      dq5 := q05 - q01; q015 := dq5/t51;
      dq23 := q013 - q012; q0123 := dq23/t32;
      dq24 := q014 - q012; q0124 := dq24/t42;
      dq25 := q015 - q012; q0125 := dq25/t52;
      dq34 := q0124 - q0123; q01234 := dq34/t43;
      dq35 := q0125 - q0123; q01235 := dq35/t53;
      dq45 := q01235 - q01234; q012345 := dq45/t54;
     -- res(i) := x0(i) + dt0*(q01 + dt1*(q012 
     --                 + dt2*(q0123 + dt3*(q01234 + dt4*q012345))));
      res(i) := dt4*q012345;
      Add(res(i),q01234); Mul(res(i),dt3);
      Add(res(i),q0123);  Mul(res(i),dt2);
      Add(res(i),q012);   Mul(res(i),dt1);
      Add(res(i),q01);    Mul(res(i),dt0);
      Add(res(i),x0(i));
      Clear(x10); Clear(x20); Clear(x30); Clear(x40); Clear(x50);
      Clear(q01); Clear(q02); Clear(q03); Clear(q04); Clear(q05);
      Clear(dq2); Clear(dq3); Clear(dq4); Clear(dq5);
      Clear(q012); Clear(q013); Clear(q014); Clear(q015);
      Clear(dq23); Clear(dq24); Clear(dq25);
      Clear(dq34); Clear(dq35);
      Clear(dq45);
      Clear(q0123); Clear(q0124); Clear(q0125);
      Clear(q01234); Clear(q01235); Clear(q012345);
    end loop;
    Clear(t10); Clear(t20); Clear(t30); Clear(t40); Clear(t50);
    Clear(t21); Clear(t31); Clear(t41); Clear(t51);
    Clear(t32); Clear(t42); Clear(t52);
    Clear(t43); Clear(t53);
    Clear(t54);
    Clear(dt0); Clear(dt1); Clear(dt2); Clear(dt3); Clear(dt4);
    return res;
  end Extrapolate;

end Multprec_Extrapolators;
