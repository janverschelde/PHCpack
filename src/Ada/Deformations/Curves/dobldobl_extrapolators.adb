package body DoblDobl_Extrapolators is

  function Extrapolate ( t,t0,t1,x0,x1 : Complex_Number )
                       return Complex_Number is

  --  This is plain linear extrapolation, via evaluation
  --  of the formula x0 + (x1-x0)*(t-t0)/(t1-t0).

    t10 : constant Complex_Number := t1 - t0;
    x10 : constant Complex_Number := x1 - x0;
    q10 : constant Complex_Number := x10/t10;
    dt0 : constant Complex_Number := t - t0;
    res : constant Complex_Number := x0 + q10*dt0;

  begin
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,x0,x1,x2 : Complex_Number )
                       return Complex_Number is

  -- Quadratic extrapolation with divided differences.

    t10 : constant Complex_Number := t1 - t0;
    t20 : constant Complex_Number := t2 - t0;
    t21 : constant Complex_Number := t2 - t1;
    x10 : constant Complex_Number := x1 - x0;
    x20 : constant Complex_Number := x2 - x0;
    dt0 : constant Complex_Number := t - t0;
    dt1 : constant Complex_Number := t - t1;
    q01 : constant Complex_Number := x10/t10;
    q02 : constant Complex_Number := x20/t20;
    dq2 : constant Complex_Number := q02 - q01;
    q012 : constant Complex_Number := dq2/t21;
    res : constant Complex_Number := x0 + dt0*(q01 + q012*dt1);

  begin
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3,x0,x1,x2,x3 : Complex_Number )
                       return Complex_Number is

  -- Cubic extrapolation with divided differences.

    t10 : constant Complex_Number := t1 - t0;
    t20 : constant Complex_Number := t2 - t0;
    t30 : constant Complex_Number := t3 - t0;
    t21 : constant Complex_Number := t2 - t1;
    t31 : constant Complex_Number := t3 - t1;
    t32 : constant Complex_Number := t3 - t2;
    x10 : constant Complex_Number := x1 - x0;
    x20 : constant Complex_Number := x2 - x0;
    x30 : constant Complex_Number := x3 - x0;
    dt0 : constant Complex_Number := t - t0;
    dt1 : constant Complex_Number := t - t1;
    dt2 : constant Complex_Number := t - t2;
    q01 : constant Complex_Number := x10/t10;
    q02 : constant Complex_Number := x20/t20;
    q03 : constant Complex_Number := x30/t30;
    dq2 : constant Complex_Number := q02 - q01;
    dq3 : constant Complex_Number := q03 - q01;
    q012 : constant Complex_Number := dq2/t21;
    q013 : constant Complex_Number := dq3/t31;
    dq23 : constant Complex_Number := q013 - q012;
    q0123 : constant Complex_Number := dq23/t32;
    res : constant Complex_Number := x0 + dt0*(q01 + dt1*(q012 + dt2*q0123));

  begin
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3,t4,x0,x1,x2,x3,x4 : Complex_Number )
                       return Complex_Number is

  -- Quartic extrapolation with divided differences.

    t10 : constant Complex_Number := t1 - t0;
    t20 : constant Complex_Number := t2 - t0;
    t30 : constant Complex_Number := t3 - t0;
    t40 : constant Complex_Number := t4 - t0;
    t21 : constant Complex_Number := t2 - t1;
    t31 : constant Complex_Number := t3 - t1;
    t41 : constant Complex_Number := t4 - t1;
    t32 : constant Complex_Number := t3 - t2;
    t42 : constant Complex_Number := t4 - t2;
    t43 : constant Complex_Number := t4 - t3;
    x10 : constant Complex_Number := x1 - x0;
    x20 : constant Complex_Number := x2 - x0;
    x30 : constant Complex_Number := x3 - x0;
    x40 : constant Complex_Number := x4 - x0;
    dt0 : constant Complex_Number := t - t0;
    dt1 : constant Complex_Number := t - t1;
    dt2 : constant Complex_Number := t - t2;
    dt3 : constant Complex_Number := t - t3;
    q01 : constant Complex_Number := x10/t10;
    q02 : constant Complex_Number := x20/t20;
    q03 : constant Complex_Number := x30/t30;
    q04 : constant Complex_Number := x40/t40;
    dq2 : constant Complex_Number := q02 - q01;
    dq3 : constant Complex_Number := q03 - q01;
    dq4 : constant Complex_Number := q04 - q01;
    q012 : constant Complex_Number := dq2/t21;
    q013 : constant Complex_Number := dq3/t31;
    q014 : constant Complex_Number := dq4/t41;
    dq23 : constant Complex_Number := q013 - q012;
    dq24 : constant Complex_Number := q014 - q012;
    q0123 : constant Complex_Number := dq23/t32;
    q0124 : constant Complex_Number := dq24/t42;
    dq34 : constant Complex_Number := q0124 - q0123;
    q01234 : constant Complex_Number := dq34/t43;
    res : constant Complex_Number
        := x0 + dt0*(q01 + dt1*(q012 + dt2*(q0123 + dt3*q01234)));

  begin
    return res;
  end Extrapolate;

  function Extrapolate
             ( t,t0,t1,t2,t3,t4,t5,x0,x1,x2,x3,x4,x5 : Complex_Number )
             return Complex_Number is

  -- Quintic extrapolation with divided differences.

    t10 : constant Complex_Number := t1 - t0;
    t20 : constant Complex_Number := t2 - t0;
    t30 : constant Complex_Number := t3 - t0;
    t40 : constant Complex_Number := t4 - t0;
    t50 : constant Complex_Number := t5 - t0;
    t21 : constant Complex_Number := t2 - t1;
    t31 : constant Complex_Number := t3 - t1;
    t41 : constant Complex_Number := t4 - t1;
    t51 : constant Complex_Number := t5 - t1;
    t32 : constant Complex_Number := t3 - t2;
    t42 : constant Complex_Number := t4 - t2;
    t52 : constant Complex_Number := t5 - t2;
    t43 : constant Complex_Number := t4 - t3;
    t53 : constant Complex_Number := t5 - t3;
    t54 : constant Complex_Number := t5 - t4;
    x10 : constant Complex_Number := x1 - x0;
    x20 : constant Complex_Number := x2 - x0;
    x30 : constant Complex_Number := x3 - x0;
    x40 : constant Complex_Number := x4 - x0;
    x50 : constant Complex_Number := x5 - x0;
    dt0 : constant Complex_Number := t - t0;
    dt1 : constant Complex_Number := t - t1;
    dt2 : constant Complex_Number := t - t2;
    dt3 : constant Complex_Number := t - t3;
    dt4 : constant Complex_Number := t - t4;
    q01 : constant Complex_Number := x10/t10;
    q02 : constant Complex_Number := x20/t20;
    q03 : constant Complex_Number := x30/t30;
    q04 : constant Complex_Number := x40/t40;
    q05 : constant Complex_Number := x50/t50;
    dq2 : constant Complex_Number := q02 - q01;
    dq3 : constant Complex_Number := q03 - q01;
    dq4 : constant Complex_Number := q04 - q01;
    dq5 : constant Complex_Number := q05 - q01;
    q012 : constant Complex_Number := dq2/t21;
    q013 : constant Complex_Number := dq3/t31;
    q014 : constant Complex_Number := dq4/t41;
    q015 : constant Complex_Number := dq5/t51;
    dq23 : constant Complex_Number := q013 - q012;
    dq24 : constant Complex_Number := q014 - q012;
    dq25 : constant Complex_Number := q015 - q012;
    q0123 : constant Complex_Number := dq23/t32;
    q0124 : constant Complex_Number := dq24/t42;
    q0125 : constant Complex_Number := dq25/t52;
    dq34 : constant Complex_Number := q0124 - q0123;
    dq35 : constant Complex_Number := q0125 - q0123;
    q01234 : constant Complex_Number := dq34/t43;
    q01235 : constant Complex_Number := dq35/t53;
    dq45 : constant Complex_Number := q01235 - q01234;
    q012345 : constant Complex_Number := dq45/t54;
    res : constant Complex_Number
        := x0 + dt0*(q01
              + dt1*(q012 + dt2*(q0123 + dt3*(q01234 + dt4*q012345))));

  begin
    return res;
  end Extrapolate;

-- VECTOR VERSIONS :

  function Extrapolate ( t,t0,t1 : Complex_Number; x0,x1 : Vector )
                       return Vector is

  -- This is plain linear extrapolation, via evaluation
  -- of the formula x0 + (x1-x0)*(t-t0)/(t1-t0).

    res : Vector(x0'range);
    t10 : constant Complex_Number := t1 - t0;
    dt0 : constant Complex_Number := t - t0;
    x10,q10 : Complex_Number;

  begin
    for i in res'range loop
      x10 := x1(i) - x0(i); q10 := x10/t10;
      res(i) := x0(i) + q10*dt0;
    end loop;
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2 : Complex_Number; x0,x1,x2 : Vector )
                       return Vector is

  -- Quadratic extrapolation with divided differences.

    res : Vector(x0'range);
    t10 : constant Complex_Number := t1 - t0;
    t20 : constant Complex_Number := t2 - t0;
    t21 : constant Complex_Number := t2 - t1;
    dt0 : constant Complex_Number := t - t0;
    dt1 : constant Complex_Number := t - t1;
    x10,x20,q01,q02,dq2,q012 : Complex_Number;

  begin
    for i in x0'range loop
      x10 := x1(i) - x0(i); q01 := x10/t10;
      x20 := x2(i) - x0(i); q02 := x20/t20;
      dq2 := q02 - q01; q012 := dq2/t21;
      res(i) := x0(i) + dt0*(q01 + q012*dt1);
    end loop;
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3 : Complex_Number;
                         x0,x1,x2,x3 : Vector ) return Vector is

  -- Cubic extrapolation with divided differences.

    res : Vector(x0'range);
    t10 : constant Complex_Number := t1 - t0;
    t20 : constant Complex_Number := t2 - t0;
    t30 : constant Complex_Number := t3 - t0;
    t21 : constant Complex_Number := t2 - t1;
    t31 : constant Complex_Number := t3 - t1;
    t32 : constant Complex_Number := t3 - t2;
    dt0 : constant Complex_Number := t - t0;
    dt1 : constant Complex_Number := t - t1;
    dt2 : constant Complex_Number := t - t2;
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
      res(i) := x0(i) + dt0*(q01 + dt1*(q012 + dt2*q0123));
    end loop;
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3,t4 : Complex_Number;
                         x0,x1,x2,x3,x4 : Vector ) return Vector is

  -- Quartic extrapolation with divided differences.

    res : Vector(x0'range);
    t10 : constant Complex_Number := t1 - t0;
    t20 : constant Complex_Number := t2 - t0;
    t30 : constant Complex_Number := t3 - t0;
    t40 : constant Complex_Number := t4 - t0;
    t21 : constant Complex_Number := t2 - t1;
    t31 : constant Complex_Number := t3 - t1;
    t41 : constant Complex_Number := t4 - t1;
    t32 : constant Complex_Number := t3 - t2;
    t42 : constant Complex_Number := t4 - t2;
    t43 : constant Complex_Number := t4 - t3;
    dt0 : constant Complex_Number := t - t0;
    dt1 : constant Complex_Number := t - t1;
    dt2 : constant Complex_Number := t - t2;
    dt3 : constant Complex_Number := t - t3;
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
      res(i) := x0(i) + dt0*(q01 + dt1*(q012 + dt2*(q0123 + dt3*q01234)));
    end loop;
    return res;
  end Extrapolate;

  function Extrapolate ( t,t0,t1,t2,t3,t4,t5 : Complex_Number;
                         x0,x1,x2,x3,x4,x5 : Vector ) return Vector is

  -- Quintic extrapolation with divided differences.

    res : Vector(x0'range);
    t10 : constant Complex_Number := t1 - t0;
    t20 : constant Complex_Number := t2 - t0;
    t30 : constant Complex_Number := t3 - t0;
    t40 : constant Complex_Number := t4 - t0;
    t50 : constant Complex_Number := t5 - t0;
    t21 : constant Complex_Number := t2 - t1;
    t31 : constant Complex_Number := t3 - t1;
    t41 : constant Complex_Number := t4 - t1;
    t51 : constant Complex_Number := t5 - t1;
    t32 : constant Complex_Number := t3 - t2;
    t42 : constant Complex_Number := t4 - t2;
    t52 : constant Complex_Number := t5 - t2;
    t43 : constant Complex_Number := t4 - t3;
    t53 : constant Complex_Number := t5 - t3;
    t54 : constant Complex_Number := t5 - t4;
    dt0 : constant Complex_Number := t - t0;
    dt1 : constant Complex_Number := t - t1;
    dt2 : constant Complex_Number := t - t2;
    dt3 : constant Complex_Number := t - t3;
    dt4 : constant Complex_Number := t - t4;
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
      res(i) := x0(i) + dt0*(q01 + dt1*(q012 
                      + dt2*(q0123 + dt3*(q01234 + dt4*q012345))));
    end loop;
    return res;
  end Extrapolate;



end DoblDobl_Extrapolators;
