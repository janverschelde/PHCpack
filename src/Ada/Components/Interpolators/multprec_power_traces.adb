package body Multprec_Power_Traces is

  function Power_Sums_to_Trace
             ( s,t : Vector; k : integer32 ) return Complex_Number is

    sum,acc : Complex_Number;

  begin
    Copy(s(k),sum);
    for i in 1..k-1 loop
      acc := s(i)*t(k-i);
      Add(sum,acc);              -- sum := sum + s(i)*t(k-i)
      Clear(acc);
    end loop;
    acc := Create(k);
    Div(sum,acc);
    Clear(acc);
    Min(sum);
    return sum;
  end Power_Sums_to_Trace;

  function Traces_to_Power_Sum 
             ( t,s : Vector; k : integer32 ) return Complex_Number is

    sum,acc : Complex_Number;

  begin
    Copy(t(k),sum);
    acc := Create(k);
    Mul(sum,acc);
    Clear(acc);
    for i in 1..k-1 loop
      acc := s(i)*t(k-i);
      Add(sum,acc);              -- sum := sum + s(i)*t(k-i)
      Clear(acc);
    end loop;
    Min(sum);
    return sum;
  end Traces_to_Power_Sum;

  function Power_Sums_to_Traces ( s : Vector ) return Vector is

    t : Vector(s'range);

  begin
    t(t'first) := Create(integer(0));        -- avoid compiler warning
    for i in s'range loop
      t(i) := Power_Sums_to_Trace(s,t,i);
    end loop;
    return t;
  end Power_Sums_to_Traces;

  function Traces_to_Power_Sums ( t : Vector ) return Vector is

    s : Vector(t'range);

  begin
    s(s'first) := Create(integer(0));       -- avoid compiler warning
    for i in t'range loop
      s(i) := Traces_to_Power_Sum(t,s,i);
    end loop;
    return s;
  end Traces_to_Power_Sums;

end Multprec_Power_Traces;
