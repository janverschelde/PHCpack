package body Standard_Power_Traces is

  function Power_Sums_to_Trace
             ( s,t : Vector; k : integer32 ) return Complex_Number is

    sum : Complex_Number := s(k);

  begin
    for i in 1..k-1 loop
      sum := sum + s(i)*t(k-i);
    end loop;
    sum := sum/Create(k);
    return -sum;
  end Power_Sums_to_Trace;

  function Traces_to_Power_Sum 
             ( t,s : Vector; k : integer32 ) return Complex_Number is

    sum : Complex_Number := t(k);

  begin
    sum := sum*Create(k);
    for i in 1..k-1 loop
      sum := sum + s(i)*t(k-i); 
    end loop;
    return -sum;
  end Traces_to_Power_Sum;

  function Power_Sums_to_Traces ( s : Vector ) return Vector is

    t : Vector(s'range);

  begin
    t(t'first) := Create(0.0);     -- avoid compiler warning
    for i in s'range loop
      t(i) := Power_Sums_to_Trace(s,t,i);
    end loop;
    return t;
  end Power_Sums_to_Traces;

  function Traces_to_Power_Sums ( t : Vector ) return Vector is

    s : Vector(t'range);

  begin
    s(s'first) := Create(0.0);     -- avoid compiler warning
    for i in t'range loop
      s(i) := Traces_to_Power_Sum(t,s,i);
    end loop;
    return s;
  end Traces_to_Power_Sums;

end Standard_Power_Traces;
