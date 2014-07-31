with Double_Double_Numbers;             use Double_Double_Numbers;

package body DoblDobl_Power_Traces is

  function Power_Sums_to_Trace
             ( s,t : Vector; k : integer32 ) return Complex_Number is

    sum : Complex_Number := s(k);
    kdd : constant double_double := create(k);

  begin
    for i in 1..k-1 loop
      sum := sum + s(i)*t(k-i);
    end loop;
    sum := sum/Create(kdd);
    return -sum;
  end Power_Sums_to_Trace;

  function Traces_to_Power_Sum 
             ( t,s : Vector; k : integer32 ) return Complex_Number is

    sum : Complex_Number := t(k);
    kdd : constant double_double := create(k);

  begin
    sum := sum*Create(kdd);
    for i in 1..k-1 loop
      sum := sum + s(i)*t(k-i); 
    end loop;
    return -sum;
  end Traces_to_Power_Sum;

  function Power_Sums_to_Traces ( s : Vector ) return Vector is

    t : Vector(s'range);
    zero : constant double_double := create(0.0);

  begin
    t(t'first) := Create(zero);     -- avoid compiler warning
    for i in s'range loop
      t(i) := Power_Sums_to_Trace(s,t,i);
    end loop;
    return t;
  end Power_Sums_to_Traces;

  function Traces_to_Power_Sums ( t : Vector ) return Vector is

    s : Vector(t'range);
    zero : constant double_double := create(0.0);

  begin
    s(s'first) := Create(zero);     -- avoid compiler warning
    for i in t'range loop
      s(i) := Traces_to_Power_Sum(t,s,i);
    end loop;
    return s;
  end Traces_to_Power_Sums;

end DoblDobl_Power_Traces;
