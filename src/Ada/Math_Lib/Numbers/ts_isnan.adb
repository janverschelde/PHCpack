with text_io;                           use text_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;

procedure ts_isnan is

-- DESCRIPTION :
--   Check if NaN can be detected.

  num,den,quo : double_float := 0.0;
  one : constant double_float := 1.0;
  tol : constant double_float := 1.0E-8;

  function Is_Valid ( x : double_float ) return boolean is

  -- DESCRIPTION :
  --   A valid number is either positive, negative, or zero,
  --   in which case true is returned, otherwise false is returned.

  begin
    if x > 0.0 then
      return true;
    elsif x < 0.0 then
      return true;
    elsif x = 0.0 then
      return true;
    else
      return false;
    end if;
  end Is_Valid;

begin
  put("give numerator (suggestion: 0) : "); get(num);
  put("give denominator (suggestion: 0) : "); get(den);
  declare
  begin
    quo := num/den;
  exception
    when others => put_line("ignoring the exception ...");
  end;
  put("quotient = "); put(quo); new_line;
  if one = quo
   then put_line("1.0 equals the quotient");
   else put_line("1.0 does not equal the quotient");
  end if;
  put("1 - quotient : "); put(one - quo); new_line;
  if one - quo < tol
   then put_line("1 - quotient < 1.0E-8");
   else put_line("1 - quotient >= 1.0E-8");
  end if;
  put("abs(1 - quotient) : "); put(abs(one - quo)); new_line;
  if abs(one - quo) < tol
   then put_line("abs(1 - quotient) < 1.0E-8");
   else put_line("abs(1 - quotient) >= 1.0E-8");
  end if;
  if abs(one - quo) = 0.0 then
    put_line("abs(1 - quotient) = 0.0");
  elsif abs(one - quo) > 0.0 then
    put_line("abs(1 - quotient) > 0.0");
  else
    put_line("abs(1 - quotient) < 0.0");
  end if;
  if abs(quo - one) = 0.0 then
    put_line("abs(quotient - 1) = 0.0");
  elsif abs(quo - one) > 0.0 then
    put_line("abs(quotient - 1) > 0.0");
  else
    put_line("abs(quotient - 1) < 0.0");
  end if;
  if quo < 0.0
   then put_line("quo < 0");
   else put_line("quo >= 0");
  end if;
  if quo > 0.0
   then put_line("quo > 0");
   else put_line("quo <= 0");
  end if;
  if quo = 0.0
   then put_line("quo = 0");
   else put_line("quo /= 0");
  end if;
  if Standard_Floating_Numbers.is_valid(double_float(1.0))
   then put_line("1 is valid");
   else put_line("1 is not valid");
  end if;
  if Standard_Floating_Numbers.is_valid(quo)
   then put_line("quo is valid");
   else put_line("quo is not valid");
  end if;
end ts_isnan;
