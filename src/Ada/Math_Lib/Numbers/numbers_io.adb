with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;

package body Numbers_io is

  procedure Read_Positive ( p : out positive ) is

    package positive_io is new text_io.integer_io(positive);

  begin
    positive_io.get(p); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a positive number, please try again : ");
      Read_Positive(p);
  end Read_Positive;

  procedure Read_Natural ( n : out natural32 ) is
  begin
    n := 0;
    get(n); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a natural number, please try again : ");
      Read_Natural(n);
  end Read_Natural;

  procedure Read_Integer ( i : out integer32 ) is
  begin
    i := 0;
    get(i); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not an integer number, please try again : ");
      Read_Integer(i);
  end Read_Integer;

  procedure Read_Single_Float ( f : out single_float ) is
  begin
    f := 0.0;
    get(f); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a floating point number, please try again : ");
      Read_Single_Float(f);
  end Read_Single_Float;

  procedure Read_Double_Float ( f : out double_float ) is
  begin
    f := 0.0;
    get(f); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a floating point number, please try again : ");
      Read_Double_Float(f);
  end Read_Double_Float;

  procedure Read_Double_Double ( f : out double_double ) is
  begin
    f := create(0.0);
    get(f); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a double double, please try again : ");
      Read_Double_Double(f);
  end Read_Double_Double;

  procedure Read_Quad_Double ( f : out quad_double ) is
  begin
    f := create(0.0);
    get(f); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a quad double, please try again : ");
      Read_Quad_Double(f);
  end Read_Quad_Double;

  function Number_of_Decimal_Places ( n : natural32 ) return natural32 is
  begin
    if n < 10 then
      return 1;
    elsif n < 100 then
      return 2;
    elsif n < 1000 then
      return 3;
    elsif n < 10000 then
      return 4;
    elsif n < 100000 then
      return 5;
    elsif n < 1000000 then
      return 6;
    elsif n < 10000000 then
      return 7;
    elsif n < 100000000 then
      return 8;
    else
      return 9;
    end if;
  end Number_of_Decimal_Places;

  procedure Read_Positive_Float ( f : in out double_float ) is
  begin
    loop
      Read_Double_Float(f);
      exit when f > 0.0;
      put("Zero or negative value not allowed, please try again : ");
    end loop;
  end Read_Positive_Float;

  procedure Read_Positive_Double_Double ( f : in out double_double ) is
  begin
    loop
      Read_Double_Double(f);
      exit when f > 0.0;
      put("Zero or negative value not allowed, please try again : ");
    end loop;
  end Read_Positive_Double_Double;

  procedure Read_Positive_Quad_Double ( f : in out quad_double ) is
  begin
    loop
      Read_Quad_Double(f);
      exit when f > 0.0;
      put("Zero or negative value not allowed, please try again : ");
    end loop;
  end Read_Positive_Quad_Double;

end Numbers_io;
