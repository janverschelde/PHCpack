with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;

procedure ts_pslaur is

-- DESCRIPTION :
--   Performs some basic tests on truncated Laurent power series,
--   defined by a leading exponent (which may be negative)
--   and a complex coefficient vector.

  procedure Write ( e : in integer32;
                    c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the power series with leading exponent e
  --   and coefficients in c.

  begin
    for i in c'range loop
      if i > c'first
       then put(" + (");
       else put("   (");
      end if;
      put(c(i)); put(")*t^"); put(e+i,1); new_line;
    end loop;
  end Write;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the degree of truncation
  --   and the leading exponent.

     d,e : integer32 := 0;

  begin
    new_line;
    put("Give the truncation degree : "); get(d);
    put("Give the leading exponent  : "); get(e);
    declare
      c : constant Standard_Complex_Vectors.Vector(0..d)
        := Standard_Random_Vectors.Random_Vector(0,d);
    begin
      new_line;
      put_line("A random series :"); Write(e,c);
    end;
  end Main;

begin
  Main;
end ts_pslaur;
