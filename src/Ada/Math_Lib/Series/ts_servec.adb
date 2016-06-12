with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Dense_Series_io;            use Standard_Dense_Series_io;
with Standard_Dense_Series_Vectors;
with Standard_Random_Series;

procedure ts_servec is

-- DESCRIPTION :
--   Test on dense series vectors.

  procedure Main is

    v : Standard_Dense_Series_Vectors.Vector(0..2)
      := Standard_Random_Series.Random_Series_Vector(0,2,4);

  begin
    put_line("A vector of random series of order 4 :");
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      put(v(i));
    end loop;
  end Main;

begin
  Main;
end ts_servec;
