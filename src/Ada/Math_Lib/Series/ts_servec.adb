with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Dense_Series;               use Standard_Dense_Series;
with Standard_Dense_Series_io;            use Standard_Dense_Series_io;
with Standard_Dense_Series_Vectors;       use Standard_Dense_Series_Vectors;
with Standard_Random_Series;
with Standard_Series_Vector_Norms;        use Standard_Series_Vector_Norms;

procedure ts_servec is

-- DESCRIPTION :
--   Test on vectors of truncated power series.

  procedure Write ( v : in Vector ) is

  -- DESCRIPTION :
  --   Writes the components of the vector to standard output.
 
  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      put(v(i));
    end loop;
  end Write;

  procedure Test_Norm ( v : in Vector ) is

  -- DESCRIPTION :
  --   Test on normalization of the vector v.

    sn : constant Series := Norm(v);
    snv : Series;
    nv : Vector(v'range);

  begin
    new_line;
    put_line("The norm of v : "); put(sn);
    nv := Normalize(v);
    put_line("The normalized vector : "); Write(nv);
    snv := Square_of_Norm(nv);
    put_line("Square of norm of normalized vector :"); put(snv);
  end Test_Norm;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension, an order
  --   and then generates a random vector.

    dim,order : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the order : "); get(order);
    declare
      v : Vector(1..dim)
        := Standard_Random_Series.Random_Series_Vector(1,dim,order);
    begin
      Write(v);
      Test_Norm(v);
    end;
  end Main;

begin
  Main;
end ts_servec;
