with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with DecaDobl_Complex_Series;
with DecaDobl_Complex_Series_io;         use DecaDobl_Complex_Series_io;
with DecaDobl_Complex_Vector_Series;
with DecaDobl_Complex_Vector_Series_io;  use DecaDobl_Complex_Vector_Series_io;
with DecaDobl_CSeries_Vector_Norms;
with DecaDobl_Random_Series_Vectors;

package body Test_DecaDobl_Vector_Series is

  procedure Write ( v : in DecaDobl_Complex_Series_Vectors.Vector ) is

  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      put(v(i));
    end loop;
  end Write;

  procedure DecaDobl_Test_Norm 
              ( v : in DecaDobl_Complex_Series_Vectors.Vector ) is

    use DecaDobl_Complex_Series;
    use DecaDobl_Complex_Series_Vectors;
    use DecaDobl_CSeries_Vector_Norms;

    sn : constant Series := Norm(v);
    snv : Series(sn.deg);
    nv : Vector(v'range);

  begin
    new_line;
    put_line("The norm of v : "); put(sn);
    nv := Normalize(v);
    put_line("The normalized vector : "); Write(nv);
    snv := Square_of_Norm(nv);
    put_line("Square of norm of normalized vector :"); put(snv);
  end DecaDobl_Test_Norm;

  procedure DecaDobl_Test ( dim,degree : in integer32 ) is

    sv : constant DecaDobl_Complex_Series_Vectors.Vector(1..dim)
       := DecaDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,degree);
    vs : constant DecaDobl_Complex_Vector_Series.Vector(degree)
       := DecaDobl_Complex_Vector_Series.Create(sv);

  begin
    Write(sv);
    put_line("The coefficients of the vector series :"); put(vs);
    DecaDobl_Test_Norm(sv);
  end DecaDobl_Test;

  procedure Main is

    dim,degree : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(degree);
    DecaDobl_Test(dim,degree);
  end Main;

end Test_DecaDobl_Vector_Series;
