with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with TripDobl_Complex_Series;
with TripDobl_Complex_Series_io;         use TripDobl_Complex_Series_io;
with TripDobl_Complex_Vector_Series;
with TripDobl_Complex_Vector_Series_io;  use TripDobl_Complex_Vector_Series_io;
with TripDobl_CSeries_Vector_Norms;
with TripDobl_Random_Series_Vectors;

package body Test_TripDobl_Vector_Series is

  procedure Write ( v : in TripDobl_Complex_Series_Vectors.Vector ) is

  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      put(v(i));
    end loop;
  end Write;

  procedure TripDobl_Test_Norm 
              ( v : in TripDobl_Complex_Series_Vectors.Vector ) is

    use TripDobl_Complex_Series;
    use TripDobl_Complex_Series_Vectors;
    use TripDobl_CSeries_Vector_Norms;

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
  end TripDobl_Test_Norm;

  procedure TripDobl_Test ( dim,degree : in integer32 ) is

    sv : constant TripDobl_Complex_Series_Vectors.Vector(1..dim)
       := TripDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,degree);
    vs : constant TripDobl_Complex_Vector_Series.Vector(degree)
       := TripDobl_Complex_Vector_Series.Create(sv);

  begin
    Write(sv);
    put_line("The coefficients of the vector series :"); put(vs);
    TripDobl_Test_Norm(sv);
  end TripDobl_Test;

  procedure Main is

    dim,degree : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(degree);
    TripDobl_Test(dim,degree);
  end Main;

end Test_TripDobl_Vector_Series;
