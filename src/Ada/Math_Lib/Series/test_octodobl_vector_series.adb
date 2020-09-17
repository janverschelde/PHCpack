with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with OctoDobl_Complex_Series;
with OctoDobl_Complex_Series_io;         use OctoDobl_Complex_Series_io;
with OctoDobl_Complex_Vector_Series;
with OctoDobl_Complex_Vector_Series_io;  use OctoDobl_Complex_Vector_Series_io;
with OctoDobl_CSeries_Vector_Norms;
with OctoDobl_Random_Series_Vectors;

package body Test_OctoDobl_Vector_Series is

  procedure Write ( v : in OctoDobl_Complex_Series_Vectors.Vector ) is

  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      put(v(i));
    end loop;
  end Write;

  procedure OctoDobl_Test_Norm 
              ( v : in OctoDobl_Complex_Series_Vectors.Vector ) is

    use OctoDobl_Complex_Series;
    use OctoDobl_Complex_Series_Vectors;
    use OctoDobl_CSeries_Vector_Norms;

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
  end OctoDobl_Test_Norm;

  procedure OctoDobl_Test ( dim,degree : in integer32 ) is

    sv : constant OctoDobl_Complex_Series_Vectors.Vector(1..dim)
       := OctoDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,degree);
    vs : constant OctoDobl_Complex_Vector_Series.Vector(degree)
       := OctoDobl_Complex_Vector_Series.Create(sv);

  begin
    Write(sv);
    put_line("The coefficients of the vector series :"); put(vs);
    OctoDobl_Test_Norm(sv);
  end OctoDobl_Test;

  procedure Main is

    dim,degree : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(degree);
    OctoDobl_Test(dim,degree);
  end Main;

end Test_OctoDobl_Vector_Series;
