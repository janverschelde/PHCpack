with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Functions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Functions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Functions;
with QuadDobl_Speelpenning_Convolutions;
with Series_Polynomial_Gradients;        use Series_Polynomial_Gradients;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;

procedure ts_hesscnv is

-- DESCRIPTION :
--   Tests the computation of the second derivatives of a polynomial
--   stored as a convolution circuit.

  procedure Standard_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate a random circuit for testing in double precision.

    use Standard_Complex_Series;
    use Standard_Complex_Series_Vectors;
    use Standard_CSeries_Polynomials;
    use Standard_Speelpenning_Convolutions;

    c : constant Circuit(nbr,dim,dim-1,dim-2)
      := Standard_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    p : constant Poly := Standard_Polynomial(c);
    q,h : Poly;
    d : Standard_Complex_Numbers.Complex_Number;
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    sx : Standard_Complex_Series_Vectors.Vector(1..dim);
    idx1,idx2 : integer32 := 0;
    eva : Standard_Complex_Series.Link_to_Series;

  begin
    put_line("The exponents : ");
    for i in c.xps'range loop
      put(c.xps(i).all); new_line;
    end loop;
    put("Give the first index : "); get(idx1);
    put("Give the second index : "); get(idx2);
    d := Diff(c,x,idx1,idx2);
    put("The value of the second derivative w.r.t. ");
    put(idx1,1); put(" and "); put(idx2,1); put_line(" :");
    put(d); new_line;
    q := Diff(p,idx1);
    h := Diff(q,idx2);
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    eva := Standard_CSeries_Poly_Functions.Eval(h,sx);
    put_line("The value computed via symbolic differentiation :");
    put(eva.cff(0)); new_line;
  end Standard_Test;

  procedure DoblDobl_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate a random circuit for testing in double double precision.

    use DoblDobl_Complex_Series;
    use DoblDobl_Complex_Series_Vectors;
    use DoblDobl_CSeries_Polynomials;
    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuit(nbr,dim,dim-1,dim-2)
      := DoblDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    p : constant Poly := DoblDobl_Polynomial(c);
    q,h : Poly;
    d : DoblDobl_Complex_Numbers.Complex_Number;
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    sx : DoblDobl_Complex_Series_Vectors.Vector(1..dim);
    idx1,idx2 : integer32 := 0;
    eva : DoblDobl_Complex_Series.Link_to_Series;

  begin
    put_line("The exponents : ");
    for i in c.xps'range loop
      put(c.xps(i).all); new_line;
    end loop;
    put("Give the first index : "); get(idx1);
    put("Give the second index : "); get(idx2);
    d := Diff(c,x,idx1,idx2);
    put("The value of the second derivative w.r.t. ");
    put(idx1,1); put(" and "); put(idx2,1); put_line(" :");
    put(d); new_line;
    q := Diff(p,idx1);
    h := Diff(q,idx2);
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    eva := DoblDobl_CSeries_Poly_Functions.Eval(h,sx);
    put_line("The value computed via symbolic differentiation :");
    put(eva.cff(0)); new_line;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate a random circuit for testing in quad double precision.

    use QuadDobl_Complex_Series;
    use QuadDobl_Complex_Series_Vectors;
    use QuadDobl_CSeries_Polynomials;
    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuit(nbr,dim,dim-1,dim-2)
      := QuadDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    p : constant Poly := QuadDobl_Polynomial(c);
    q,h : Poly;
    d : QuadDobl_Complex_Numbers.Complex_Number;
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    sx : QuadDobl_Complex_Series_Vectors.Vector(1..dim);
    idx1,idx2 : integer32 := 0;
    eva : QuadDobl_Complex_Series.Link_to_Series;

  begin
    put_line("The exponents : ");
    for i in c.xps'range loop
      put(c.xps(i).all); new_line;
    end loop;
    put("Give the first index : "); get(idx1);
    put("Give the second index : "); get(idx2);
    d := Diff(c,x,idx1,idx2);
    put("The value of the second derivative w.r.t. ");
    put(idx1,1); put(" and "); put(idx2,1); put_line(" :");
    put(d); new_line;
    q := Diff(p,idx1);
    h := Diff(q,idx2);
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    eva := QuadDobl_CSeries_Poly_Functions.Eval(h,sx);
    put_line("The value computed via symbolic differentiation :");
    put(eva.cff(0)); new_line;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension, degree, number of terms,
  --   and largest power to generate a random circuit.

    dim,deg,nbr,pwr : integer32 := 0;
    precision : character;

  begin
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    put("Give the number of terms : "); get(nbr);
    put("Give the largest power : "); get(pwr);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    case precision is
      when '0' => Standard_Test(dim,deg,nbr,pwr);
      when '1' => DoblDobl_Test(dim,deg,nbr,pwr);
      when '2' => QuadDobl_Test(dim,deg,nbr,pwr);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hesscnv;
