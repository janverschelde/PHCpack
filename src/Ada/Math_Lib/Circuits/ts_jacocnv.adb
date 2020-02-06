with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
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
with Jacobian_Convolution_Circuits;

procedure ts_jacocnv is

-- DESCRIPTION :
--   Tests the computation of the first derivatives of polynomials
--   stored as a convolution circuit.

  procedure Standard_Test_System ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate random circuits for testing in double precision.
  --   Computes the Jacobian for the circuit
  --   and compares with the evaluation at the symbolic derivative.

    use Standard_Complex_Numbers;
    use Standard_Complex_Matrices;
    use Standard_Complex_Series;
    use Standard_Complex_Series_Vectors;
    use Standard_CSeries_Polynomials;
    use Standard_Speelpenning_Convolutions;

    c : constant Circuits(1..dim)
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    jm,pj : Matrix(1..dim,1..dim);
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    p,dp : Poly;
    sx : Standard_Complex_Series_Vectors.Vector(1..dim);
    eva : Standard_Complex_Series.Link_to_Series;
    err : Complex_Number;
    val,errsum : double_float := 0.0;

  begin
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    jm := Jacobian_Convolution_Circuits.Jacobian(c,x);
    for i in 1..dim loop
      p := Standard_Polynomial(c(i).all);
      for j in 1..dim loop
        dp := Diff(p,j);
        eva := Standard_CSeries_Poly_Functions.Eval(dp,sx);
        pj(i,j) := eva.cff(0);
        Standard_CSeries_Polynomials.Clear(dp);
        err := pj(i,j) - jm(i,j); val := AbsVal(err);
        errsum := errsum + val;
      end loop;
      Standard_CSeries_Polynomials.Clear(p);
    end loop;
    put("Sum of errors :"); put(errsum,3); new_line;
  end Standard_Test_System;

  procedure DoblDobl_Test_System ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate random circuits for testing in double double precision.
  --   Computes the Jacobian for the circuit
  --   and compares with the evaluation at the symbolic derivative.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Series;
    use DoblDobl_Complex_Series_Vectors;
    use DoblDobl_CSeries_Polynomials;
    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits(1..dim)
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    jm,pj : Matrix(1..dim,1..dim);
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    p,dp : Poly;
    sx : DoblDobl_Complex_Series_Vectors.Vector(1..dim);
    eva : DoblDobl_Complex_Series.Link_to_Series;
    err : Complex_Number;
    val,errsum : double_double := create(0.0);

  begin
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    jm := Jacobian_Convolution_Circuits.Jacobian(c,x);
    for i in 1..dim loop
      p := DoblDobl_Polynomial(c(i).all);
      for j in 1..dim loop
        dp := Diff(p,j);
        eva := DoblDobl_CSeries_Poly_Functions.Eval(dp,sx);
        pj(i,j) := eva.cff(0);
        DoblDobl_CSeries_Polynomials.Clear(dp);
        err := pj(i,j) - jm(i,j); val := AbsVal(err);
        errsum := errsum + val;
      end loop;
      DoblDobl_CSeries_Polynomials.Clear(p);
    end loop;
    put("Sum of errors : "); put(errsum,3); new_line;
  end DoblDobl_Test_System;

  procedure QuadDobl_Test_System ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate random circuits for testing in quad double precision.
  --   Computes the Jacobian for the circuit
  --   and compares with the evaluation at the symbolic derivative.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Series;
    use QuadDobl_Complex_Series_Vectors;
    use QuadDobl_CSeries_Polynomials;
    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits(1..dim)
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    jm,pj : Matrix(1..dim,1..dim);
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    p,dp : Poly;
    sx : QuadDobl_Complex_Series_Vectors.Vector(1..dim);
    eva : QuadDobl_Complex_Series.Link_to_Series;
    err : Complex_Number;
    val,errsum : quad_double := create(0.0);

  begin
    for k in sx'range loop
      sx(k) := new Series'(Create(x(k),deg));
    end loop;
    jm := Jacobian_Convolution_Circuits.Jacobian(c,x);
    for i in 1..dim loop
      p := QuadDobl_Polynomial(c(i).all);
      for j in 1..dim loop
        dp := Diff(p,j);
        eva := QuadDobl_CSeries_Poly_Functions.Eval(dp,sx);
        pj(i,j) := eva.cff(0);
        QuadDobl_CSeries_Polynomials.Clear(dp);
        err := pj(i,j) - jm(i,j); val := AbsVal(err);
        errsum := errsum + val;
      end loop;
      QuadDobl_CSeries_Polynomials.Clear(p);
    end loop;
    put("Sum of errors : "); put(errsum,3); new_line;
  end QuadDobl_Test_System;

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
    new_line;
    case precision is
      when '0' => Standard_Test_System(dim,deg,nbr,pwr);
      when '1' => DoblDobl_Test_System(dim,deg,nbr,pwr);
      when '2' => QuadDobl_Test_System(dim,deg,nbr,pwr);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_jacocnv;
