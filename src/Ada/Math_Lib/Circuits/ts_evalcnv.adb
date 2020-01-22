with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;

procedure ts_evalcnv is

-- DESCRIPTION :
--   Basic test at the plain evaluation of convolution circuits
--   at some random point.

  procedure Test ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Makes a system s of convolution circuits of p, generates a random
  --   vector x and then evaluates both p and s at x.  
  --   Finally the difference of both evaluations is computed.

    use Standard_Complex_Numbers;
    use Standard_Speelpenning_Convolutions;

    dim : constant natural32
        := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    x : constant Standard_Complex_Vectors.Vector(1..integer32(dim))
      := Standard_Random_Vectors.Random_Vector(1,integer32(dim));
    y : constant Standard_Complex_Vectors.Vector(p'range)
      := Standard_Complex_Poly_SysFun.Eval(p,x);
    c : constant Convolution_Circuits(p'range)
      := Make_Convolution_Circuits(p,0);
    z : constant Standard_Complex_Vectors.Vector(p'range) := Eval(c,x);
    val : Complex_Number;
    err : double_float := 0.0;

  begin
    put("Evaluation at "); put(dim,1); put_line(" random numbers ...");
    put_line(x);
    put_line("The value at the system : ");
    put_line(y);
    put_line("The value at the convolution circuits : ");
    put_line(z);
    for k in y'range loop
      val := y(k) - z(k);
      err := err + Standard_Complex_Numbers.AbsVal(val);
    end loop;
    put("The sum of errors :"); put(err,2); new_line;
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system p,
  --   and then launches the test.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    Test(lp.all);
  end Main;

begin
  Main;
end ts_evalcnv;
