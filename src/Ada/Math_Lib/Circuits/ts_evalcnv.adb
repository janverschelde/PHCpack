with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
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
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;

procedure ts_evalcnv is

-- DESCRIPTION :
--   Basic test at the plain evaluation of convolution circuits
--   at some random point.

  procedure Standard_Test ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Makes a system s of convolution circuits of p, generates a random
  --   vector x and then evaluates both p and s at x.  
  --   Finally the difference of both evaluations is computed,
  --   in standard double precision.

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
    z : Standard_Complex_Vectors.Vector(p'range) := Eval(c,x);
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
    Set_Solution_Constant(c,x);
    z := Eval(c,x);
    put_line("Evaluation after setting the solution constant :");
    put_line(z);
  end Standard_Test;

  procedure DoblDobl_Test ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Makes a system s of convolution circuits of p, generates a random
  --   vector x and then evaluates both p and s at x.  
  --   Finally the difference of both evaluations is computed,
  --   all in double double precision.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Speelpenning_Convolutions;

    dim : constant natural32
        := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    x : constant DoblDobl_Complex_Vectors.Vector(1..integer32(dim))
      := DoblDobl_Random_Vectors.Random_Vector(1,integer32(dim));
    y : constant DoblDobl_Complex_Vectors.Vector(p'range)
      := DoblDobl_Complex_Poly_SysFun.Eval(p,x);
    c : constant Convolution_Circuits(p'range)
      := Make_Convolution_Circuits(p,0);
    z : DoblDobl_Complex_Vectors.Vector(p'range) := Eval(c,x);
    val : Complex_Number;
    err : double_double := create(0.0);

  begin
    put("Evaluation at "); put(dim,1); put_line(" random numbers ...");
    put_line(x);
    put_line("The value at the system : ");
    put_line(y);
    put_line("The value at the convolution circuits : ");
    put_line(z);
    for k in y'range loop
      val := y(k) - z(k);
      err := err + DoblDobl_Complex_Numbers.AbsVal(val);
    end loop;
    put("The sum of errors : "); put(err,2); new_line;
    Set_Solution_Constant(c,x);
    z := Eval(c,x);
    put_line("Evaluation after setting the solution constant :");
    put_line(z);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Makes a system s of convolution circuits of p, generates a random
  --   vector x and then evaluates both p and s at x.  
  --   Finally the difference of both evaluations is computed,
  --   all in double double precision.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Speelpenning_Convolutions;

    dim : constant natural32
        := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    x : constant QuadDobl_Complex_Vectors.Vector(1..integer32(dim))
      := QuadDobl_Random_Vectors.Random_Vector(1,integer32(dim));
    y : constant QuadDobl_Complex_Vectors.Vector(p'range)
      := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    c : constant Convolution_Circuits(p'range)
      := Make_Convolution_Circuits(p,0);
    z : QuadDobl_Complex_Vectors.Vector(p'range) := Eval(c,x);
    val : Complex_Number;
    err : quad_double := create(0.0);

  begin
    put("Evaluation at "); put(dim,1); put_line(" random numbers ...");
    put_line(x);
    put_line("The value at the system : ");
    put_line(y);
    put_line("The value at the convolution circuits : ");
    put_line(z);
    for k in y'range loop
      val := y(k) - z(k);
      err := err + QuadDobl_Complex_Numbers.AbsVal(val);
    end loop;
    put("The sum of errors : "); put(err,2); new_line;
    Set_Solution_Constant(c,x);
    z := Eval(c,x);
    put_line("Evaluation after setting the solution constant :");
    put_line(z);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system p,
  --   and then launches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    put_line("Reading a polynomial system ...");
    case ans is 
      when '0' =>
        declare
          lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
        begin
          get(lp);
          Standard_Test(lp.all);
        end;
      when '1' =>
        declare
          lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
        begin
          get(lp);
          DoblDobl_Test(lp.all);
        end;
      when '2' =>
        declare
          lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
        begin
          get(lp);
          QuadDobl_Test(lp.all);
        end;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_evalcnv;
