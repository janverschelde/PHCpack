with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
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
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Newton_Convolutions;
with Newton_Power_Convolutions;          use Newton_Power_Convolutions;

procedure ts_sernewcnv is

-- DESCRIPTION :
--   Procedure to develop the linearized Newton's method for power series,
--   on convolution circuits.

  procedure Add_Parameter_to_Constant
              ( c : in Standard_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 ) is

  -- DESCRIPTION :
  --   Adds the continuation parameter t to the constant of c.

  -- REQUIRED : the coefficients in the power series of c
  --   have degree at least deg.

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(1.0);
    else
      c.cst := new Standard_Complex_Vectors.Vector'(0..deg => Create(0.0));
      c.cst(1) := Create(1.0);
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 ) is

  -- DESCRIPTION :
  --   Adds the continuation parameter t to the constant of c.

  -- REQUIRED : the coefficients in the power series of c
  --   have degree at least deg.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new DoblDobl_Complex_Vectors.Vector'
                     (0..deg => Create(integer32(0)));
      c.cst(1) := Create(integer32(1));
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 ) is

  -- DESCRIPTION :
  --   Adds the continuation parameter t to the constant of c.

  -- REQUIRED : the coefficients in the power series of c
  --   have degree at least deg.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new QuadDobl_Complex_Vectors.Vector'
                     (0..deg => Create(integer32(0)));
      c.cst(1) := Create(integer32(1));
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System ) is

  -- DESCRIPTION :
  --   Adds the continuation parameter to every circuit in s.

  -- REQUIRED : s /= null.

  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System ) is

  -- DESCRIPTION :
  --   Adds the continuation parameter to every circuit in s.

  -- REQUIRED : s /= null.

  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System ) is

  -- DESCRIPTION :
  --   Adds the continuation parameter to every circuit in s.

  -- REQUIRED : s /= null.

  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Standard_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Given the degree of the power series, prompts the user
  --   for a polynomial system with a parameter,
  --   runs Newton's method in standard double precision.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim,nbr : natural32;

    use Standard_Speelpenning_Convolutions;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(lp,sols);
    dim := natural32(lp'last);
    nbr := Standard_Complex_Solutions.Length_Of(sols);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    declare
      c : constant Convolution_Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      s : constant Link_to_System := Create(c,lp'last,deg);
      sol : constant Standard_Complex_Solutions.Link_to_Solution
          := Standard_Complex_Solutions.Head_Of(sols);
      scf : constant Standard_Complex_VecVecs.VecVec(1..sol.n)
          := Newton_Convolutions.Series_Coefficients(sol.v,deg);
      info,nbrit,maxit : integer32 := 0;
      ipvt : Standard_Integer_Vectors.Vector(1..sol.n);
      wrk : Standard_Complex_Vectors.Link_to_Vector
          := new Standard_Complex_Vectors.Vector(1..sol.n);
      absdx : double_float;
      tol : constant double_float := 1.0E-14;
      fail : boolean;
    begin
      Add_Parameter_to_Constant(s);
      new_line;
      put("Give the number of iterations : "); get(maxit);
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk);
      if fail then
        put("Failed to reach"); put(tol,3);
        put("  absdx : "); put(absdx,3); new_line;
      else
        put("Reached"); put(tol,3);
        put("  absdx : "); put(absdx,3); new_line;
      end if;
      put("after "); put(nbrit,1); put_line(" iterations.");
      Standard_Complex_Vectors.Clear(wrk);
    end;
  end Standard_Test;

  procedure DoblDobl_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Given the degree of the power series, prompts the user
  --   for a polynomial system with a parameter,
  --   runs Newton's method in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    dim,nbr : natural32;

    use DoblDobl_Speelpenning_Convolutions;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    dim := natural32(lp'last);
    nbr := DoblDobl_Complex_Solutions.Length_Of(sols);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    declare
      c : constant Convolution_Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      s : constant Link_to_System := Create(c,lp'last,deg);
      sol : constant DoblDobl_Complex_Solutions.Link_to_Solution
          := DoblDobl_Complex_Solutions.Head_Of(sols);
      scf : constant DoblDobl_Complex_VecVecs.VecVec(1..sol.n)
          := Newton_Convolutions.Series_Coefficients(sol.v,deg);
      info,nbrit,maxit : integer32 := 0;
      ipvt : Standard_Integer_Vectors.Vector(1..sol.n);
      wrk : DoblDobl_Complex_Vectors.Link_to_Vector
          := new DoblDobl_Complex_Vectors.Vector(1..sol.n);
      absdx : double_double;
      tol : constant double_float := 1.0E-14;
      fail : boolean;
    begin
      Add_Parameter_to_Constant(s);
      new_line;
      put("Give the number of iterations : "); get(maxit);
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk);
      if fail then
        put("Failed to reach"); put(tol,3);
        put("  absdx : "); put(absdx,3); new_line;
      else
        put("Reached"); put(tol,3);
        put("  absdx : "); put(absdx,3); new_line;
      end if;
      put("after "); put(nbrit,1); put_line(" iterations.");
      DoblDobl_Complex_Vectors.Clear(wrk);
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Given the degree of the power series, prompts the user
  --   for a polynomial system with a parameter,
  --   runs Newton's method in double double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    dim,nbr : natural32;

    use QuadDobl_Speelpenning_Convolutions;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    dim := natural32(lp'last);
    nbr := QuadDobl_Complex_Solutions.Length_Of(sols);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    declare
      c : constant Convolution_Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      s : constant Link_to_System := Create(c,lp'last,deg);
      sol : constant QuadDobl_Complex_Solutions.Link_to_Solution
          := QuadDobl_Complex_Solutions.Head_Of(sols);
      scf : constant QuadDobl_Complex_VecVecs.VecVec(1..sol.n)
          := Newton_Convolutions.Series_Coefficients(sol.v,deg);
      info,nbrit,maxit : integer32 := 0;
      ipvt : Standard_Integer_Vectors.Vector(1..sol.n);
      wrk : QuadDobl_Complex_Vectors.Link_to_Vector
          := new QuadDobl_Complex_Vectors.Vector(1..sol.n);
      absdx : quad_double;
      tol : constant double_float := 1.0E-14;
      fail : boolean;
    begin
      Add_Parameter_to_Constant(s);
      new_line;
      put("Give the number of iterations : "); get(maxit);
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk);
      if fail then
        put("Failed to reach"); put(tol,3);
        put("  absdx : "); put(absdx,3); new_line;
      else
        put("Reached"); put(tol,3);
        put("  absdx : "); put(absdx,3); new_line;
      end if;
      put("after "); put(nbrit,1); put_line(" iterations.");
      QuadDobl_Complex_Vectors.Clear(wrk);
    end;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree of the power series
  --   and then launches the test.

    deg : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Linearized Newton on power series with convolution circuits.");
    new_line;
    put("Give the degree of the power series : "); get(deg); skip_line;
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test(deg);
      when '1' => DoblDobl_Test(deg);
      when '2' => QuadDobl_Test(deg);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sernewcnv;
