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
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Standard_Complex_Series;
with Standard_Complex_Series_io;         use Standard_Complex_Series_io;
with Standard_Complex_Series_Functions;
with Standard_Complex_Random_Series;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_io;         use DoblDobl_Complex_Series_io;
with DoblDobl_Complex_Series_Functions;
with DoblDobl_Complex_Random_Series;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_io;         use QuadDobl_Complex_Series_io;
with QuadDobl_Complex_Series_Functions;
with QuadDobl_Complex_Random_Series;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Shift_Convolution_Circuits;         use Shift_Convolution_Circuits;

procedure ts_shiftcnv is

-- DESCRIPTION :
--   Tests the development of the procedures to shift convolution circuits.

  procedure Standard_Test ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the shifting of the series parameter on a random series
  --   of the given degree, in double precision.

    use Standard_Complex_Numbers;
    use Standard_Complex_Series;
    use Standard_Complex_Series_Functions;

    s : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    rc : double_float := 0.0;
    shifteds : Series(degree);
    cc,y,z : Complex_Number;
    cff,wrk : Standard_Complex_Vectors.Vector(0..degree);

  begin
    put_line("Shifting the series parameter ...");
    put_line("on a random series s :"); put(s);
    put("Give a real constant for the shift : "); get(rc);
    shifteds := Shift(s,rc);
    y := Eval(s,-rc);
    z := Eval(shifteds,0.0);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("The shifted coefficients :"); put(shifteds);
    cff := s.cff;
    Shift(cff,wrk,rc);
    put_line("The shifted coefficient vector :"); put_line(cff);
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,-cc);
    z := Eval(shifteds,0.0);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    put_line("The shifted coefficients :"); put(shifteds);
    cff := s.cff;
    Shift(cff,wrk,cc);
    put_line("The shifted coefficient vector :"); put_line(cff);
  end Standard_Test;

  procedure DoblDobl_Test ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the shifting of the series parameter on a random series
  --   of the given degree, in double double precision.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Series;
    use DoblDobl_Complex_Series_Functions;

    s : constant Series(degree)
      := DoblDobl_Complex_Random_Series.Random_Series(degree);
    zero : constant double_double := create(0.0);
    rc : double_double := zero;
    shifteds : Series(degree);
    cc,y,z : Complex_Number;
    cff,wrk : DoblDobl_Complex_Vectors.Vector(0..degree);

  begin
    put_line("Shifting the series parameter ...");
    put_line("on a random series s :"); put(s);
    put("Give a real constant for the shift : "); get(rc);
    shifteds := Shift(s,rc);
    y := Eval(s,-rc);
    z := Eval(shifteds,zero);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("The shifted coefficients :"); put_line(shifteds.cff);
    cff := s.cff;
    Shift(cff,wrk,rc);
    put_line("The shifted coefficient vector :"); put_line(cff);
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,-cc);
    z := Eval(shifteds,zero);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("The shifted coefficients :"); put_line(shifteds.cff);
    cff := s.cff;
    Shift(cff,wrk,cc);
    put_line("The shifted coefficient vector :"); put_line(cff);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the shifting of the series parameter on a random series
  --   of the given degree, in quaddouble precision.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Series;
    use QuadDobl_Complex_Series_Functions;

    s : constant Series(degree)
      := QuadDobl_Complex_Random_Series.Random_Series(degree);
    zero : constant quad_double := create(0.0);
    rc : quad_double := zero;
    shifteds : Series(degree);
    cc,y,z : Complex_Number;
    cff,wrk : QuadDobl_Complex_Vectors.Vector(0..degree);

  begin
    put_line("Shifting the series parameter ...");
    put_line("on a random series s :"); put(s);
    put("Give a real constant for the shift : "); get(rc);
    shifteds := Shift(s,rc);
    y := Eval(s,-rc);
    z := Eval(shifteds,zero);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("The shifted coefficients :"); put_line(shifteds.cff); 
    cff := s.cff;
    Shift(cff,wrk,rc);
    put_line("The shifted coefficient vector :"); put_line(cff);
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,-cc);
    z := Eval(shifteds,zero);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("The shifted coefficients :"); put_line(shifteds.cff);
    Shift(cff,wrk,cc);
    put_line("The shifted coefficient vector :"); put_line(cff);
  end QuadDobl_Test;

  procedure Standard_System_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of exponents and coefficients.
  --   Run tests in double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Complex_Numbers;
    use Standard_Complex_Series;
    use Standard_Complex_Series_Functions;
    use Standard_Speelpenning_Convolutions;

    c : constant Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : constant Link_to_System := Create(c,dim,deg);
    rc : double_float := 0.0;
    zero : constant Complex_Number := Create(0.0);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector'(0..deg => zero);
    y,z : Complex_Number;
    scst : Series(deg);
    xpt : constant Standard_Complex_Vectors.Vector(1..s.dim)
        := Standard_Random_Vectors.Random_Vector(1,s.dim);
    sy,sz : Standard_Complex_Vectors.Vector(1..s.neq);
    ans : character;

  begin
    new_line;
    put("Give a real constant for the shift : "); get(rc);
    put("Test shifting of a series ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      scst.cff := s.crc(1).cst.all;
      y := Eval(scst,-rc);
      scst.cff := s.crc(1).cst.all;
      z := Eval(scst,zero);
      put("s(-shift constant) : "); put(y); new_line;
      put(" shifted series(0) : "); put(z); new_line;
    else -- testing shifting the entire system
      sy := Eval(s.crc,xpt,Standard_Complex_Numbers.Create(-rc));
      Shift(s,wrk,rc);
      sz := Eval(s.crc,xpt,zero);
      put_line("s(-shift constant) : "); put_line(sy); new_line;
      put_line(" shifted system(0) : "); put_line(sz); new_line;
    end if;
  end Standard_System_Test;

  procedure DoblDobl_System_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of exponents and coefficients.
  --   Run tests in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Series;
    use DoblDobl_Complex_Series_Functions;
    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : constant Link_to_System := Create(c,dim,deg);
    rc : double_double := create(0.0);
    zero : constant Complex_Number := Create(integer(0));
    wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector'(0..deg => zero);
    y,z : Complex_Number;
    xpt : constant DoblDobl_Complex_Vectors.Vector(1..s.dim)
        := DoblDobl_Random_Vectors.Random_Vector(1,s.dim);
    sy,sz : DoblDobl_Complex_Vectors.Vector(1..s.neq);
    scst : Series(deg);
    ans : character;

  begin
    new_line;
    put("Give a real constant for the shift : "); get(rc);
    put("Test shifting of a series ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      scst.cff := s.crc(1).cst.all;
      y := Eval(scst,-rc);
      Shift(s,wrk,rc);
      scst.cff := s.crc(1).cst.all;
      z := Eval(scst,zero);
      put("s(-shift constant) : "); put(y); new_line;
      put(" shifted series(0) : "); put(z); new_line;
    else -- testing shifting the entire system
      sy := Eval(s.crc,xpt,DoblDobl_Complex_Numbers.Create(-rc));
      Shift(s,wrk,rc);
      sz := Eval(s.crc,xpt,zero);
      put_line("s(-shift constant) : "); put_line(sy); new_line;
      put_line(" shifted system(0) : "); put_line(sz); new_line;
    end if;
  end DoblDobl_System_Test;

  procedure QuadDobl_System_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of exponents and coefficients.
  --   Run tests in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Series;
    use QuadDobl_Complex_Series_Functions;
    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : constant Link_to_System := Create(c,dim,deg);
    rc : quad_double := create(0.0);
    zero : constant Complex_Number := Create(integer(0));
    wrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector'(0..deg => zero);
    y,z : Complex_Number;
    xpt : constant QuadDobl_Complex_Vectors.Vector(1..s.dim)
        := QuadDobl_Random_Vectors.Random_Vector(1,s.dim);
    sy,sz : QuadDobl_Complex_Vectors.Vector(1..s.neq);
    scst : Series(deg);
    ans : character;

  begin
    new_line;
    put("Give a real constant for the shift : "); get(rc);
    put("Test shifting of a series ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      scst.cff := s.crc(1).cst.all;
      y := Eval(scst,-rc);
      Shift(s,wrk,rc);
      scst.cff := s.crc(1).cst.all;
      z := Eval(scst,zero);
      put("s(-shift constant) : "); put(y); new_line;
      put(" shifted series(0) : "); put(z); new_line;
    else -- testing shifting the entire system
      sy := Eval(s.crc,xpt,QuadDobl_Complex_Numbers.Create(-rc));
      Shift(s,wrk,rc);
      sz := Eval(s.crc,xpt,zero);
      put_line("s(-shift constant) : "); put_line(sy); new_line;
      put_line(" shifted system(0) : "); put_line(sz); new_line;
    end if;
  end QuadDobl_System_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a degree and then launches the test.

    deg,dim,nbr,pwr : integer32 := 0;
    systest,precision : character;

  begin
    new_line;
    put_line("Testing the shifting of coefficients of circuits ...");
    new_line;
    put("Give the degree : "); get(deg);
    new_line;
    put("Test system ? (y/n) "); Ask_Yes_or_No(systest);
    if systest = 'y' then
      put("Give the dimension of system : "); get(dim);
      put("Give the number of terms : "); get(nbr);
      put("Give the highest power : "); get(pwr);
    end if;
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    if systest = 'y' then
      case precision is
        when '0' => Standard_System_Test(dim,deg,nbr,pwr);
        when '1' => DoblDobl_System_Test(dim,deg,nbr,pwr);
        when '2' => QuadDobl_System_Test(dim,deg,nbr,pwr);
        when others => null;
      end case;
    else
      case precision is
        when '0' => Standard_Test(deg);
        when '1' => DoblDobl_Test(deg);
        when '2' => QuadDobl_Test(deg);
        when others => null;
      end case;
    end if;
  end Main;


begin
  Main;
end ts_shiftcnv;
