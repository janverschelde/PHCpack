with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with Standard_Integer_Vectors;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;

procedure ts_ratapp is

-- DESCRIPTION :
--    This test procedure checks the results of reorganized construction
--    of a vector of Pade approximants, comparing with the previous code.

  procedure Standard_Test ( numdeg,dendeg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Runs a test on a random coefficient vector of series,
  --   in double precision.

    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := Standard_Random_Vectors.Random_Vector(0,dim);
    numcff1 : Standard_Complex_Vectors.Vector(0..numdeg);
    numcff2 : Standard_Complex_Vectors.Vector(0..numdeg);
    dencff1 : Standard_Complex_Vectors.Vector(0..dendeg);
    dencff2 : Standard_Complex_Vectors.Vector(0..dendeg);
    mat : Standard_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : Standard_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use Standard_Rational_Approximations;

  begin
    Pade(numdeg,dendeg,cff,numcff1,dencff1,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff1);
    put_line("The coefficients of the denominator :"); put_line(dencff1);
    Pade(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff2);
    put_line("The coefficients of the denominator :"); put_line(dencff2);
  end Standard_Test;

  procedure DoblDobl_Test ( numdeg,dendeg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Runs a test on a random coefficient vector of series,
  --   in double double precision.

    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := DoblDobl_Random_Vectors.Random_Vector(0,dim);
    numcff1 : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    numcff2 : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    dencff1 : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    dencff2 : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    mat : DoblDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : DoblDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use DoblDobl_Rational_Approximations;

  begin
    Pade(numdeg,dendeg,cff,numcff1,dencff1,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff1);
    put_line("The coefficients of the denominator :"); put_line(dencff1);
    Pade(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff2);
    put_line("The coefficients of the denominator :"); put_line(dencff2);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( numdeg,dendeg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Runs a test on a random coefficient vector of series,
  --   in quad double precision.

    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := QuadDobl_Random_Vectors.Random_Vector(0,dim);
    numcff1 : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    numcff2 : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    dencff1 : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    dencff2 : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    mat : QuadDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use QuadDobl_Rational_Approximations;

  begin
    Pade(numdeg,dendeg,cff,numcff1,dencff1,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff1);
    put_line("The coefficients of the denominator :"); put_line(dencff1);
    Pade(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff2);
    put_line("The coefficients of the denominator :"); put_line(dencff2);
  end QuadDobl_Test;

  procedure Standard_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Runs a test on a vector of range 1..nbr of random coefficient 
  --   vectors of series, in double precision.

    cff : Standard_Complex_VecVecs.VecVec(1..nbr);
    numcff1 : Standard_Complex_VecVecs.VecVec(0..numdeg);
    numcff2 : Standard_Complex_VecVecs.VecVec(0..numdeg);
    dencff1 : Standard_Complex_VecVecs.VecVec(0..dendeg);
    dencff2 : Standard_Complex_VecVecs.VecVec(0..dendeg);
    mat : Standard_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : Standard_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use Standard_Rational_Approximations;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant Standard_Complex_Vectors.Vector(0..dim)
             := Standard_Random_Vectors.Random_Vector(0,dim);
      begin
        cff(k) := new Standard_Complex_Vectors.Vector'(kcff);
      end;
    end loop;
    for k in 1..nbr loop
      declare
        knum : Standard_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => Standard_Complex_Numbers.Create(integer(0)));
        kden : Standard_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => Standard_Complex_Numbers.Create(integer(0)));
        kcff : constant Standard_Complex_Vectors.Link_to_Vector := cff(k);
      begin
        Pade(numdeg,dendeg,kcff.all,knum,kden,info,true);
        numcff1(k) := new Standard_Complex_Vectors.Vector'(knum);
        dencff1(k) := new Standard_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    put_line("The coefficients of the numerators :"); put_line(numcff1);
    put_line("The coefficients of the denominators :"); put_line(dencff1);
    for k in 1..nbr loop
      declare
        knum : constant Standard_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => Standard_Complex_Numbers.Create(integer(0)));
        kden : constant Standard_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => Standard_Complex_Numbers.Create(integer(0)));
      begin
        numcff2(k) := new Standard_Complex_Vectors.Vector'(knum);
        dencff2(k) := new Standard_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    Pade_Vector(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerators :"); put_line(numcff2);
    put_line("The coefficients of the denominators :"); put_line(dencff2);
  end Standard_Test_Vector;

  procedure DoblDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Runs a test on a vector of range 1..nbr of random coefficient 
  --   vectors of series, in double double precision.

    cff : DoblDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1 : DoblDobl_Complex_VecVecs.VecVec(0..numdeg);
    numcff2 : DoblDobl_Complex_VecVecs.VecVec(0..numdeg);
    dencff1 : DoblDobl_Complex_VecVecs.VecVec(0..dendeg);
    dencff2 : DoblDobl_Complex_VecVecs.VecVec(0..dendeg);
    mat : DoblDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : DoblDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use DoblDobl_Rational_Approximations;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
             := DoblDobl_Random_Vectors.Random_Vector(0,dim);
      begin
        cff(k) := new DoblDobl_Complex_Vectors.Vector'(kcff);
      end;
    end loop;
    for k in 1..nbr loop
      declare
        knum : DoblDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => DoblDobl_Complex_Numbers.Create(integer(0)));
        kden : DoblDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => DoblDobl_Complex_Numbers.Create(integer(0)));
        kcff : constant DoblDobl_Complex_Vectors.Link_to_Vector := cff(k);
      begin
        Pade(numdeg,dendeg,kcff.all,knum,kden,info,true);
        numcff1(k) := new DoblDobl_Complex_Vectors.Vector'(knum);
        dencff1(k) := new DoblDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    put_line("The coefficients of the numerators :"); put_line(numcff1);
    put_line("The coefficients of the denominators :"); put_line(dencff1);
    for k in 1..nbr loop
      declare
        knum : constant DoblDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => DoblDobl_Complex_Numbers.Create(integer(0)));
        kden : constant DoblDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => DoblDobl_Complex_Numbers.Create(integer(0)));
      begin
        numcff2(k) := new DoblDobl_Complex_Vectors.Vector'(knum);
        dencff2(k) := new DoblDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    Pade_Vector(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerators :"); put_line(numcff2);
    put_line("The coefficients of the denominators :"); put_line(dencff2);
  end DoblDobl_Test_Vector;

  procedure QuadDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Runs a test on a vector of range 1..nbr of random coefficient 
  --   vectors of series, in quad double precision.

    cff : QuadDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1 : QuadDobl_Complex_VecVecs.VecVec(0..numdeg);
    numcff2 : QuadDobl_Complex_VecVecs.VecVec(0..numdeg);
    dencff1 : QuadDobl_Complex_VecVecs.VecVec(0..dendeg);
    dencff2 : QuadDobl_Complex_VecVecs.VecVec(0..dendeg);
    mat : QuadDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use QuadDobl_Rational_Approximations;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
             := QuadDobl_Random_Vectors.Random_Vector(0,dim);
      begin
        cff(k) := new QuadDobl_Complex_Vectors.Vector'(kcff);
      end;
    end loop;
    for k in 1..nbr loop
      declare
        knum : QuadDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => QuadDobl_Complex_Numbers.Create(integer(0)));
        kden : QuadDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => QuadDobl_Complex_Numbers.Create(integer(0)));
        kcff : constant QuadDobl_Complex_Vectors.Link_to_Vector := cff(k);
      begin
        Pade(numdeg,dendeg,kcff.all,knum,kden,info,true);
        numcff1(k) := new QuadDobl_Complex_Vectors.Vector'(knum);
        dencff1(k) := new QuadDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    put_line("The coefficients of the numerators :"); put_line(numcff1);
    put_line("The coefficients of the denominators :"); put_line(dencff1);
    for k in 1..nbr loop
      declare
        knum : constant QuadDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => QuadDobl_Complex_Numbers.Create(integer(0)));
        kden : constant QuadDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => QuadDobl_Complex_Numbers.Create(integer(0)));
      begin
        numcff2(k) := new QuadDobl_Complex_Vectors.Vector'(knum);
        dencff2(k) := new QuadDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    Pade_Vector(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerators :"); put_line(numcff2);
    put_line("The coefficients of the denominators :"); put_line(dencff2);
  end QuadDobl_Test_Vector;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the problem
  --   and then launches the test.

    numdeg,dendeg,dim,nbr : integer32 := 0;
    ans,precision : character;

  begin
    new_line;
    put_line("Testing the construction of rational approximations ...");
    new_line;
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    dim := numdeg+dendeg;
    new_line;
    precision := Prompt_for_Precision;
    new_line;
    put("Test vector ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then
      case precision is
        when '0' => Standard_Test(numdeg,dendeg,dim);
        when '1' => DoblDobl_Test(numdeg,dendeg,dim);
        when '2' => QuadDobl_Test(numdeg,dendeg,dim);
        when others => null;
      end case;
    else
      new_line;
      put("Give the number of components : "); get(nbr);
      case precision is
        when '0' => Standard_Test_Vector(nbr,numdeg,dendeg,dim);
        when '1' => DoblDobl_Test_Vector(nbr,numdeg,dendeg,dim);
        when '2' => QuadDobl_Test_Vector(nbr,numdeg,dendeg,dim);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_ratapp;
