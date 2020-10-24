with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Vectors_io;        use TripDobl_Complex_Vectors_io;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_VecVecs_io;        use TripDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Vectors_io;        use PentDobl_Complex_Vectors_io;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_VecVecs_io;        use PentDobl_Complex_VecVecs_io;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Vectors_io;        use OctoDobl_Complex_Vectors_io;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_VecVecs_io;        use OctoDobl_Complex_VecVecs_io;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Vectors_io;        use DecaDobl_Complex_Vectors_io;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_VecVecs_io;        use DecaDobl_Complex_VecVecs_io;
with Standard_Integer_Vectors;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with TripDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with PentDobl_Random_Vectors;
with OctoDobl_Random_Vectors;
with DecaDobl_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with TripDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with PentDobl_Complex_Matrices;
with OctoDobl_Complex_Matrices;
with DecaDobl_Complex_Matrices;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with TripDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with PentDobl_Rational_Approximations;
with OctoDobl_Rational_Approximations;
with DecaDobl_Rational_Approximations;

package body Test_Rational_Approximations is

  procedure Standard_Test ( numdeg,dendeg,dim : in integer32 ) is

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

  procedure TripDobl_Test ( numdeg,dendeg,dim : in integer32 ) is

    cff : constant TripDobl_Complex_Vectors.Vector(0..dim)
        := TripDobl_Random_Vectors.Random_Vector(0,dim);
    numcff1 : TripDobl_Complex_Vectors.Vector(0..numdeg);
    numcff2 : TripDobl_Complex_Vectors.Vector(0..numdeg);
    dencff1 : TripDobl_Complex_Vectors.Vector(0..dendeg);
    dencff2 : TripDobl_Complex_Vectors.Vector(0..dendeg);
    mat : TripDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : TripDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use TripDobl_Rational_Approximations;

  begin
    Pade(numdeg,dendeg,cff,numcff1,dencff1,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff1);
    put_line("The coefficients of the denominator :"); put_line(dencff1);
    Pade(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff2);
    put_line("The coefficients of the denominator :"); put_line(dencff2);
  end TripDobl_Test;

  procedure QuadDobl_Test ( numdeg,dendeg,dim : in integer32 ) is

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

  procedure PentDobl_Test ( numdeg,dendeg,dim : in integer32 ) is

    cff : constant PentDobl_Complex_Vectors.Vector(0..dim)
        := PentDobl_Random_Vectors.Random_Vector(0,dim);
    numcff1 : PentDobl_Complex_Vectors.Vector(0..numdeg);
    numcff2 : PentDobl_Complex_Vectors.Vector(0..numdeg);
    dencff1 : PentDobl_Complex_Vectors.Vector(0..dendeg);
    dencff2 : PentDobl_Complex_Vectors.Vector(0..dendeg);
    mat : PentDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : PentDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use PentDobl_Rational_Approximations;

  begin
    Pade(numdeg,dendeg,cff,numcff1,dencff1,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff1);
    put_line("The coefficients of the denominator :"); put_line(dencff1);
    Pade(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff2);
    put_line("The coefficients of the denominator :"); put_line(dencff2);
  end PentDobl_Test;

  procedure OctoDobl_Test ( numdeg,dendeg,dim : in integer32 ) is

    cff : constant OctoDobl_Complex_Vectors.Vector(0..dim)
        := OctoDobl_Random_Vectors.Random_Vector(0,dim);
    numcff1 : OctoDobl_Complex_Vectors.Vector(0..numdeg);
    numcff2 : OctoDobl_Complex_Vectors.Vector(0..numdeg);
    dencff1 : OctoDobl_Complex_Vectors.Vector(0..dendeg);
    dencff2 : OctoDobl_Complex_Vectors.Vector(0..dendeg);
    mat : OctoDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : OctoDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use OctoDobl_Rational_Approximations;

  begin
    Pade(numdeg,dendeg,cff,numcff1,dencff1,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff1);
    put_line("The coefficients of the denominator :"); put_line(dencff1);
    Pade(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff2);
    put_line("The coefficients of the denominator :"); put_line(dencff2);
  end OctoDobl_Test;

  procedure DecaDobl_Test ( numdeg,dendeg,dim : in integer32 ) is

    cff : constant DecaDobl_Complex_Vectors.Vector(0..dim)
        := DecaDobl_Random_Vectors.Random_Vector(0,dim);
    numcff1 : DecaDobl_Complex_Vectors.Vector(0..numdeg);
    numcff2 : DecaDobl_Complex_Vectors.Vector(0..numdeg);
    dencff1 : DecaDobl_Complex_Vectors.Vector(0..dendeg);
    dencff2 : DecaDobl_Complex_Vectors.Vector(0..dendeg);
    mat : DecaDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : DecaDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use DecaDobl_Rational_Approximations;

  begin
    Pade(numdeg,dendeg,cff,numcff1,dencff1,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff1);
    put_line("The coefficients of the denominator :"); put_line(dencff1);
    Pade(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerator :"); put_line(numcff2);
    put_line("The coefficients of the denominator :"); put_line(dencff2);
  end DecaDobl_Test;

  procedure Standard_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

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

  procedure TripDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

    cff : TripDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1 : TripDobl_Complex_VecVecs.VecVec(0..numdeg);
    numcff2 : TripDobl_Complex_VecVecs.VecVec(0..numdeg);
    dencff1 : TripDobl_Complex_VecVecs.VecVec(0..dendeg);
    dencff2 : TripDobl_Complex_VecVecs.VecVec(0..dendeg);
    mat : TripDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : TripDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use TripDobl_Rational_Approximations;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant TripDobl_Complex_Vectors.Vector(0..dim)
             := TripDobl_Random_Vectors.Random_Vector(0,dim);
      begin
        cff(k) := new TripDobl_Complex_Vectors.Vector'(kcff);
      end;
    end loop;
    for k in 1..nbr loop
      declare
        knum : TripDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => TripDobl_Complex_Numbers.Create(integer(0)));
        kden : TripDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => TripDobl_Complex_Numbers.Create(integer(0)));
        kcff : constant TripDobl_Complex_Vectors.Link_to_Vector := cff(k);
      begin
        Pade(numdeg,dendeg,kcff.all,knum,kden,info,true);
        numcff1(k) := new TripDobl_Complex_Vectors.Vector'(knum);
        dencff1(k) := new TripDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    put_line("The coefficients of the numerators :"); put_line(numcff1);
    put_line("The coefficients of the denominators :"); put_line(dencff1);
    for k in 1..nbr loop
      declare
        knum : constant TripDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => TripDobl_Complex_Numbers.Create(integer(0)));
        kden : constant TripDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => TripDobl_Complex_Numbers.Create(integer(0)));
      begin
        numcff2(k) := new TripDobl_Complex_Vectors.Vector'(knum);
        dencff2(k) := new TripDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    Pade_Vector(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerators :"); put_line(numcff2);
    put_line("The coefficients of the denominators :"); put_line(dencff2);
  end TripDobl_Test_Vector;

  procedure QuadDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

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

  procedure PentDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

    cff : PentDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1 : PentDobl_Complex_VecVecs.VecVec(0..numdeg);
    numcff2 : PentDobl_Complex_VecVecs.VecVec(0..numdeg);
    dencff1 : PentDobl_Complex_VecVecs.VecVec(0..dendeg);
    dencff2 : PentDobl_Complex_VecVecs.VecVec(0..dendeg);
    mat : PentDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : PentDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use PentDobl_Rational_Approximations;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant PentDobl_Complex_Vectors.Vector(0..dim)
             := PentDobl_Random_Vectors.Random_Vector(0,dim);
      begin
        cff(k) := new PentDobl_Complex_Vectors.Vector'(kcff);
      end;
    end loop;
    for k in 1..nbr loop
      declare
        knum : PentDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => PentDobl_Complex_Numbers.Create(integer(0)));
        kden : PentDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => PentDobl_Complex_Numbers.Create(integer(0)));
        kcff : constant PentDobl_Complex_Vectors.Link_to_Vector := cff(k);
      begin
        Pade(numdeg,dendeg,kcff.all,knum,kden,info,true);
        numcff1(k) := new PentDobl_Complex_Vectors.Vector'(knum);
        dencff1(k) := new PentDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    put_line("The coefficients of the numerators :"); put_line(numcff1);
    put_line("The coefficients of the denominators :"); put_line(dencff1);
    for k in 1..nbr loop
      declare
        knum : constant PentDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => PentDobl_Complex_Numbers.Create(integer(0)));
        kden : constant PentDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => PentDobl_Complex_Numbers.Create(integer(0)));
      begin
        numcff2(k) := new PentDobl_Complex_Vectors.Vector'(knum);
        dencff2(k) := new PentDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    Pade_Vector(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerators :"); put_line(numcff2);
    put_line("The coefficients of the denominators :"); put_line(dencff2);
  end PentDobl_Test_Vector;

  procedure OctoDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

    cff : OctoDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1 : OctoDobl_Complex_VecVecs.VecVec(0..numdeg);
    numcff2 : OctoDobl_Complex_VecVecs.VecVec(0..numdeg);
    dencff1 : OctoDobl_Complex_VecVecs.VecVec(0..dendeg);
    dencff2 : OctoDobl_Complex_VecVecs.VecVec(0..dendeg);
    mat : OctoDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : OctoDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use OctoDobl_Rational_Approximations;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant OctoDobl_Complex_Vectors.Vector(0..dim)
             := OctoDobl_Random_Vectors.Random_Vector(0,dim);
      begin
        cff(k) := new OctoDobl_Complex_Vectors.Vector'(kcff);
      end;
    end loop;
    for k in 1..nbr loop
      declare
        knum : OctoDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => OctoDobl_Complex_Numbers.Create(integer(0)));
        kden : OctoDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => OctoDobl_Complex_Numbers.Create(integer(0)));
        kcff : constant OctoDobl_Complex_Vectors.Link_to_Vector := cff(k);
      begin
        Pade(numdeg,dendeg,kcff.all,knum,kden,info,true);
        numcff1(k) := new OctoDobl_Complex_Vectors.Vector'(knum);
        dencff1(k) := new OctoDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    put_line("The coefficients of the numerators :"); put_line(numcff1);
    put_line("The coefficients of the denominators :"); put_line(dencff1);
    for k in 1..nbr loop
      declare
        knum : constant OctoDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => OctoDobl_Complex_Numbers.Create(integer(0)));
        kden : constant OctoDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => OctoDobl_Complex_Numbers.Create(integer(0)));
      begin
        numcff2(k) := new OctoDobl_Complex_Vectors.Vector'(knum);
        dencff2(k) := new OctoDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    Pade_Vector(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerators :"); put_line(numcff2);
    put_line("The coefficients of the denominators :"); put_line(dencff2);
  end OctoDobl_Test_Vector;

  procedure DecaDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 ) is

    cff : DecaDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1 : DecaDobl_Complex_VecVecs.VecVec(0..numdeg);
    numcff2 : DecaDobl_Complex_VecVecs.VecVec(0..numdeg);
    dencff1 : DecaDobl_Complex_VecVecs.VecVec(0..dendeg);
    dencff2 : DecaDobl_Complex_VecVecs.VecVec(0..dendeg);
    mat : DecaDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : DecaDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;

    use DecaDobl_Rational_Approximations;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant DecaDobl_Complex_Vectors.Vector(0..dim)
             := DecaDobl_Random_Vectors.Random_Vector(0,dim);
      begin
        cff(k) := new DecaDobl_Complex_Vectors.Vector'(kcff);
      end;
    end loop;
    for k in 1..nbr loop
      declare
        knum : DecaDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => DecaDobl_Complex_Numbers.Create(integer(0)));
        kden : DecaDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => DecaDobl_Complex_Numbers.Create(integer(0)));
        kcff : constant DecaDobl_Complex_Vectors.Link_to_Vector := cff(k);
      begin
        Pade(numdeg,dendeg,kcff.all,knum,kden,info,true);
        numcff1(k) := new DecaDobl_Complex_Vectors.Vector'(knum);
        dencff1(k) := new DecaDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    put_line("The coefficients of the numerators :"); put_line(numcff1);
    put_line("The coefficients of the denominators :"); put_line(dencff1);
    for k in 1..nbr loop
      declare
        knum : constant DecaDobl_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => DecaDobl_Complex_Numbers.Create(integer(0)));
        kden : constant DecaDobl_Complex_Vectors.Vector(0..dendeg)
             := (0..numdeg => DecaDobl_Complex_Numbers.Create(integer(0)));
      begin
        numcff2(k) := new DecaDobl_Complex_Vectors.Vector'(knum);
        dencff2(k) := new DecaDobl_Complex_Vectors.Vector'(kden);
      end;
    end loop;
    Pade_Vector(numdeg,dendeg,cff,numcff2,dencff2,mat,rhs,ipvt,info,true);
    put_line("The coefficients of the numerators :"); put_line(numcff2);
    put_line("The coefficients of the denominators :"); put_line(dencff2);
  end DecaDobl_Test_Vector;

  procedure Main is

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
    put_line("MENU for the working precision :");
    put_line("  1. double precision");
    put_line("  2. double double precision");
    put_line("  3. triple double precision");
    put_line("  4. quad double precision");
    put_line("  5. penta double precision");
    put_line("  6. octo double precision");
    put_line("  7. deca double precision");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(precision,"1234567");
    new_line;
    put("Test vector ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then
      case precision is
        when '1' => Standard_Test(numdeg,dendeg,dim);
        when '2' => DoblDobl_Test(numdeg,dendeg,dim);
        when '3' => TripDobl_Test(numdeg,dendeg,dim);
        when '4' => QuadDobl_Test(numdeg,dendeg,dim);
        when '5' => PentDobl_Test(numdeg,dendeg,dim);
        when '6' => OctoDobl_Test(numdeg,dendeg,dim);
        when '7' => DecaDobl_Test(numdeg,dendeg,dim);
        when others => null;
      end case;
    else
      new_line;
      put("Give the number of components : "); get(nbr);
      case precision is
        when '1' => Standard_Test_Vector(nbr,numdeg,dendeg,dim);
        when '2' => DoblDobl_Test_Vector(nbr,numdeg,dendeg,dim);
        when '3' => TripDobl_Test_Vector(nbr,numdeg,dendeg,dim);
        when '4' => QuadDobl_Test_Vector(nbr,numdeg,dendeg,dim);
        when '5' => PentDobl_Test_Vector(nbr,numdeg,dendeg,dim);
        when '6' => OctoDobl_Test_Vector(nbr,numdeg,dendeg,dim);
        when '7' => DecaDobl_Test_Vector(nbr,numdeg,dendeg,dim);
        when others => null;
      end case;
    end if;
  end Main;

end Test_Rational_Approximations;
