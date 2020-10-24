with text_io;                           use text_io;
with PentDobl_Complex_Matrices_io;      use PentDobl_Complex_Matrices_io;
with PentDobl_Complex_Linear_Solvers;   use PentDobl_Complex_Linear_Solvers;

package body PentDobl_Rational_Approximations is

  procedure Denominator_System
              ( numdeg,dendeg : in integer32; 
                cff : in PentDobl_Complex_Vectors.Vector;
                mat : out PentDobl_Complex_Matrices.Matrix;
                rhs : out PentDobl_Complex_Vectors.Vector ) is

    dim : constant integer32 := numdeg+dendeg;
    idx : integer32 := 0;
    zero : constant penta_double := create(0.0);

  begin
    for i in 1..dendeg loop
      idx := numdeg - dendeg + i;
      for j in 1..dendeg loop
        if idx < 1
         then mat(i,j) := Create(zero);
         else mat(i,j) := cff(idx);
        end if;
        idx := idx + 1;
      end loop;
    end loop;
    idx := 0;
    for k in numdeg+1..dim loop
      idx := idx + 1;
      rhs(idx) := -cff(k);
    end loop;
  end Denominator_System;

  function Numerator_Coefficients
              ( numdeg,dendeg : integer32;
                dencff,sercff : PentDobl_Complex_Vectors.Vector )
              return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(0..numdeg);
    mindeg : integer32;

  begin
    res(0) := sercff(0);
    if dendeg <= numdeg
     then mindeg := dendeg;
     else mindeg := numdeg;
    end if;
    for i in 1..numdeg loop
      res(i) := sercff(i);
      for j in 1..i loop
        exit when (j > mindeg);
        res(i) := res(i) + dencff(j)*sercff(i-j);
      end loop; 
    end loop;
    return res;
  end Numerator_Coefficients;

  procedure Assign_Numerator_Coefficients
              ( numdeg,dendeg : in integer32;
                dencff,sercff : in PentDobl_Complex_Vectors.Vector;
                cff : out PentDobl_Complex_Vectors.Vector ) is

    mindeg : integer32;

  begin
    cff(0) := sercff(0);
    if dendeg <= numdeg
     then mindeg := dendeg;
     else mindeg := numdeg;
    end if;
    for i in 1..numdeg loop
      cff(i) := sercff(i);
      for j in 1..i loop
        exit when (j > mindeg);
        cff(i) := cff(i) + dencff(j)*sercff(i-j);
      end loop; 
    end loop;
  end Assign_Numerator_Coefficients;

  procedure Pade ( numdeg,dendeg : in integer32;
                   cff : in PentDobl_Complex_Vectors.Vector;
                   numcff,dencff : out PentDobl_Complex_Vectors.Vector;
                   mat : in out PentDobl_Complex_Matrices.Matrix;
                   rhs : in out PentDobl_Complex_Vectors.Vector;
                   ipvt : in out Standard_Integer_Vectors.Vector;
                   info : out integer32; verbose : in boolean := false ) is

    zero : constant penta_double := create(0.0);
    cmplx_zero : constant Complex_Number := create(zero);
    one : constant penta_double := create(1.0);

  begin
    Denominator_System(numdeg,dendeg,cff,mat,rhs);
    if verbose then
      put_line("The matrix of the denominator system :");
      put(mat,3);
    end if;
    lufac(mat,dendeg,ipvt,info);
    if info = 0 then
      lusolve(mat,dendeg,ipvt,rhs);
      dencff(0) := Create(one);
      for i in 1..dendeg loop
        dencff(i) := rhs(dendeg-i+1);
      end loop;
      numcff := Numerator_Coefficients(numdeg,dendeg,dencff,cff);
    else -- singular coefficient matrix detected
      numcff := (0..numdeg => Create(zero));
      dencff := (0..dendeg => Create(zero));
      dencff(0) := Create(one); -- denominator is just one
      for k in rhs'range loop -- check if right hand side is zero
        if rhs(k) /= cmplx_zero
         then return;
        end if;
      end loop;
     -- numcff := Numerator_Coefficients(numdeg,dendeg,dencff,cff);
      Assign_Numerator_Coefficients(numdeg,dendeg,dencff,cff,numcff);
      info := 0; -- ignore the singular coefficient matrix
    end if;
  end Pade;

  procedure Pade ( numdeg,dendeg : in integer32;
                   cff : in PentDobl_Complex_Vectors.Vector;
                   numcff,dencff : out PentDobl_Complex_Vectors.Vector;
                   info : out integer32; verbose : in boolean := false ) is

    mat : PentDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : PentDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);

  begin
    Pade(numdeg,dendeg,cff,numcff,dencff,mat,rhs,ipvt,info,verbose);
  end Pade;

  procedure Pade_Vector
              ( numdeg,dendeg : in integer32;
                cff : in PentDobl_Complex_VecVecs.VecVec;
                numcff,dencff : in PentDobl_Complex_VecVecs.VecVec;
                mat : in out PentDobl_Complex_Matrices.Matrix;
                rhs : in out PentDobl_Complex_Vectors.Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := false ) is

    lnkcff,lnknum,lnkden : PentDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in cff'range loop
      lnkcff := cff(i); lnknum := numcff(i); lnkden := dencff(i);
      Pade(numdeg,dendeg,lnkcff.all,lnknum.all,lnkden.all,
           mat,rhs,ipvt,info,verbose);
    end loop;
  end Pade_Vector;

  function Evaluate
              ( p : PentDobl_Complex_Vectors.Vector;
                x : Complex_Number ) return Complex_Number is

    res : Complex_Number := p(p'last);

  begin
    for k in reverse 0..p'last-1 loop
      res := res*x + p(k);
    end loop;
    return res;
  end Evaluate;

  function Evaluate
              ( num,den : PentDobl_Complex_Vectors.Vector;
                x : Complex_Number ) return Complex_Number is

    numval : constant Complex_Number := Evaluate(num,x);
    denval : constant Complex_Number := Evaluate(den,x);
    res : constant Complex_Number := numval/denval;

  begin
    return res;
  end Evaluate;

  function Evaluate
              ( num,den : PentDobl_Complex_Vectors.Vector;
                x : penta_double ) return Complex_Number is

    cx : constant Complex_Number := create(x);

  begin
    return Evaluate(num,den,cx);
  end Evaluate;

  procedure Evaluate
              ( num,den : in PentDobl_Complex_VecVecs.VecVec;
                x : in Complex_Number;
                eva : out PentDobl_Complex_Vectors.Vector ) is
  begin
    for k in eva'range loop
      eva(k) := Evaluate(num(k).all,den(k).all,x);
    end loop;
  end Evaluate;

  procedure Evaluate
              ( num,den : in PentDobl_Complex_VecVecs.VecVec;
                x : in penta_double;
                eva : out PentDobl_Complex_Vectors.Vector ) is

    cx : constant Complex_Number := Create(x);

  begin
    Evaluate(num,den,cx,eva);
  end Evaluate;

end PentDobl_Rational_Approximations;
