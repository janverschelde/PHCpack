with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Linear_Solvers;   use DoblDobl_Complex_Linear_Solvers;

package body DoblDobl_Rational_Approximations is

  procedure Denominator_System
              ( numdeg,dendeg : in integer32; 
                cff : in DoblDobl_Complex_Vectors.Vector;
                mat : out DoblDobl_Complex_Matrices.Matrix;
                rhs : out DoblDobl_Complex_Vectors.Vector ) is

    dim : constant integer32 := numdeg+dendeg;
    idx : integer32 := 0;
    zero : constant double_double := create(0.0);

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
                dencff,sercff : DoblDobl_Complex_Vectors.Vector )
              return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(0..numdeg);
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

  procedure Pade
              ( numdeg,dendeg : in integer32;
                cff : in DoblDobl_Complex_Vectors.Vector;
                numcff,dencff : out DoblDobl_Complex_Vectors.Vector ) is

    dim : constant integer32 := numdeg + dendeg;
    mat : DoblDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : DoblDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);

  begin
    Denominator_System(numdeg,dendeg,cff,mat,rhs);
    lufac(mat,dendeg,ipvt,info);
    if info = 0 then
      lusolve(mat,dendeg,ipvt,rhs);
      dencff(0) := Create(one);
      for i in 1..dendeg loop
        dencff(i) := rhs(dendeg-i+1);
      end loop;
      numcff := Numerator_Coefficients(numdeg,dendeg,dencff,cff);
    else
      numcff := (0..numdeg => Create(zero));
      dencff := (0..dendeg => Create(zero));
    end if;
  end Pade;

  function Evaluate
              ( p : DoblDobl_Complex_Vectors.Vector;
                x : Complex_Number ) return Complex_Number is

    res : Complex_Number := p(p'last);

  begin
    for k in reverse 0..p'last-1 loop
      res := res*x + p(k);
    end loop;
    return res;
  end Evaluate;

  function Evaluate
              ( num,den : DoblDobl_Complex_Vectors.Vector;
                x : Complex_Number ) return Complex_Number is

    numval : constant Complex_Number := Evaluate(num,x);
    denval : constant Complex_Number := Evaluate(den,x);
    res : constant Complex_Number := numval/denval;

  begin
    return res;
  end Evaluate;

end DoblDobl_Rational_Approximations;
