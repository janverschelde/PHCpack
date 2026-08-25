with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Linear_Solvers;
with Double_rpSeries_Operations;
with Test_Real_Powered_Series;

package body Test_rpSeries_LU_Solver is

  procedure Random_rpSeries_Vector
              ( nbt : in integer32;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec ) is
  begin
    for i in cff'range loop
      declare
        icff : Standard_Complex_Vectors.Vector(0..nbt);
        ipwr : Standard_Floating_Vectors.Vector(1..nbt);
      begin
        Test_Real_Powered_Series.Random_Series(nbt,icff,ipwr);
        cff(i) := new Standard_Complex_Vectors.Vector'(icff);
        pwr(i) := new Standard_Floating_Vectors.Vector'(ipwr);
      end;
    end loop;
  end Random_rpSeries_Vector;

  procedure Random_rpSeries_Matrix
              ( nbt : in integer32;
                cff : out Double_Complex_MatVecs.MatVec;
                pwr : out Double_Real_MatVecs.MatVec ) is
  begin
    for i in cff'range(1) loop
      for j in cff'range(2) loop
        declare
          ijcff : Standard_Complex_Vectors.Vector(0..nbt);
          ijpwr : Standard_Floating_Vectors.Vector(1..nbt);
        begin
          Test_Real_Powered_Series.Random_Series(nbt,ijcff,ijpwr);
          cff(i,j) := new Standard_Complex_Vectors.Vector'(ijcff);
          pwr(i,j) := new Standard_Floating_Vectors.Vector'(ijpwr);
        end;
      end loop;
    end loop;
  end Random_rpSeries_Matrix;

  procedure Write ( cff : in Standard_Complex_VecVecs.VecVec;
                    pwr : in Standard_Floating_VecVecs.VecVec ) is
  begin
    for i in cff'range loop
      put("["); put(i,1); put_line("] :");
      Test_Real_Powered_Series.Write(cff(i).all,pwr(i).all);
    end loop;
  end Write;

  procedure Write ( cff : in Double_Complex_MatVecs.MatVec;
                    pwr : in Double_Real_MatVecs.MatVec ) is
  begin
    for i in cff'range(1) loop
      for j in cff'range(2) loop
        put("["); put(i,1); put(",");
        put(j,1); put_line("] :");
        Test_Real_Powered_Series.Write(cff(i,j).all,pwr(i,j).all);
      end loop;
    end loop;
  end Write;

  function Flatten ( v : Standard_Floating_VecVecs.VecVec ) 
                   return Standard_Floating_Vectors.Vector is

    dim : constant integer32 := v'last*v(v'first)'last;
    res : Standard_Floating_Vectors.Vector(1..dim);
    idx : integer32 := 0;

  begin
    for i in v'range loop
      declare
        iv : constant Standard_Floating_Vectors.Link_to_Vector := v(i);
      begin
        for j in iv'range loop
          idx := idx + 1;
          res(idx) := iv(j);
        end loop;
      end;
    end loop;
    return res;
  end Flatten;

  function Index ( v : Standard_Floating_Vectors.Vector;
                   x : double_float ) return integer32 is
  begin
    for i in v'range loop
      if v(i) = x
       then return i;
      end if;
    end loop;
    return 0;
  end Index;

  procedure Matrix_Vector_Multiply
              ( Acff : in Double_Complex_MatVecs.MatVec;
                Xcff : in Standard_Complex_VecVecs.VecVec;
                Xpwr : in Standard_Floating_VecVecs.VecVec;
                Bcff : out Standard_Complex_VecVecs.VecVec;
                Bpwr : out Standard_Floating_VecVecs.VecVec ) is

    dim : constant integer32 := Xcff'last;
    nbt : constant integer32 := Xcff(Xcff'first)'last;
    Bsize : constant integer32 := dim*nbt;
    pwrs : Standard_Floating_Vectors.Vector(1..Bsize) := Flatten(Xpwr);
    cnst : Complex_Number;
    iBcff : Standard_Complex_Vectors.Link_to_Vector;
    idx : integer32;

  begin
    Double_rpSeries_Operations.Sort(pwrs);
    for i in Acff'range(1) loop
      Bcff(i) := new Standard_Complex_Vectors.Vector'(0..Bsize => create(0.0));
      Bpwr(i) := new Standard_Floating_Vectors.Vector'(pwrs);
      cnst := Standard_Complex_Numbers.create(0.0);
      iBcff := Bcff(i);
      for j in Acff'range(2) loop
        declare
          acf : constant Standard_Complex_Vectors.Link_to_Vector := Acff(i,j);
          bcf : constant Standard_Complex_Vectors.Link_to_Vector := Xcff(j);
          bpw : constant Standard_Floating_Vectors.Link_to_Vector := Xpwr(j);
        begin
          cnst := cnst + acf(0)*bcf(0);
          for k in bpw'range loop
            idx := Index(pwrs,bpw(k));
            put("index : "); put(idx,1); new_line;
            if idx > 0
             then iBcff(idx) := iBcff(idx) + acf(0)*bcf(k);
             else put("Zero index : "); put(bpw(k)); put_line(" not found!?");
            end if;
          end loop;
        end;
      end loop;
      iBcff(0) := cnst;
    end loop;
  end Matrix_Vector_Multiply;

  function Constant_Coefficients
             ( A : Double_Complex_MatVecs.MatVec ) 
             return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));
    Aij : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        Aij := A(i,j);
        res(i,j) := Aij(0);
      end loop;
    end loop;
    return res;
  end Constant_Coefficients;

  procedure Solve_Linear_System
              ( Acf0 : in Standard_Complex_Matrices.Matrix;
                Bcff : in Standard_Complex_VecVecs.VecVec;
                Bpwr : in Standard_Floating_VecVecs.VecVec;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec ) is

    tolzero : constant double_float := 1.0e-12;
    dim : constant integer32 := Acf0'last;
    wrk : Standard_Complex_Matrices.Matrix(Acf0'range(1),Acf0'range(2))
        := Acf0;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    rhs : Standard_Complex_Vectors.Vector(1..dim);
    Bsize : constant integer32 := Bcff(Bcff'first)'last;
    nbt : constant integer32 := Bsize/dim;
    Xpwr : double_float;
    next : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);

  begin
    for i in cff'range loop
      cff(i) := new Standard_Complex_Vectors.Vector'(0..nbt => create(0.0));
      pwr(i) := new Standard_Floating_Vectors.Vector'(1..nbt => 0.0);
    end loop;
    Standard_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
    for i in 1..dim loop
      rhs(i) := Bcff(i)(0);
    end loop;
    Standard_Complex_Linear_Solvers.lusolve(wrk,dim,ipvt,rhs);
    put_line("the constant coefficients of the solution :");
    put_line(rhs);
    for i in 1..dim loop
      cff(i)(0) := rhs(i);
    end loop;
    for k in 1..Bsize loop
      Xpwr := Bpwr(Bpwr'first)(k);
      for i in 1..dim loop
        rhs(i) := Bcff(i)(k);
      end loop;
      Standard_Complex_Linear_Solvers.lusolve(wrk,dim,ipvt,rhs);
      put("the coefficients "); put(k,1); put_line(" of the solution :");
      put_line(rhs);
      for i in 1..dim loop
        if absVal(rhs(i)) > tolzero then
          next(i) := next(i) + 1;
          cff(i)(next(i)) := rhs(i);
          pwr(i)(next(i)) := Xpwr;
        end if;
      end loop;
    end loop;
  end Solve_Linear_System;

  function Difference_Sum
             ( x,y : Standard_Complex_VecVecs.VecVec ) return double_float is

    res : double_float := 0.0;
    ix,iy : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      ix := x(i); iy := y(i);
      for k in ix'range loop
        res := res + AbsVal(ix(k) - iy(k));
      end loop;
    end loop;
    return res;
  end Difference_Sum;

  function Difference_Sum
             ( x,y : Standard_Floating_VecVecs.VecVec ) return double_float is

    res : double_float := 0.0;
    ix,iy : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      ix := x(i); iy := y(i);
      for k in ix'range loop
        res := res + abs(ix(k) - iy(k));
      end loop;
    end loop;
    return res;
  end Difference_Sum;

  procedure Test ( dim,nbt : in integer32 ) is

    Acff : Double_Complex_MatVecs.MatVec(1..dim,1..dim);
    Apwr : Double_Real_MatVecs.MatVec(1..dim,1..dim);
    Xcff : Standard_Complex_VecVecs.VecVec(1..dim);
    Xpwr : Standard_Floating_VecVecs.VecVec(1..dim);
    Bcff : Standard_Complex_VecVecs.VecVec(1..dim);
    Bpwr : Standard_Floating_VecVecs.VecVec(1..dim);
    Ycff : Standard_Complex_VecVecs.VecVec(1..dim);
    Ypwr : Standard_Floating_VecVecs.VecVec(1..dim);
    Acf0 : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    errcff,errpwr : double_float;

  begin
    Random_rpSeries_Matrix(nbt,Acff,Apwr);
    put_line("A random real powered series matrix :"); Write(Acff,Apwr);
    Random_rpSeries_Vector(nbt,Xcff,Xpwr);
    put_line("A random real powered series vector :"); Write(Xcff,Xpwr);
    Matrix_Vector_Multiply(Acff,Xcff,Xpwr,Bcff,Bpwr);
    put_line("The right hand side vector :"); Write(Bcff,Bpwr);
    Acf0 := Constant_Coefficients(Acff);
    Solve_Linear_System(Acf0,Bcff,Bpwr,Ycff,Ypwr);
    put_line("The computed solution :"); Write(Ycff,Ypwr);
    put_line("The test solution :"); Write(Xcff,Xpwr);
    errcff := Difference_Sum(Ycff,Xcff);
    put("Error on coefficients :"); put(errcff,2); new_line;
    errpwr := Difference_Sum(Ypwr,Xpwr);
    put("Error on powers :"); put(errpwr,2); new_line;
  end Test;

  procedure Main is

    dim,nbt : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    put("Give the number of terms : "); get(nbt);
    Test(dim,nbt);
  end Main;

end Test_rpSeries_LU_Solver;
