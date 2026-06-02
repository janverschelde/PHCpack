with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Double_Real_Power_Series_IO;
with Test_Real_Powered_Series;
with Double_Linear_rpSeries_Solver;     use Double_Linear_rpSeries_Solver;

package body Test_Linear_Series_Solver is

  function Random_Series
             ( size : integer32 )
             return Double_Real_Power_Series.Link_to_Series is

    rep : Double_Real_Power_Series.Series(size);
    res : Double_Real_Power_Series.Link_to_Series;

  begin
    Test_Real_Powered_Series.Random_Series(size,rep.cff,rep.pwt);
    res := new Double_Real_Power_Series.Series'(rep);
    return res;
  end Random_Series;

  function Random_Series_Vector
             ( dim,size : integer32 )
             return Double_rpSeries_Vectors.Vector is

    res : Double_rpSeries_Vectors.Vector(1..dim);

  begin
    for i in res'range loop
      res(i) := Random_Series(size);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Series_Matrix
             ( dim,size : integer32 )
             return Double_rpSeries_Matrices.Matrix is

    res : Double_rpSeries_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random_Series(size);
      end loop;
    end loop;
    return res;
  end Random_Series_Matrix;

  procedure Test ( dim,nbr : in integer32 ) is

    x : constant Double_rpSeries_Vectors.Vector(1..dim)
      := Random_Series_Vector(dim,nbr);
    A : constant Double_rpSeries_Matrices.Matrix(1..dim,1..dim)
      := Random_Series_Matrix(dim,nbr);
    b : constant Double_rpSeries_Vectors.Vector(1..dim)
      := Right_Hand_Side(A,x);
    z0,z1c : Standard_Complex_Vectors.Vector(1..dim);
    z1p : Standard_Floating_Vectors.Vector(1..dim);
    sumerr : double_float := 0.0;

  begin
    put_line("a random vector x :"); Double_Real_Power_Series_IO.Write(x);
    put_line("a random matrix A :"); Double_Real_Power_Series_IO.Write(A);
    put_line("the right hand side vector b :");
    Double_Real_Power_Series_IO.Write(b);
    Real_Power_Series_Solver(A,b,z0,z1c,z1p,2);
    new_line;
    put_line("the solution series x :"); Double_Real_Power_Series_IO.Write(x);
    put_line("the computed series :");
    for i in z0'range loop
      put(z0(i)); new_line;
      sumerr := sumerr + AbsVal(z0(i) - x(i).cff(0));
      put(z1c(i)); put(" * t**"); put(z1p(i),1,14,3); new_line;
      sumerr := sumerr + AbsVal(z1c(i) - x(i).cff(1));
      sumerr := sumerr + abs(z1p(i) - x(i).pwt(1));
    end loop;
    put("sum of errors :"); put(sumerr,3); new_line;
  end Test;

  procedure Main is

    dim,nbr : integer32 := 0;
 
  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the number of terms : "); get(nbr);
    Test(dim,nbr);
  end Main;

end Test_Linear_Series_Solver;
