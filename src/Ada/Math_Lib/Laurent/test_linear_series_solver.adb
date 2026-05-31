with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Double_Puiseux_Operations;
with Test_Real_Powered_Series;

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

  procedure Write ( x : in Double_rpSeries_Vectors.Vector ) is
  begin
    for i in x'range loop
      Test_Real_Powered_Series.Write(x(i).cff,x(i).pwt);
    end loop;
  end Write;

  procedure Write ( A : in Double_rpSeries_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        Test_Real_Powered_Series.Write(A(i,j).cff,A(i,j).pwt);
      end loop;
    end loop;
  end Write;

  function Right_Hand_Side
             ( A : Double_rpSeries_Matrices.Matrix;
               x : in Double_rpSeries_Vectors.Vector ) 
             return Double_rpSeries_Vectors.Vector is

    res : Double_rpSeries_Vectors.Vector(A'range);

    use Double_Real_Power_Series;

  begin
    for i in A'range(1) loop
      res(i) := A(i,A'first(2))*x(x'first);
      for j in A'first(2)+1..A'last(2) loop
        declare
          prd : Link_to_Series := A(i,j)*x(j);
          sum : Link_to_Series := res(i) + prd;
        begin
          Copy(sum,res(i)); Clear(prd); clear(sum);
        end;
      end loop;
    end loop;
    return res;
  end Right_Hand_Side;

  function Extract_Constants
             ( A : Double_rpSeries_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j).cff(0);
      end loop;
    end loop;
    return res;
  end Extract_Constants;

  function Extract_Constants 
             ( v : Double_rpSeries_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := v(i).cff(0);
    end loop;
    return res;
  end Extract_Constants;

  procedure Test_Series_Solver 
              ( A : in Double_rpSeries_Matrices.Matrix;
                x,b : in Double_rpSeries_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    A0 : constant Standard_Complex_Matrices.Matrix(A'range(1),A'range(2))
       := Extract_Constants(A);
    b0 : constant Standard_Complex_Vectors.Vector(b'range)
       := Extract_Constants(b);
    x0 : constant Standard_Complex_Vectors.Vector(x'range)
       := Extract_Constants(x);
    z0 : Standard_Complex_Vectors.Vector(A'range(2));
    rcond : double_float;

  begin
    Double_Puiseux_Operations.Solve_Constant_Linear_System
      (A0,b0,z0,rcond,vrblvl-1);
    if vrblvl > 0 then
      put_line("the constant coefficients of the solution :");
      put_line(x0);
      put_line("the computed constant coefficients :");
      put_line(z0);
      put("rcond :"); put(rcond,3); new_line;
    end if;
  end Test_Series_Solver;

  procedure Test ( dim,nbr : in integer32 ) is

    x : constant Double_rpSeries_Vectors.Vector(1..dim)
      := Random_Series_Vector(dim,nbr);
    A : constant Double_rpSeries_Matrices.Matrix(1..dim,1..dim)
      := Random_Series_Matrix(dim,nbr);
    b : constant Double_rpSeries_Vectors.Vector(1..dim)
      := Right_Hand_Side(A,x);

  begin
    put_line("a random vector x :"); Write(x);
    put_line("a random matrix A :"); Write(A);
    put_line("the right hand side vector b :"); Write(b);
    Test_Series_Solver(A,x,b,2);
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
