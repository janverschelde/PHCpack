with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Test_Leading_Terms;

package body Test_Linear_Series_Solver is

  procedure Test ( dim,nbr : in integer32 ) is

    X : Standard_Floating_Matrices.Matrix(1..dim,1..nbr);
    cX : Standard_Complex_Matrices.Matrix(1..dim,1..nbr);

  begin
    Test_Leading_Terms.random_series(dim,nbr,X,cX);
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
