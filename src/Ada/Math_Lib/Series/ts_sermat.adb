with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Dense_Series_io;            use Standard_Dense_Series_io;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Matrices;
with Standard_Random_Series;

procedure ts_sermat is

-- DESCRIPTION :
--   Test on matrices of truncated dense power series.

  procedure Main is

    A : Standard_Dense_Series_Matrices.Matrix(0..2,0..2)
      := Standard_Random_Series.Random_Series_Matrix(0,2,0,2,4);

  begin
    put_line("A 2-by-2 matrix of random series of order 4 :");
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        put(A(i,j));
      end loop;
    end loop;
  end Main;

begin
  Main;
end ts_sermat;
