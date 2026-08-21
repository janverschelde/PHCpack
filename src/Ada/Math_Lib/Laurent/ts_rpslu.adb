with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Double_Real_MatVecs;
with Double_Complex_MatVecs;
with Test_Real_Powered_Series;

procedure ts_rpslu is

-- DESCRIPTION :
--   Solves a linear system of real powered series via row reduction.

  procedure Random_rpSeries_Matrix
              ( nbt : in integer32;
                Acff : out Double_Complex_MatVecs.MatVec;
                Apwr : out Double_Real_MatVecs.MatVec ) is

  -- DESCRIPTION :
  --   Generates a random matrix represented by coefficients
  --   and powers in (Acff, Apwr), with series of size nbt.

  begin
    for i in Acff'range(1) loop
      for j in Acff'range(2) loop
        declare
          cff : Standard_Complex_Vectors.Vector(0..nbt);
          pwr : Standard_Floating_Vectors.Vector(1..nbt);
        begin
          Test_Real_Powered_Series.Random_Series(nbt,cff,pwr);
          Acff(i,j) := new Standard_Complex_Vectors.Vector'(cff);
          Apwr(i,j) := new Standard_Floating_Vectors.Vector'(pwr);
        end;
      end loop;
    end loop;
  end Random_rpSeries_Matrix;

  procedure Write ( Acff : in Double_Complex_MatVecs.MatVec;
                    Apwr : in Double_Real_MatVecs.MatVec ) is

  -- DESCRIPTION :
  --   Writes the matrix with coefficients and powers in (Acff, Apwr).

  begin
    for i in Acff'range(1) loop
      for j in Acff'range(2) loop
        put("["); put(i,1); put(",");
        put(j,1); put_line("] :");
        Test_Real_Powered_Series.Write(Acff(i,j).all,Apwr(i,j).all);
      end loop;
    end loop;
  end Write;

  procedure Test ( dim,nbt : in integer32 ) is

  -- DESCRIPTION :
  --   Tests a system of dimension dim with number of terms equal to nbt.

    Acff : Double_Complex_MatVecs.MatVec(1..dim,1..dim);
    Apwr : Double_Real_MatVecs.MatVec(1..dim,1..dim);

  begin
    Random_rpSeries_Matrix(nbt,Acff,Apwr);
    put_line("A random real powered series matrix :");
    Write(Acff,Apwr);
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the dimensions and then generates a test problem.

    dim,nbt : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    put("Give the number of terms : "); get(nbt);
    Test(dim,nbt);
  end Main;

begin
  Main;
end ts_rpslu;
