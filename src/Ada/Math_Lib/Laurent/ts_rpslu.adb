with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Double_Real_MatVecs;
with Double_Complex_MatVecs;
with Double_rpSeries_Operations;
with Test_Real_Powered_Series;

procedure ts_rpslu is

-- DESCRIPTION :
--   Solves a linear system of real powered series via row reduction.

  procedure Random_rpSeries_Vector
              ( nbt : in integer32;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Generates a random vector of real powered series,
  --   returned in (cff, pwr), with series of size nbt.

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

  -- DESCRIPTION :
  --   Generates a random matrix returned by coefficients
  --   and powers in (cff, pwr), with series of size nbt.

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

  -- DESCRIPTION :
  --   Writes the vector with coefficients and powers in (cff, pwr).

  begin
    for i in cff'range loop
      put("["); put(i,1); put_line("] :");
      Test_Real_Powered_Series.Write(cff(i).all,pwr(i).all);
    end loop;
  end Write;

  procedure Write ( cff : in Double_Complex_MatVecs.MatVec;
                    pwr : in Double_Real_MatVecs.MatVec ) is

  -- DESCRIPTION :
  --   Writes the matrix with coefficients and powers in (cff, pwr).

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

  -- DESCRIPTION :
  --   Collects all numbers in v, assuming all vectors have the same size.

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

  -- DESCRIPTION :
  --   Returns the index of x in v, or zero if x is not in v.

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

  -- DESCRIPTION :
  --   Multiplies the matrix in (Acff, Apwr) with the vector in (Xcff, Xpwr)
  --   and returns the result in (Bcff, Bpwr).

  -- ASSUMED : all series have the same size.

  -- REQUIRED : the dimensions must be compatible, that is:
  --    Bcff'range = Acff'range(1), Bpwr'range = Apwr'range(1)
  --    Xcff'range = Acff'range(2), Xpwr'range = Apwr'range(2).

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

  procedure Test ( dim,nbt : in integer32 ) is

  -- DESCRIPTION :
  --   Tests a system of dimension dim with number of terms equal to nbt.

    Acff : Double_Complex_MatVecs.MatVec(1..dim,1..dim);
    Apwr : Double_Real_MatVecs.MatVec(1..dim,1..dim);
    Xcff : Standard_Complex_VecVecs.VecVec(1..dim);
    Xpwr : Standard_Floating_VecVecs.VecVec(1..dim);
    Bcff : Standard_Complex_VecVecs.VecVec(1..dim);
    Bpwr : Standard_Floating_VecVecs.VecVec(1..dim);

  begin
    Random_rpSeries_Matrix(nbt,Acff,Apwr);
    put_line("A random real powered series matrix :"); Write(Acff,Apwr);
    Random_rpSeries_Vector(nbt,Xcff,Xpwr);
    put_line("A random real powered series vector :"); Write(Xcff,Xpwr);
    Matrix_Vector_Multiply(Acff,Xcff,Xpwr,Bcff,Bpwr);
    put_line("The right hand side vector :"); Write(Bcff,Bpwr);
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
