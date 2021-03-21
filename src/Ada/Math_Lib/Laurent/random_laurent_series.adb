with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Random_Matrices;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;

package body Random_Laurent_Series is

  procedure Random_Series_Coefficients
              ( dim,deg : in integer32;
                cff : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    res : Standard_Complex_VecVecs.VecVec(1..dim);

  begin
    for i in 1..dim loop
      declare
        val : constant Standard_Complex_Vectors.Vector(0..deg)
            := Standard_Random_Vectors.Random_Vector(0,deg);
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(val);
      end;
    end loop;
    cff := new Standard_Complex_VecVecs.VecVec'(res);
  end Random_Series_Coefficients;

  procedure Random_Vector
              ( dim,deg,low,upp : in integer32;
                e : out Standard_Integer_Vectors.Vector;
                c : out Standard_Complex_VecVecs.Link_to_VecVec ) is
  begin
    e := Standard_Random_Vectors.Random_Vector(1,dim,low,upp);
    Random_Series_Coefficients(dim,deg,c);
  end Random_Vector;

  procedure Random_VecVecVec
              ( v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is
  begin
    for i in v'range loop
      declare
        vi : constant Standard_Complex_VecVecs.Link_to_VecVec := v(i);
      begin
        for j in vi'range loop
          declare
            vij : constant Standard_Complex_Vectors.Link_to_Vector := vi(j);
          begin
            for k in vij'range loop
              vij(k) := Standard_Random_Numbers.Random1;
            end loop;
          end;
        end loop;
      end;
    end loop;
  end Random_VecVecVec;

  procedure Random_Lower_VecVecVec
              ( v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is
  begin
    for i in v'range loop -- row i
      declare
        vi : constant Standard_Complex_VecVecs.Link_to_VecVec := v(i);
        vii : constant Standard_Complex_Vectors.Link_to_Vector := vi(i);
      begin
        for j in vi'first..(i-1) loop -- column j
          declare
            vij : constant Standard_Complex_Vectors.Link_to_Vector := vi(j);
          begin
            for k in vij'range loop
              vij(k) := Standard_Random_Numbers.Random1;
            end loop;
          end;
        end loop;
        vii(0) := Standard_Complex_Numbers.Create(1.0);
      end;
    end loop;
  end Random_Lower_VecVecVec;

  procedure Random_Upper_VecVecVec
              ( v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is
  begin
    for i in v'range loop -- row i
      declare
        vi : constant Standard_Complex_VecVecs.Link_to_VecVec := v(i);
      begin
        for j in i..vi'last loop -- column j
          declare
            vij : constant Standard_Complex_Vectors.Link_to_Vector := vi(j);
          begin
            for k in vij'range loop
              vij(k) := Standard_Random_Numbers.Random1;
            end loop;
          end;
        end loop;
      end;
    end loop;
  end Random_Upper_VecVecVec;

  procedure Random_Matrix
              ( nbrows,nbcols : in natural32; low,upp : in integer32;
                e : out Standard_Integer_Matrices.Matrix;
                v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is
  begin
    e := Standard_Random_Matrices.Random_Matrix(nbrows,nbcols,low,upp);
    Random_VecVecVec(v);
  end Random_Matrix;

  procedure Random_Lower_Matrix
              ( nbrows,nbcols : in natural32; low,upp : in integer32;
                e : out Standard_Integer_Matrices.Matrix;
                v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is
  begin
    e := Standard_Random_Matrices.Random_Matrix(nbrows,nbcols,low,upp);
    for i in e'range(1) loop
      for j in i..e'last(2) loop
        e(i,j) := 0;
      end loop;
    end loop;
    Random_Lower_VecVecVec(v);
  end Random_Lower_Matrix;

  procedure Random_Upper_Matrix
              ( nbrows,nbcols : in natural32; low,upp : in integer32;
                e : out Standard_Integer_Matrices.Matrix;
                v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is
  begin
    e := Standard_Random_Matrices.Random_Matrix(nbrows,nbcols,low,upp);
    for i in e'range(1) loop
      for j in e'first(2)..i loop
        e(i,j) := 0;
      end loop;
    end loop;
    Random_Upper_VecVecVec(v);
  end Random_Upper_Matrix;

end Random_Laurent_Series;
