with Standard_Dense_Series2_Vectors;
with Random_Series_Generators;
with Random_Series_Vectors;

package body Random_Series_Matrices is

  function Random_Series_VecVec
             ( vvfirst,vvlast,first,last,degree : integer32 ) return VecVec is

    res : VecVec(vvfirst..vvlast);

  begin
    for i in res'range loop
      declare
        v : constant Standard_Dense_Series2_Vectors.Vector(first..last)
          := Random_Series_Vectors.Random_Series_Vector(first,last,degree);
      begin
        res(i) := new Standard_Dense_Series2_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Random_Series_VecVec;

  function Random_Series_Matrix
             ( rowfirst,rowlast,columnfirst,columnlast,degree : integer32 )
             return Matrix is

    res : Matrix(rowfirst..rowlast,columnfirst..columnlast);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random_Series_Generators.Random_Series(degree);
      end loop;
    end loop;
    return res;
  end Random_Series_Matrix;

end Random_Series_Matrices;
