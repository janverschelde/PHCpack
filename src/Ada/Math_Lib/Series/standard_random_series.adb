with Standard_Complex_Vectors;
with Standard_Random_Vectors;

package body Standard_Random_Series is

  function Random_Series ( order : integer32 ) return Series is

    cff : Standard_Complex_Vectors.Vector(0..order)
        := Standard_Random_Vectors.Random_Vector(0,order);

  begin
    return Create(cff);
  end Random_Series;

  function Random_Series_Vector
             ( first,last,order : integer32 )
             return Standard_Dense_Series_Vectors.Vector is

    res : Standard_Dense_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := Random_Series(order);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Series_Matrix
             ( rowfirst,rowlast,columnfirst,columnlast,order : integer32 )
             return Standard_Dense_Series_Matrices.Matrix is

    use Standard_Dense_Series_Matrices;

    res : Matrix(rowfirst..rowlast,columnfirst..columnlast);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random_Series(order);
      end loop;
    end loop;
    return res;
  end Random_Series_Matrix;

end Standard_Random_Series;
