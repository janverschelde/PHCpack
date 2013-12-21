with QuadDobl_Random_Numbers;   

package body QuadDobl_Random_Vectors is

  function Random_Vector ( first,last : integer32 )
                         return Quad_Double_Vectors.Vector is

    res : Quad_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Random_Numbers.Random;
    end loop;
    return res;
  end Random_Vector;

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return Quad_Double_Vectors.Vector is

    res : Quad_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Random_Numbers.Random_Magnitude(m);
    end loop;
    return res;
  end Random_Vector;

  function Random_Vector ( first,last : integer32 )
                         return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Random_Numbers.Random1;
    end loop;
    return res;
  end Random_Vector;

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Random_Numbers.Random_Magnitude(m);
    end loop;
    return res;
  end Random_Vector;

end QuadDobl_Random_Vectors;
