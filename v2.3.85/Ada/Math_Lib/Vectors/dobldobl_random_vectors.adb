with DoblDobl_Random_Numbers;

package body DoblDobl_Random_Vectors is

  function Random_Vector ( first,last : integer32 )
                         return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Random_Numbers.Random;
    end loop;
    return res;
  end Random_Vector;

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Random_Numbers.Random_Magnitude(m);
    end loop;
    return res;
  end Random_Vector;

  function Random_Vector ( first,last : integer32 )
                         return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Random_Numbers.Random1;
    end loop;
    return res;
  end Random_Vector;

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Random_Numbers.Random_Magnitude(m);
    end loop;
    return res;
  end Random_Vector;

end DoblDobl_Random_Vectors;
