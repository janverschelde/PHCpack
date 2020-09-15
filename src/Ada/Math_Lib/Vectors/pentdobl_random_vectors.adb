with PentDobl_Random_Numbers;

package body PentDobl_Random_Vectors is

  function Random_Vector ( first,last : integer32 )
                         return Penta_Double_Vectors.Vector is

    res : Penta_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := PentDobl_Random_Numbers.Random;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Penta_Double_Vectors.Vector ) is
  begin
    for i in v'range loop
      PentDobl_Random_Numbers.Random_Penta_Double(seed,v(i));
    end loop;
  end Random_Vector;

  function Random_Vector ( first,last : integer32 )
                         return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := PentDobl_Random_Numbers.Random1;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out PentDobl_Complex_Vectors.Vector ) is
  begin
    for i in v'range loop
      PentDobl_Random_Numbers.Random1_Complex_Number(seed,v(i));
    end loop;
  end Random_Vector;

end PentDobl_Random_Vectors;
