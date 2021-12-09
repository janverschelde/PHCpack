with HexaDobl_Random_Numbers;

package body HexaDobl_Random_Vectors is

  function Random_Vector ( first,last : integer32 )
                         return Hexa_Double_Vectors.Vector is

    res : Hexa_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := HexaDobl_Random_Numbers.Random;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Hexa_Double_Vectors.Vector ) is
  begin
    for i in v'range loop
      HexaDobl_Random_Numbers.Random_Hexa_Double(seed,v(i));
    end loop;
  end Random_Vector;

  function Random_Vector ( first,last : integer32 )
                         return HexaDobl_Complex_Vectors.Vector is

    res : HexaDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := HexaDobl_Random_Numbers.Random1;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out HexaDobl_Complex_Vectors.Vector ) is
  begin
    for i in v'range loop
      HexaDobl_Random_Numbers.Random1_Complex_Number(seed,v(i));
    end loop;
  end Random_Vector;

end HexaDobl_Random_Vectors;
