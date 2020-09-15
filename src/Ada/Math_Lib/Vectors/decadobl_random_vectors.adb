with DecaDobl_Random_Numbers;

package body DecaDobl_Random_Vectors is

  function Random_Vector ( first,last : integer32 )
                         return Deca_Double_Vectors.Vector is

    res : Deca_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := DecaDobl_Random_Numbers.Random;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Deca_Double_Vectors.Vector ) is
  begin
    for i in v'range loop
      DecaDobl_Random_Numbers.Random_Deca_Double(seed,v(i));
    end loop;
  end Random_Vector;

  function Random_Vector ( first,last : integer32 )
                         return DecaDobl_Complex_Vectors.Vector is

    res : DecaDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := DecaDobl_Random_Numbers.Random1;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out DecaDobl_Complex_Vectors.Vector ) is
  begin
    for i in v'range loop
      DecaDobl_Random_Numbers.Random1_Complex_Number(seed,v(i));
    end loop;
  end Random_Vector;

end DecaDobl_Random_Vectors;
