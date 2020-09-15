with OctoDobl_Random_Numbers;

package body OctoDobl_Random_Vectors is

  function Random_Vector ( first,last : integer32 )
                         return Octo_Double_Vectors.Vector is

    res : Octo_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := OctoDobl_Random_Numbers.Random;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Octo_Double_Vectors.Vector ) is
  begin
    for i in v'range loop
      OctoDobl_Random_Numbers.Random_Octo_Double(seed,v(i));
    end loop;
  end Random_Vector;

  function Random_Vector ( first,last : integer32 )
                         return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := OctoDobl_Random_Numbers.Random1;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out OctoDobl_Complex_Vectors.Vector ) is
  begin
    for i in v'range loop
      OctoDobl_Random_Numbers.Random1_Complex_Number(seed,v(i));
    end loop;
  end Random_Vector;

end OctoDobl_Random_Vectors;
