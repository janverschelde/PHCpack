with TripDobl_Random_Numbers;

package body TripDobl_Random_Vectors is

  function Random_Vector ( first,last : integer32 )
                         return Triple_Double_Vectors.Vector is

    res : Triple_Double_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := TripDobl_Random_Numbers.Random;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Triple_Double_Vectors.Vector ) is
  begin
    for i in v'range loop
      TripDobl_Random_Numbers.Random_Triple_Double(seed,v(i));
    end loop;
  end Random_Vector;

  function Random_Vector ( first,last : integer32 )
                         return TripDobl_Complex_Vectors.Vector is

    res : TripDobl_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := TripDobl_Random_Numbers.Random1;
    end loop;
    return res;
  end Random_Vector;

  procedure Random_Vector
              ( seed : in out integer32;
                v : out TripDobl_Complex_Vectors.Vector ) is
  begin
    for i in v'range loop
      TripDobl_Random_Numbers.Random1_Complex_Number(seed,v(i));
    end loop;
  end Random_Vector;

end TripDobl_Random_Vectors;
