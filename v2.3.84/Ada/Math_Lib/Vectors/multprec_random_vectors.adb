with Multprec_Random_Numbers;            use Multprec_Random_Numbers;

package body Multprec_Random_Vectors is

  function Random_Vector ( first,last : integer32; sz : natural32 )
                         return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := Random(sz);
    end loop;
    return res;
  end Random_Vector;

  function Random_Vector ( first,last : integer32; sz : natural32 )
                         return Multprec_Floating_Vectors.Vector is

    res : Multprec_Floating_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := Random(sz);
    end loop;
    return res;
  end Random_Vector;

  function Random_Vector ( first,last : integer32; sz : natural32 )
                         return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(first..last);

  begin
    for i in res'range loop
      res(i) := Random(sz);
    end loop;
    return res;
  end Random_Vector;

end Multprec_Random_Vectors;
