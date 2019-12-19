with Standard_Integer_Numbers;            use Standard_Integer_Numbers;

package body Exponent_Indices is

  function Exponent_Index
             ( xp : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    deg : constant integer32 := Standard_Integer_Vectors.Sum(xp);
    res : Standard_Integer_Vectors.Vector(1..deg); 
    idx : integer32 := 0;

  begin
    for k in xp'range loop
      if xp(k) = 1 then
        idx := idx + 1;
        res(idx) := k;
      end if;
    end loop;
    return res;
  end Exponent_Index;

  function Exponent_Index
             ( xp : Standard_Integer_VecVecs.VecVec )
             return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(xp'range);

  begin
    for i in xp'range loop
      declare
        idx : constant Standard_Integer_Vectors.Vector
            := Exponent_Index(xp(i).all);
      begin
        res(i) := new Standard_Integer_Vectors.Vector'(idx);
      end;
    end loop;
    return res;
  end Exponent_Index;

end Exponent_Indices;
