package body Lexicographic_Root_Enumeration is

  procedure Lexicographic_Enumeration
             ( k,n : in natural32;
               d : in Standard_Natural_Vectors.Vector;
               acc : in out Standard_Natural_Vectors.Vector;
               continue : out boolean ) is
  begin
    if k <= n then
      for i in 1..d(integer32(k)) loop
         acc(integer32(k)) := i;
         Lexicographic_Enumeration(k+1,n,d,acc,continue);
         exit when not continue;
       end loop;
    else
      Process(acc,continue);
    end if;
  end Lexicographic_Enumeration;

  function Consecutive_Products
              ( d : Standard_Natural_Vectors.Vector )
              return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(d'first..d'last-1);

  begin
    res(res'last) := d(d'last);
    for i in reverse res'first+1..res'last loop
      res(i-1) := res(i)*d(i);
    end loop;
    return res;
  end Consecutive_Products;

  function Root_Map ( n,k : natural32;
                      d : Standard_Natural_Vectors.Vector )
                    return Standard_Natural_Vectors.Vector is
  begin
    return Root_Map(n,k,d,Consecutive_Products(d));
  end Root_Map;

  function Root_Map ( n,k : natural32;
                      d,cp : Standard_Natural_Vectors.Vector )
                    return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(1..integer32(n));
    acc : natural32 := k;

  begin
    for i in cp'range loop
      res(i) := acc/cp(i) + 1;
      acc := acc mod cp(i);
    end loop;
    if acc > 0 then
      res(res'last) := acc;
    else
      res(res'last) := d(d'last);
      for i in reverse res'first..res'last-1 loop
        res(i) := res(i) - 1;
        exit when (res(i) > 0);
        res(i) := d(i);
      end loop;
    end if;
    return res;
  end Root_Map;

  procedure Add_One ( d : in Standard_Natural_Vectors.Vector;
                      x : in out Standard_Natural_Vectors.Vector;
                      ind : out integer32; fail : out boolean ) is
  begin
    ind := 0;
    fail := true;
    for i in reverse x'range loop
      if x(i) < d(i) then
        x(i) := x(i) + 1;
        for j in i+1..x'last loop
          x(j) := 1;
        end loop;
        ind := i;
        fail := false;
      end if;
      exit when not fail;
    end loop;
  end Add_One;

end Lexicographic_Root_Enumeration;
