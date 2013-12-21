package body C_Double_Arrays is

  function Concat ( a,b : C_Double_Array ) return C_Double_Array is

    res : C_Double_Array(0..Interfaces.C.size_T(a'length+b'length-1));
    ind : Interfaces.C.size_T := 0;

    use Interfaces.C;

  begin
    for i in a'range loop
      res(ind) := a(i);
      ind := ind+1;
    end loop;
    for i in b'range loop
      res(ind) := b(i);
      ind := ind+1;
    end loop;
    return res;
  end Concat;

end C_Double_Arrays;
