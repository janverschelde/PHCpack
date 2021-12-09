with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Parse_Numbers;            
with Multprec_Floating_Numbers_io;
with Hexa_Double_Numbers_io;

-- for testing:
--with text_io; use text_io;

package body Multprec_HexaDobl_Convertors is

  function to_floating_number ( d : hexa_double ) return Floating_Number is

    res : Floating_Number;
    s : string(1..264);
    ends : integer;
    p : integer := 1;

  begin
    Hexa_Double_Numbers_io.to_string(d,256,0,false,false,true,' ',s,ends);
    Multprec_Parse_Numbers.Parse(s(1..ends),p,res);
    return res;
  end to_floating_number;

  function to_hexa_double ( i : Integer_Number ) return hexa_double is

    res : hexa_double;
    f : Floating_Number := Create(i);

  begin
    res := to_hexa_double(f);
    Clear(f);
    return res;
  end to_hexa_double; 

  function to_hexa_double ( f : Floating_Number ) return hexa_double is

    res : hexa_double;
    sz : constant natural32 := Multprec_Floating_Numbers_io.Character_Size(f);
    s : string(1..integer(sz));
    fail : boolean;

  begin
   -- put("converting ");
   -- Multprec_Floating_Numbers_io.put(f); new_line;
    Multprec_Floating_Numbers_io.put(s,f);
    Hexa_Double_Numbers_io.read(s,res,fail);
   -- put("res : ");
   -- Hexa_Double_Numbers_io.put(res);
   -- new_line;
    return res;
  end to_hexa_double;

end Multprec_HexaDobl_Convertors;
