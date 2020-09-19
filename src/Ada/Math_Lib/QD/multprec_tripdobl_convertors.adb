with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Parse_Numbers;            
with Multprec_Floating_Numbers_io;
with Triple_Double_Numbers_io;

-- for testing:
--with text_io; use text_io;

package body Multprec_TripDobl_Convertors is

  function to_floating_number ( d : triple_double ) return Floating_Number is

    res : Floating_Number;
    s : string(1..58);
    ends : integer;
    p : integer := 1;

  begin
    Triple_Double_Numbers_io.to_string(d,50,0,false,false,true,' ',s,ends);
    Multprec_Parse_Numbers.Parse(s(1..ends),p,res);
    return res;
  end to_floating_number;

  function to_triple_double ( i : Integer_Number ) return triple_double is

    res : triple_double;
    f : Floating_Number := Create(i);

  begin
    res := to_triple_double(f);
    Clear(f);
    return res;
  end to_triple_double; 

  function to_triple_double ( f : Floating_Number ) return triple_double is

    res : triple_double;
    sz : constant natural32 := Multprec_Floating_Numbers_io.Character_Size(f);
    s : string(1..integer(sz));
    fail : boolean;

  begin
   -- put("converting ");
   -- Multprec_Floating_Numbers_io.put(f); new_line;
    Multprec_Floating_Numbers_io.put(s,f);
    Triple_Double_Numbers_io.read(s,res,fail);
   -- put("res : ");
   -- Triple_Double_Numbers_io.put(res);
   -- new_line;
    return res;
  end to_triple_double;

end Multprec_TripDobl_Convertors;
