with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Parse_Numbers;            
with Multprec_Floating_Numbers_io;
with Quad_Double_Numbers_io;

package body Multprec_QuadDobl_Convertors is

  function to_floating_number ( d : quad_double ) return Floating_Number is

    res : Floating_Number;
    s : string(1..80);
    ends : integer;
    p : integer := 1;

  begin
    Quad_Double_Numbers_io.to_string(d,64,0,false,false,true,' ',s,ends);
    Multprec_Parse_Numbers.Parse(s(1..ends),p,res);
    return res;
  end to_floating_number;

  function to_quad_double ( i : Integer_Number ) return quad_double is

    res : quad_double;
    f : Floating_Number := Create(i);

  begin
    res := to_quad_double(f);
    Clear(f);
    return res;
  end to_quad_double;

  function to_quad_double ( f : Floating_Number ) return quad_double is

    res : quad_double;
    sz : constant natural32 := Multprec_Floating_Numbers_io.Character_Size(f);
    s : string(1..natural(sz));
    fail : boolean;

  begin
    Multprec_Floating_Numbers_io.put(s,f);
    Quad_Double_Numbers_io.read(s,res,fail);
    return res;
  end to_quad_double;

end Multprec_QuadDobl_Convertors;
