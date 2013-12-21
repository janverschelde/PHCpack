with Multprec_Parse_Numbers;

package body Multprec_Floating_Constants is

  function Pi ( d : natural32 ) return Floating_Number is

    res : Floating_Number;
    p : natural := 1;
    s : constant string := strpi(1..integer(d)+1) & " ";

  begin
    Multprec_Parse_Numbers.Parse(s,p,res);
    return res;
  end Pi;

  function TwoPi ( d : natural32 ) return Floating_Number is

    res : Floating_Number;
    p : natural := 1;
    s : constant string := strtwopi(1..integer(d)+1) & " ";

  begin
    Multprec_Parse_Numbers.Parse(s,p,res);
    return res;
  end TwoPi;

end Multprec_Floating_Constants;
