with text_io,integer_io;               use text_io,integer_io;
with Interfaces.C.Strings;             use Interfaces.C.Strings;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Symbol_Table;
with Standard_Complex_Polynomials;     use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;  use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Strings;    use Standard_Complex_Poly_Strings;

procedure print_poly ( n : in integer; s : in chars_ptr ) is

-- IMPORTANT NOTICE :
--   For one weird reason, the string returned by the Value function
--   starts with the newline symbol which is not recognized as a space
--   by the Ada parser.

  vs : constant string := Value(s);
  p : Poly;

begin
  put("The number of variables : "); put(n,1); new_line;
  put("The string on input : "); put(vs); new_line;
 -- remove comments below to unveil the first symbol of vs :
 -- put("vs'first : "); put(vs'first,1); new_line;
 -- put("vs("); put(vs'first,1); put(") = "); put(vs(vs'first)); new_line;
 -- put("vs("); put(vs'first+1,1); put(") = "); put(vs(vs'first+1)); new_line;
  Symbol_Table.Init(natural32(n));
  p := Parse(natural32(n),vs(vs'first+1..vs'last));
  put("The parsed polynomial : ");
  put(p);
  new_line;
end print_poly;
