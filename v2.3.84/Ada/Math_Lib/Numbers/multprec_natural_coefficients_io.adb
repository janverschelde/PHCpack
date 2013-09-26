with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;

package body Multprec_Natural_Coefficients_io is

  procedure put ( n : Array_of_Naturals ) is
  begin
    for i in reverse n'range loop
      if i < n'last
       then put(" + ");
      end if;
      put(n(i),8);
      put("*b^");
      put(i,1);
    end loop;
  end put;

end Multprec_Natural_Coefficients_io;
