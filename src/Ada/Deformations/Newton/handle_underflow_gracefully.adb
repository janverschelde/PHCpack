with text_io; use text_io;

package body Handle_Underflow_Gracefully is

  function Underflow_to_Zero ( f : double_float ) return double_float is
  begin
    if f + 1.0 = 1.0
     then return 0.0;
     else return f;
    end if;
  exception 
    when others => return 0.0;
  end Underflow_to_Zero;

  function Get_Real_Part ( x : Complex_Number ) return double_float is

  -- DESCRIPTION :
  --   Exception handlers wrap Standard_Complex_Numbers.REAL_PART,
  --   taking care of underflow.
 
    f : constant double_float := REAL_PART(x);

  begin
    return Underflow_to_Zero(f);
  exception 
    when others => return 0.0;
  end Get_Real_Part;

  function Get_Imaginary_Part ( x : Complex_Number ) return double_float is

  -- DESCRIPTION :
  --   Exception handlers wrap Standard_Complex_Numbers.IMAG_PART,
  --   taking care of underflow.

    f : constant double_float := IMAG_PART(x);

  begin
    return Underflow_to_Zero(f);
  exception 
    when others => return 0.0;
  end Get_Imaginary_Part;

  procedure Underflow_to_Zero ( x : in out Complex_Number ) is

    re : constant double_float := Get_Real_Part(x);
    im : constant double_float := Get_Imaginary_Part(x);

  begin
    x := Create(re,im);
  exception
    when others => put_line("exception when creating a complex number");
                   raise;
  end Underflow_to_Zero;

  procedure Underflow_to_Zero ( x : in out Vector ) is
  begin
    for i in x'range loop
      Underflow_to_Zero(x(i));
    end loop;
  exception
    when others => put_line("exception raised in Underflow_to_Zero"); raise;
  end Underflow_to_Zero;

end Handle_Underflow_Gracefully;
