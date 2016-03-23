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

  function Underflow_to_Zero ( f : double_double ) return double_double is

    one : constant double_double := create(1.0);
    zero : constant double_double := create(0.0);

  begin
    if f + one = one 
     then return zero;
     else return f;
    end if;
  exception 
    when others => return zero;
  end Underflow_to_Zero;

  function Underflow_to_Zero ( f : quad_double ) return quad_double is

    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    if f + one = one 
     then return zero;
     else return f;
    end if;
  exception 
    when others => return zero;
  end Underflow_to_Zero;

  function Get_Real_Part
             ( x : Standard_Complex_Numbers.Complex_Number ) 
             return double_float is

  -- DESCRIPTION :
  --   Exception handlers wrap Standard_Complex_Numbers.REAL_PART,
  --   taking care of underflow.
 
    f : constant double_float := Standard_Complex_Numbers.REAL_PART(x);

  begin
    return Underflow_to_Zero(f);
  exception 
    when others => return 0.0;
  end Get_Real_Part;

  function Get_Real_Part
             ( x : DoblDobl_Complex_Numbers.Complex_Number ) 
             return double_double is

  -- DESCRIPTION :
  --   Exception handlers wrap DoblDobl_Complex_Numbers.REAL_PART,
  --   taking care of underflow.
 
    zero : constant double_double := create(0.0);
    f : constant double_double := DoblDobl_Complex_Numbers.REAL_PART(x);

  begin
    return Underflow_to_Zero(f);
  exception 
    when others => return zero;
  end Get_Real_Part;

  function Get_Real_Part
             ( x : QuadDobl_Complex_Numbers.Complex_Number ) 
             return quad_double is

  -- DESCRIPTION :
  --   Exception handlers wrap QuadDobl_Complex_Numbers.REAL_PART,
  --   taking care of underflow.
 
    zero : constant quad_double := create(0.0);
    f : constant quad_double := QuadDobl_Complex_Numbers.REAL_PART(x);

  begin
    return Underflow_to_Zero(f);
  exception 
    when others => return zero;
  end Get_Real_Part;

  function Get_Imaginary_Part
             ( x : Standard_Complex_Numbers.Complex_Number )
             return double_float is

  -- DESCRIPTION :
  --   Exception handlers wrap Standard_Complex_Numbers.IMAG_PART,
  --   taking care of underflow.

    f : constant double_float := Standard_Complex_Numbers.IMAG_PART(x);

  begin
    return Underflow_to_Zero(f);
  exception 
    when others => return 0.0;
  end Get_Imaginary_Part;

  function Get_Imaginary_Part
             ( x : DoblDobl_Complex_Numbers.Complex_Number )
             return double_double is

  -- DESCRIPTION :
  --   Exception handlers wrap DoblDobl_Complex_Numbers.IMAG_PART,
  --   taking care of underflow.

    zero : constant double_double := create(0.0);
    f : constant double_double := DoblDobl_Complex_Numbers.IMAG_PART(x);

  begin
    return Underflow_to_Zero(f);
  exception 
    when others => return zero;
  end Get_Imaginary_Part;

  function Get_Imaginary_Part
             ( x : QuadDobl_Complex_Numbers.Complex_Number )
             return quad_double is

  -- DESCRIPTION :
  --   Exception handlers wrap QuadDobl_Complex_Numbers.IMAG_PART,
  --   taking care of underflow.

    zero : constant quad_double := create(0.0);
    f : constant quad_double := QuadDobl_Complex_Numbers.IMAG_PART(x);

  begin
    return Underflow_to_Zero(f);
  exception 
    when others => return zero;
  end Get_Imaginary_Part;

  procedure Underflow_to_Zero
              ( x : in out Standard_Complex_Numbers.Complex_Number ) is

    re : constant double_float := Get_Real_Part(x);
    im : constant double_float := Get_Imaginary_Part(x);

  begin
    x := Standard_Complex_Numbers.Create(re,im);
  exception
    when others => put_line("exception when creating a complex number");
                   raise;
  end Underflow_to_Zero;

  procedure Underflow_to_Zero
              ( x : in out DoblDobl_Complex_Numbers.Complex_Number ) is

    re : constant double_double := Get_Real_Part(x);
    im : constant double_double := Get_Imaginary_Part(x);

  begin
    x := DoblDobl_Complex_Numbers.Create(re,im);
  exception
    when others => put_line("exception when creating a complex number");
                   raise;
  end Underflow_to_Zero;

  procedure Underflow_to_Zero
              ( x : in out QuadDobl_Complex_Numbers.Complex_Number ) is

    re : constant quad_double := Get_Real_Part(x);
    im : constant quad_double := Get_Imaginary_Part(x);

  begin
    x := QuadDobl_Complex_Numbers.Create(re,im);
  exception
    when others => put_line("exception when creating a complex number");
                   raise;
  end Underflow_to_Zero;

  procedure Underflow_to_Zero
              ( x : in out Standard_Complex_Vectors.Vector ) is
  begin
    for i in x'range loop
      Underflow_to_Zero(x(i));
    end loop;
  exception
    when others => put_line("exception raised in Underflow_to_Zero"); raise;
  end Underflow_to_Zero;

  procedure Underflow_to_Zero
              ( x : in out DoblDobl_Complex_Vectors.Vector ) is
  begin
    for i in x'range loop
      Underflow_to_Zero(x(i));
    end loop;
  exception
    when others => put_line("exception raised in Underflow_to_Zero"); raise;
  end Underflow_to_Zero;

  procedure Underflow_to_Zero
              ( x : in out QuadDobl_Complex_Vectors.Vector ) is
  begin
    for i in x'range loop
      Underflow_to_Zero(x(i));
    end loop;
  exception
    when others => put_line("exception raised in Underflow_to_Zero"); raise;
  end Underflow_to_Zero;

end Handle_Underflow_Gracefully;
