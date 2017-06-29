with Standard_Floating_Numbers;
with Standard_Complex_Numbers;

package body Valid_Vector_Checks is

  function is_valid ( v : Standard_Floating_Vectors.Vector ) return boolean is
  begin
    for i in v'range loop
      if not Standard_Floating_Numbers.is_valid(v(i))
       then return false;
      end if;
    end loop;
    return true;
  end is_valid;

  function is_valid ( v : Standard_Complex_Vectors.Vector ) return boolean is
  begin
    for i in v'range loop
      if not Standard_Complex_Numbers.is_valid(v(i))
       then return false;
      end if;
    end loop;
    return true;
  end is_valid;

end Valid_Vector_Checks;
