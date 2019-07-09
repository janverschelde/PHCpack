with text_io;                           use text_io;
with Interfaces.C;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;

function use_tabform ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  procedure Extract_Dimensions ( neq,nvr,nbt : out integer32 ) is

  -- DESCRIPTION :
  --   Extracts the dimensions out of the three numbers of a.

  -- ON RETURN :
  --   neq       the number of equations is in a[0];
  --   nvr       the number of variables is in a[1];
  --   nbt       the total number of terms is in a[2].

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));

    use Interfaces.C;

  begin
    neq := integer32(v(v'first));
    nvr := integer32(v(v'first+1));
    nbt := integer32(v(v'first+2));
  end Extract_Dimensions;

  function Job0 return integer32 is -- store tableau form

    neq,nvr,nbt : integer32;

  begin
    Extract_Dimensions(neq,nvr,nbt);
    put("-> the number of equations : "); put(neq,1); new_line;
    put("-> the number of variables : "); put(nvr,1); new_line;
    put("-> total number of terms : "); put(nbt,1); new_line;
    return 0;
  end Job0;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0;
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_tabform;
