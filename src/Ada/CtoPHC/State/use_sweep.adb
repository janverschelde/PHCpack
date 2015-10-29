with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Symbol_Table;
-- with Symbol_Table_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_Integer_Vectors;
with Parameter_Homotopy_State;

function use_sweep ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is -- define parameters numerically

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nbequ : constant integer32 := integer32(v_a(v_a'first));
    nbvar : constant integer32 := integer32(v_a(v_a'first+1));
    nbpar : constant integer32 := integer32(v_a(v_a'first+2));
    idx : Standard_Integer_Vectors.Vector(1..nbpar);

  begin
    Parameter_Homotopy_State.Set_Number_of_Equations(nbequ);
    Parameter_Homotopy_State.Set_Number_of_Variables(nbvar);
    Parameter_Homotopy_State.Set_Number_of_Parameters(nbpar);
    Assign(natural32(nbpar),b,idx);
    Parameter_Homotopy_State.Set_Indices(idx);
    return 0;
  end Job0;

  function Job1 return integer32 is -- define parameters symbolically

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    nbequ : constant integer32 := integer32(v_a(v_a'first));
    nbvar : constant integer32 := integer32(v_a(v_a'first+1));
    nbpar : constant integer32 := integer32(v_a(v_a'first+2));
    nbchr : constant natural32 := natural32(v_a(v_a'first+3));
    idx : Standard_Integer_Vectors.Vector(1..nbpar);
    v_b : constant C_Integer_Array 
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(nbchr));
    s : constant string
      := C_Integer_Array_to_String(nbchr,v_b);
    ind : integer;

  begin
    Parameter_Homotopy_State.Set_Number_of_Equations(nbequ);
    Parameter_Homotopy_State.Set_Number_of_Variables(nbvar);
    Parameter_Homotopy_State.Set_Number_of_Parameters(nbpar);
   -- put_line("The string in the interface : " & s);
    ind := s'first-1;
    for k in 1..nbpar loop
      ind := ind + 1;
      declare
        sb : Symbol_Table.Symbol;
        sx : integer := sb'first-1;
        idxsb : natural32;
      begin
        sb := (sb'range => ' ');
        while ind <= s'last loop
          exit when (s(ind) = ' ');
          sx := sx + 1;
          sb(sx) := s(ind);
          ind := ind + 1;
        end loop;
       -- put("Read symbol : "); Symbol_Table_io.put(sb); new_line;
        idxsb := Symbol_Table.Get(sb);
        if idxsb > 0
         then idx(k) := integer32(idxsb);
        end if;
      end;
    end loop;
    Parameter_Homotopy_State.Set_Indices(idx);
    return 0;
  end Job1;

  function Job2 return integer32 is -- return number of equations

    n : constant integer32
      :=  Parameter_Homotopy_State.Get_Number_of_Equations;
  
  begin
    Assign(n,a);
    return 0;
  end Job2;

  function Job3 return integer32 is -- return number of variables

    n : constant integer32
      :=  Parameter_Homotopy_State.Get_Number_of_Variables;
  
  begin
    Assign(n,a);
    return 0;
  end Job3;

  function Job4 return integer32 is -- return number of parameters

    n : constant integer32
      :=  Parameter_Homotopy_State.Get_Number_of_Parameters;
  
  begin
    Assign(n,a);
    return 0;
  end Job4;

  function Job5 return integer32 is -- return parameters numerically

    idx : constant Standard_Integer_Vectors.Link_to_Vector
        := Parameter_Homotopy_State.Get_Indices;

    use Standard_Integer_Vectors;

  begin
    if idx /= null
     then Assign(idx.all,a);
    end if;
    return 0;
  end Job5;

  function Job6 return integer32 is -- return parameters symbolically

    idx : constant Standard_Integer_Vectors.Link_to_Vector
        := Parameter_Homotopy_State.Get_Indices;

    use Standard_Integer_Vectors;
    use String_Splitters;

    ls : Link_to_String;
    sb : Symbol_Table.Symbol;
    len : natural32;

  begin
    if idx /= null then
      for i in idx'range loop
        sb := Symbol_Table.Get(natural32(idx(i)));
       -- put("Retrieved symbol : "); Symbol_Table_io.put(sb); new_line;
        len := Symbol_Table.Length(sb);
        declare
          strsb : string(1..integer(len+1));
        begin
          for k in 1..len loop
            strsb(integer(k)) := sb(integer(k));
          end loop;
          strsb(strsb'last) := ' ';
          Append(ls,strsb);
        end;
      end loop;
      len := natural32(ls'last) - 1; -- skip the last space
      Assign(integer32(len),a);
      if len > 0 then
        declare
          v : constant Standard_Integer_Vectors.Vector
            := String_to_Integer_Vector(ls(ls'first..integer(len)));
        begin
          Assign(v,b);
        end;
      end if;
    end if;
    return 0;
  end Job6;

  function Job7 return integer32 is -- clear parameter definitions
  begin
    Parameter_Homotopy_State.Clear;
    return 0;
  end Job7;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- define parameters numerically
      when 1 => return Job1; -- define parameters symbolically
      when 2 => return Job2; -- return the number of equations
      when 3 => return Job3; -- return the number of variables
      when 4 => return Job4; -- return the number of parameters
      when 5 => return Job5; -- return parameters numerically
      when 6 => return Job6; -- return parameters symbolically
      when 7 => return Job7; -- clear parameter definitions
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_sweep;
