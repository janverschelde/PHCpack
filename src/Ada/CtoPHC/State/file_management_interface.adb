with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Interfaces.C;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with File_Management;

package body File_Management_Interface is

  function File_Management_Prompt_Input_File
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in file_management_interface.");
      put_line("File_Management_Prompt_Input_File.");
    end if;
    File_Management.Open_Input_File;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in file_management_interface.");
        put_line("File_Management_Prompt_Input_File.");
      end if;
      return 130;
  end File_Management_Prompt_Input_File;

  function File_Management_Prompt_Output_File
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in file_management_interface.");
      put_line("File_Management_Prompt_Output_File.");
    end if;
    File_Management.Create_Output_File;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in file_management_interface.");
        put_line("File_Management_Prompt_Output_File.");
      end if;
      return 131;
  end File_Management_Prompt_Output_File;

  function File_Management_Set_Output
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..integer(n)) := C_Integer_Array_to_String(n,v_b);

  begin
    if vrblvl > 0 then
      put("-> in file_management_interface.");
      put_line("File_Management_Set_Output.");
    end if;
    PHCpack_Operations.Define_Output_File(s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in file_management_interface.");
        put_line("File_Management_Set_Output.");
      end if;
      return 191;
  end File_Management_Set_Output;

  function File_Management_Close_Output
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in file_management_interface.");
      put_line("File_Management_Close_Output.");
    end if;
    PHCpack_Operations.Close_Output_File;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in file_management_interface.");
        put_line("File_Management_Close_Output.");
      end if;
      return 192;
  end File_Management_Close_Output;

  function File_Management_Write_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant integer := integer(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..n) := C_Integer_Array_to_String(natural32(n),v_b);

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("File_Management_Write_String.");
    end if;
    if PHCpack_Operations.Is_File_Defined
     then put(PHCpack_Operations.output_file,s);
     else put(standard_output,s);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in file_management_interface.");
        put_line("File_Management_Write_String.");
      end if;
      return 158;
  end File_Management_Write_String;

  function File_Management_Write_Integers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(n));

    use Interfaces.C;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("File_Management_Write_Integers.");
    end if;
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,integer32(v_b(v_b'first)),1);
      for i in v_b'first+1..v_b'last loop
        put(PHCpack_Operations.output_file," ");
        put(PHCpack_Operations.output_file,integer32(v_b(i)),1);
      end loop;
    else
      put(standard_output,integer32(v_b(v_b'first)),1);
      for i in v_b'first+1..v_b'last loop
        put(standard_output," ");
        put(standard_output,integer32(v_b(i)),1);
      end loop;
    end if;
    return 0;
  end File_Management_Write_Integers;

  function File_Management_Write_Doubles
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_c : constant C_Double_Array(0..n1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(n));
    d : double_float;

    use Interfaces.C;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("File_Management_Write_Doubles.");
    end if;
    d := double_float(v_c(v_c'first));
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,d);
      for i in v_c'first+1..v_c'last loop
        put(PHCpack_Operations.output_file," ");
        d := double_float(v_c(i));
        put(PHCpack_Operations.output_file,d);
      end loop;
    else
      put(standard_output,d);
      for i in v_c'first+1..v_c'last loop
        put(standard_output," ");
        d := double_float(v_c(i));
        put(standard_output,d);
      end loop;
    end if;
    return 0;
  end File_Management_Write_Doubles;

end File_Management_Interface;
