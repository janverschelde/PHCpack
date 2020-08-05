with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Solution_Strings;
with Total_Degree_Start_Systems;        use Total_Degree_Start_Systems;
with Lexicographic_Root_Enumeration;    use Lexicographic_Root_Enumeration;
with Drivers_to_Track_Standard_Paths;   use Drivers_to_Track_Standard_Paths;
with PHCpack_Operations;
with File_Management;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with Standard_Solutions_Interface;
with DoblDobl_Solutions_Interface;
with QuadDobl_Solutions_Interface;
with Multprec_Solutions_Interface;

function use_solcon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Job10 return integer32 is -- prompts for input file and opens it
  begin
    File_Management.Open_Input_File;
    return 0;
  end Job10;

  function Job11 return integer32 is -- prompts for output file and creates it
  begin
    File_Management.Create_Output_File;
    return 0;
  end Job11;

  function Job12 return integer32 is -- scans input file for "SOLUTIONS"

    found : boolean;

  begin
    Scan_and_Skip(File_Management.Link_to_Input.all,"SOLUTIONS",found);
    if found
     then return 0;
     else return 132;
    end if;
  exception
    when others => return 132;
  end Job12;

  function Job13 return integer32 is -- reads length and dimensions

    use Standard_Complex_Solutions_io;
    len,dim : natural32;

  begin
    Read_First(File_Management.Link_to_Input.all,len,dim);
    Assign(integer32(len),a);
    Assign(integer32(dim),b);
    return 0;
  exception
    when others => Assign(0,a); Assign(0,b); return 133;
  end Job13;

  function Job14 return integer32 is -- writes length and dimensions

    use Standard_Complex_Solutions_io;
    len,dim : natural32;

  begin
    Assign(a,integer32(len));
    Assign(b,integer32(dim));
    Write_First(File_Management.Link_to_Output.all,len,dim);
    return 0;
  exception
    when others => return 134;
  end Job14;

  function Job15 return integer32 is -- reads next solution from file

    use Standard_Complex_Solutions,Standard_Complex_Solutions_io;
    dim : natural32;
    ls : Link_to_Solution;

  begin
    Assign(a,integer32(dim));
   -- put("Dimension : "); put(dim,1); put_line(", calling Read_Next ...");
    Read_Next(File_Management.Link_to_Input.all,dim,ls);
   -- put_line("The solution read : "); put(ls.all); new_line;
    Assign_Solution(ls,b,c);
    Clear(ls);
    return 0;
  exception
    when others => return 135;
  end Job15;

  function Job16 return integer32 is -- writes next solution to file

    use Standard_Complex_Solutions,Standard_Complex_Solutions_io;
    cnt : natural32;
    ls : Link_to_Solution := Convert_to_Solution(b,c);

  begin
    Assign(a,integer32(cnt));
    Write_Next(File_Management.Link_to_Output.all,cnt,ls);
   -- put_line("Written solution : "); put(ls.all); new_line;
    Assign(integer32(cnt),a);
    Clear(ls);
    return 0;
  exception
    when others => return 136;
  end Job16;

  function Job17 return integer32 is -- close a solution input file

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));

  begin
    if k = 0
     then File_Management.Close_Input_File;
     else File_Management.Close_Input_File(k);
    end if;
    return 0;
  exception 
    when others => return 137;
  end Job17;

  function Job19 return integer32 is -- writes solutions banner to defined file
  begin
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE SOLUTIONS");
    else
      new_line(standard_output);
      put_line(standard_output,"THE SOLUTIONS");
    end if;
    return 0;
  end Job19;

  function Job20 return integer32 is -- writes sols dimensions to defined file

    use Standard_Complex_Solutions_io;
    len,dim : natural32;

  begin
    Assign(a,integer32(len));
    Assign(b,integer32(dim));
    if PHCpack_Operations.Is_File_Defined
     then Write_First(PHCpack_Operations.output_file,len,dim);
     else Write_First(standard_output,len,dim);
    end if;
    return 0;
  exception
    when others => return 140;
  end Job20;

  function Job21 return integer32 is -- writes next solution to defined file

    use Standard_Complex_Solutions,Standard_Complex_Solutions_io;
    cnt : natural32;
    ls : Link_to_Solution := Convert_to_Solution(b,c);

  begin
    Assign(a,integer32(cnt));
    if PHCpack_Operations.Is_File_Defined
     then Write_Next(PHCpack_Operations.output_file,cnt,ls);
     else Write_Next(standard_output,cnt,ls);
    end if;
   -- put_line("Written solution : "); put(ls.all); new_line;
    Assign(integer32(cnt),a);
    Clear(ls);
    return 0;
  exception
    when others => return 141;
  end Job21;

  function Job22 return integer32 is -- computes next total degree solution

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant natural32 := natural32(v_a(v_a'first));
    i : constant natural32 := natural32(v_a(v_a'first+1));
    lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    d,s : Standard_Natural_Vectors.Vector(1..integer32(n));
    cff,sol : Standard_Complex_Vectors.Vector(1..integer32(n));
    ls : Link_to_Solution;

  begin
    PHCpack_Operations.Retrieve_Start_System(lq);
    d := Degrees(lq.all);
    cff := Coefficients(lq.all);
    s := Root_Map(n,i,d);
    sol := Root(d,s,cff);
    ls := new Solution'(Create(sol));
    Assign_Solution(ls,b,c);
    Clear(ls);
    return 0;
  exception
    when others => return 142;
  end Job22;

  function Job23 return integer32 is -- computes next linear product solution

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant natural32 := natural32(v_a(v_a'first));
    i : constant natural32 := natural32(v_a(v_a'first+1));
    lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    d : Standard_Natural_Vectors.Vector(1..integer32(n));
   -- cp : Standard_Natural_Vectors.Vector(1..n-1);
    ls : Link_to_Solution;
    cnt,len : natural32;
    tol : constant double_float := 1.0E-10;
    fail : boolean;
    new_a : Standard_Natural_Vectors.Vector(1..2);

  begin
    PHCpack_Operations.Retrieve_Start_System(lq);
    d := Degrees(lq.all);
    len := Product(d);
   -- cp := Consecutive_Products(d);
    cnt := i;
   -- Next_Lex_Linear_Product_Start_Solution
   --   (false,n,d,cp,cnt,len,5,tol,ls,fail);
    Get_Next_Linear_Product_Start_Solution(false,n,cnt,len,5,tol,ls,fail);
    new_a(1) := n;   
    new_a(2) := cnt;
    Assign(new_a,a);
    if fail then
      return 143;
    else
      Assign_Solution(ls,b,c);
      Clear(ls);
      return 0;
    end if;
  exception
    when others => return 143;
  end Job23;

  function Job24 return integer32 is -- retrieve next linear product solution

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant natural32 := natural32(v_a(v_a'first));
    i : constant natural32 := natural32(v_a(v_a'first+1));
    lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    d : Standard_Natural_Vectors.Vector(1..integer32(n));
    cp : Standard_Natural_Vectors.Vector(1..integer32(n)-1);
    ls : Link_to_Solution;
    tol : constant double_float := 1.0E-10;
    fail : boolean;

  begin
    PHCpack_Operations.Retrieve_Start_System(lq);
    d := Degrees(lq.all);
    cp := Consecutive_Products(d);
    Get_Lex_Linear_Product_Start_Solution(false,n,d,cp,i,5,tol,ls,fail);
    if fail then
      return 144;
    else
      Assign_Solution(ls,b,c);
      Clear(ls);
      return 0;
    end if;
  exception
    when others => return 144;
  end Job24;

  function Job25 return integer32 is -- reads next witness point 

    use Standard_Complex_Solutions,Standard_Complex_Solutions_io;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;

  begin
   -- put("reading next witness point from set "); put(k,1); new_line;
   -- put("  solution vector has length "); put(n,1); new_line;
    Read_Next(File_Management.Link_to_Input(k).all,n,ls,
              Standard_Solutions_Container.Retrieve_Symbol_Table(0).all);
   -- was the following:
   --           Standard_Solutions_Container.Retrieve_Symbol_Table(k).all);
   -- put_line("The solution vector read : "); put(ls.all); new_line;
    Assign_Solution(ls,b,c);
    Clear(ls);
    return 0;
  exception 
    when others
      => put("Reading next witness point from set "); put(k,1);
         put_line(" raised exception."); -- raise;
         return 145;
  end Job25;

  function Job32 return integer32 is -- returns size of solution intro

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Link_to_Solution;
    fail : boolean;
    n : natural32;

  begin
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 202;
    else
      n := Standard_Solution_Strings.Length_Intro(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => return 202;
  end Job32;

  function Job33 return integer32 is -- returns solution intro string

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    s : string(1..integer(n));
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 203;
    else
      s := Standard_Solution_Strings.Write_Intro(ls.all);
      sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => return 203;
  end Job33;

  function Job34 return integer32 is -- returns size of solution vector

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Link_to_Solution;
    fail : boolean;
    n : natural32;

  begin
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 204;
    else
      n := Standard_Solution_Strings.Length_Vector(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => return 204;
  end Job34;

  function Job35 return integer32 is -- returns solution vector string

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    s : string(1..integer(n));
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 205;
    else
      s := Standard_Solution_Strings.Write_Vector(ls.all);
      sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => return 205;
  end Job35;

  function Job36 return integer32 is -- returns size of solution diagnostics

    n : natural32;

  begin
    n := Standard_Solution_Strings.Length_Diagnostics;
    Assign(integer32(n),b);
    return 0;
  exception
    when others => return 206;
  end Job36;

  function Job37 return integer32 is -- returns solution vector diagnostics

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    s : string(1..integer(n));
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 207;
    else
      s := Standard_Solution_Strings.Write_Diagnostics(ls.all);
      sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => return 207;
  end Job37;

  function Job875 return integer32 is -- set t = 0 for standard solutions

    use Standard_Complex_Solutions;
    sols : Solution_List := Standard_Solutions_Container.Retrieve;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    if not Is_Null(sols)
     then Set_Continuation_Parameter(sols,zero);
    end if;
    return 0;
  end Job875;

  function Job876 return integer32 is -- set t = 0 for dobldobl solutions

    use DoblDobl_Complex_Solutions;
    sols : Solution_List := DoblDobl_Solutions_Container.Retrieve;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer32(0));

  begin
    if not Is_Null(sols)
     then Set_Continuation_Parameter(sols,zero);
    end if;
    return 0;
  end Job876;

  function Job877 return integer32 is -- set t = 0 for quaddobl solutions

    use QuadDobl_Complex_Solutions;
    sols : Solution_List := QuadDobl_Solutions_Container.Retrieve;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer32(0));

  begin
    if not Is_Null(sols)
     then Set_Continuation_Parameter(sols,zero);
    end if;
    return 0;
  end Job877;

  function Handle_Jobs return integer32 is

    use Standard_Solutions_Interface;
    use DoblDobl_Solutions_Interface;
    use QuadDobl_Solutions_Interface;
    use Multprec_Solutions_Interface;

  begin
    case job is
      when 0 => return Standard_Solutions_Read(vrblvl);
      when 1 => return Standard_Solutions_Write(vrblvl);
      when 2 => return Standard_Solutions_Size(b,vrblvl);
      when 3 => return Standard_Solutions_Dimension(b,vrblvl);
      when 4 => return Standard_Solutions_Get(a,b,c,vrblvl);
      when 5 => return Standard_Solutions_Update(a,b,c,vrblvl);
      when 6 => return Standard_Solutions_Add(b,c,vrblvl);
      when 7 => return Standard_Solutions_Clear(vrblvl);
      when 8 => return Standard_Solutions_Drop_by_Index(a,vrblvl);
      when 9 => return Standard_Solutions_Drop_by_Name(a,b,vrblvl);
      when 10 => return Job10; -- prompts for input file and opens it
      when 11 => return Job11; -- prompts for output file and creates it
      when 12 => return Job12; -- scans input for "SOLUTIONS"
      when 13 => return Job13; -- reads length and dimensions
      when 14 => return Job14; -- writes length and dimensions
      when 15 => return Job15; -- reads next solution from file
      when 16 => return Job16; -- writes next solution to file
      when 17 => return Job17; -- close a solution input file
      when 18 => File_Management.Close_Output_File; return 0;
      when 19 => return Job19; -- writes solutions banner to defined file
      when 20 => return Job20; -- writes solution dimensions defined file
      when 21 => return Job21; -- writes next solution to defined output file
      when 22 => return Job22; -- computes next total solution
      when 23 => return Job23; -- computes next linear product solution
      when 24 => return Job24; -- retrieves next linear product solution
      when 25 => return Job25; -- reads next witness point from file
      when 30 => return Standard_Solutions_String_Size(a,b,vrblvl);
      when 31 => return Standard_Solutions_Get_String(a,b,vrblvl);
      when 32 => return Job32; -- returns size of solution intro
      when 33 => return Job33; -- returns solution intro string
      when 34 => return Job34; -- returns size of solution vector
      when 35 => return Job35; -- returns solution vector string
      when 36 => return Job36; -- returns size of solution diagnostics
      when 37 => return Job37; -- returns solution diagnostics string
      when 38 => return Standard_Solutions_Add_String(a,b,vrblvl);
      when 39 => return Standard_Solutions_Replace_String(a,b,vrblvl);
     -- corresponding operations for double double solutions
      when 40 => return DoblDobl_Solutions_Read(vrblvl);
      when 41 => return DoblDobl_Solutions_Write(vrblvl);
      when 42 => return DoblDobl_Solutions_Size(b,vrblvl);
      when 43 => return DoblDobl_Solutions_Dimension(b,vrblvl);
      when 44 => return DoblDobl_Solutions_Get(a,b,c,vrblvl);
      when 45 => return DoblDobl_Solutions_Update(a,b,c,vrblvl);
      when 46 => return DoblDobl_Solutions_Add(b,c,vrblvl);
      when 47 => return DoblDobl_Solutions_Clear(vrblvl);
      when 48 => return DoblDobl_Solutions_Drop_by_Index(a,vrblvl);
      when 49 => return DoblDobl_Solutions_Drop_by_Name(a,b,vrblvl);
      when 70 => return DoblDobl_Solutions_String_Size(a,b,vrblvl);
      when 71 => return DoblDobl_Solutions_Get_String(a,b,vrblvl);
      when 78 => return DoblDobl_Solutions_Add_String(a,b,vrblvl);
     -- corresponding operations for quad double solutions
      when 80 => return QuadDobl_Solutions_Read(vrblvl);
      when 81 => return QuadDobl_Solutions_Write(vrblvl);
      when 82 => return QuadDobl_Solutions_Size(b,vrblvl);
      when 83 => return QuadDobl_Solutions_Dimension(b,vrblvl);
      when 84 => return QuadDobl_Solutions_Get(a,b,c,vrblvl);
      when 85 => return QuadDobl_Solutions_Update(a,b,c,vrblvl);
      when 86 => return QuadDobl_Solutions_Add(b,c,vrblvl);
      when 87 => return QuadDobl_Solutions_Clear(vrblvl);
      when 88 => return QuadDobl_Solutions_Drop_by_Index(a,vrblvl);
      when 89 => return QuadDobl_Solutions_Drop_by_Name(a,b,vrblvl);
      when 110 => return QuadDobl_Solutions_String_Size(a,b,vrblvl);
      when 111 => return QuadDobl_Solutions_Get_String(a,b,vrblvl);
      when 118 => return QuadDobl_Solutions_Add_String(a,b,vrblvl);
     -- corresponding operations for multiprecision solutions
      when 120 => return Multprec_Solutions_Read(vrblvl);
      when 121 => return Multprec_Solutions_Write(vrblvl);
      when 122 => return Multprec_Solutions_Size(b,vrblvl);
      when 123 => return Multprec_Solutions_Dimension(b,vrblvl);
      when 127 => return Multprec_Solutions_Clear(vrblvl);
      when 150 => return Multprec_Solutions_String_Size(a,b,vrblvl);
      when 151 => return Multprec_Solutions_Get_String(a,b,vrblvl);
      when 158 => return Multprec_Solutions_Add_String(a,b,vrblvl);
     -- retrieve next solution
      when 276 => return Standard_Solutions_Retrieve_Next(a,b,c,vrblvl);
      when 277 => return DoblDobl_Solutions_Retrieve_Next(a,b,c,vrblvl);
      when 278 => return QuadDobl_Solutions_Retrieve_Next(a,b,c,vrblvl);
      when 279 => return Multprec_Solutions_Set_Pointer(vrblvl);
     -- move pointer from current to next solution
      when 300 => return Standard_Solutions_Move_Pointer(a,vrblvl);
      when 301 => return DoblDobl_Solutions_Move_Pointer(a,vrblvl);
      when 302 => return QuadDobl_Solutions_Move_Pointer(a,vrblvl);
      when 303 => return Multprec_Solutions_Move_Pointer(a,vrblvl);
     -- return length of current solution string
      when 304 => return Standard_Solutions_Current_Size(a,b,vrblvl);
      when 305 => return DoblDobl_Solutions_Current_Size(a,b,vrblvl);
      when 306 => return QuadDobl_Solutions_Current_Size(a,b,vrblvl);
      when 307 => return Multprec_Solutions_Current_Size(a,b,vrblvl);
     -- return current solution string
      when 308 => return Standard_Solutions_Current_String(a,b,vrblvl);
      when 309 => return DoblDobl_Solutions_Current_String(a,b,vrblvl);
      when 310 => return QuadDobl_Solutions_Current_String(a,b,vrblvl);
      when 311 => return Multprec_Solutions_Current_String(a,b,vrblvl);
     -- reading system and solutions from given file name
      when 544 => return Standard_System_Solutions_Read_from_File(a,b,vrblvl);
      when 545 => return DoblDobl_System_Solutions_Read_from_File(a,b,vrblvl);
      when 546 => return QuadDobl_System_Solutions_Read_from_File(a,b,vrblvl);
      when 547 => return Multprec_System_Solutions_Read_from_File(a,b,vrblvl);
     -- setting value of continuation parameter of the solutions to zero
      when 875 => return Job875; -- for standard solutions, set t = 0
      when 876 => return Job876; -- for dobldobl solutions, set t = 0
      when 877 => return Job877; -- for quaddobl solutions, set t = 0
     -- add one to each solution in 1-homogeneous coordinate transformation
      when 894 => return Standard_Solutions_Make_Homogeneous(vrblvl);
      when 895 => return DoblDobl_Solutions_Make_Homogeneous(vrblvl);
      when 896 => return QuadDobl_Solutions_Make_Homogeneous(vrblvl);
     -- add ones to each solution in m-homogeneous coordinate transformation
      when 910 => return Standard_Solutions_Multi_Homogeneous(a,vrblvl);
      when 911 => return DoblDobl_Solutions_Multi_Homogeneous(a,vrblvl);
      when 912 => return QuadDobl_Solutions_Multi_Homogeneous(a,vrblvl);
     -- affine transformations of the solutions
      when 898 => return Standard_Solutions_1Hom2Affine(vrblvl);
      when 899 => return DoblDobl_Solutions_1Hom2Affine(vrblvl);
      when 900 => return QuadDobl_Solutions_1Hom2Affine(vrblvl);
      when 913 => return Standard_Solutions_mHom2Affine(a,b,vrblvl);
      when 914 => return DoblDobl_Solutions_mHom2Affine(a,b,vrblvl);
      when 915 => return QuadDobl_Solutions_mHom2Affine(a,b,vrblvl);
      when 916 => return Standard_Solutions_Read_from_File(a,b,vrblvl);
      when 917 => return Dobldobl_Solutions_Read_from_File(a,b,vrblvl);
      when 918 => return QuadDobl_Solutions_Read_from_File(a,b,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_solcon;
