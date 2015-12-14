with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;
with Standard_Solution_Strings;
with DoblDobl_Solution_Strings;
with QuadDobl_Solution_Strings;
with Multprec_Solution_Strings;
with Symbol_Table;
with Solution_Drops;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Multprec_System_and_Solutions_io;
with Total_Degree_Start_Systems;        use Total_Degree_Start_Systems;
with Lexicographic_Root_Enumeration;    use Lexicographic_Root_Enumeration;
with Drivers_to_Track_Standard_Paths;   use Drivers_to_Track_Standard_Paths;
with PHCpack_Operations;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Multprec_PolySys_Container;
with File_Management;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Multprec_Solutions_Container;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;

function use_solcon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is -- read from file into container

    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Standard_Complex_Solutions_io.Read(sols);
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  end Job0;

  function Job40 return integer32 is -- read from file into container

    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    DoblDobl_Complex_Solutions_io.Read(sols);
    DoblDobl_Solutions_Container.Initialize(sols);
    return 0;
  end Job40;

  function Job80 return integer32 is -- read from file into container

    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    QuadDobl_Complex_Solutions_io.Read(sols);
    QuadDobl_Solutions_Container.Initialize(sols);
    return 0;
  end Job80;

  function Job120 return integer32 is -- read from file into container

    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    Multprec_Complex_Solutions_io.Read(sols);
    Multprec_Solutions_Container.Initialize(sols);
    return 0;
  end Job120;

  function Job1 return integer32 is -- write container to file

    use Standard_Complex_Solutions,Standard_Complex_Solutions_io;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if not Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put(PHCpack_Operations.output_file,
            Length_Of(sols),natural32(Head_Of(sols).n),sols);
      else
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
    return 0;
  end Job1;

  function Job41 return integer32 is -- write container to file

    use DoblDobl_Complex_Solutions,DoblDobl_Complex_Solutions_io;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if not Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put(PHCpack_Operations.output_file,
            Length_Of(sols),natural32(Head_Of(sols).n),sols);
      else
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
    return 0;
  end Job41;

  function Job81 return integer32 is -- write container to file

    use QuadDobl_Complex_Solutions,QuadDobl_Complex_Solutions_io;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if not Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put(PHCpack_Operations.output_file,
            Length_Of(sols),natural32(Head_Of(sols).n),sols);
      else
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
    return 0;
  end Job81;

  function Job121 return integer32 is -- write container to file

    use Multprec_Complex_Solutions,Multprec_Complex_Solutions_io;
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;

  begin
    if not Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put(PHCpack_Operations.output_file,
            Length_Of(sols),natural32(Head_Of(sols).n),sols);
      else
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
    return 0;
  end Job121;

  function Job4 return integer32 is -- return solution in container

    use Standard_Complex_Solutions;
    ls : Link_to_Solution;
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail
     then return 34;
     else Assign_Solution(ls,b,c); return 0;
    end if;
  end Job4;

  function Job44 return integer32 is -- return solution in container

    use DoblDobl_Complex_Solutions;
    ls : Link_to_Solution;
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    DoblDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail
     then return 44;
     else Assign_Solution(ls,b,c); return 0;
    end if;
  end Job44;

  function Job84 return integer32 is -- return solution in container

    use QuadDobl_Complex_Solutions;
    ls : Link_to_Solution;
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    QuadDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail
     then return 84;
     else Assign_Solution(ls,b,c); return 0;
    end if;
  end Job84;

  function Job5 return integer32 is -- change standard solution in container

    use Standard_Complex_Solutions;
    ls : Link_to_Solution := Convert_to_Solution(b,c);
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    Standard_Solutions_Container.Replace(k,ls.all,fail);
    Clear(ls);
    if fail
     then return 35;
     else return 0;
    end if;
  end Job5;

  function Job45 return integer32 is -- change dobldobl solution in container

    use DoblDobl_Complex_Solutions;
    ls : Link_to_Solution := Convert_to_Solution(b,c);
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    DoblDobl_Solutions_Container.Replace(k,ls.all,fail);
    Clear(ls);
    if fail
     then return 345;
     else return 0;
    end if;
  end Job45;

  function Job85 return integer32 is -- change quaddobl solution in container

    use QuadDobl_Complex_Solutions;
    ls : Link_to_Solution := Convert_to_Solution(b,c);
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    QuadDobl_Solutions_Container.Replace(k,ls.all,fail);
    Clear(ls);
    if fail
     then return 395;
     else return 0;
    end if;
  end Job85;

  function Job6 return integer32 is -- append standard solution to container

    use Standard_Complex_Solutions;
    ls : constant Link_to_Solution := Convert_to_Solution(b,c);

  begin
    Standard_Solutions_Container.Append(ls);
    return 0;
  end Job6;

  function Job46 return integer32 is -- append dobldobl solution to container

    use DoblDobl_Complex_Solutions;
    ls : constant Link_to_Solution := Convert_to_Solution(b,c);

  begin
    DoblDobl_Solutions_Container.Append(ls);
    return 0;
  end Job46;

  function Job86 return integer32 is -- append quaddobl solution to container

    use QuadDobl_Complex_Solutions;
    ls : constant Link_to_Solution := Convert_to_Solution(b,c);

  begin
    QuadDobl_Solutions_Container.Append(ls);
    return 0;
  end Job86;

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

  function Job30 return integer32 is -- returns size of solution string

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
      return 200;
    else
      n := Standard_Solution_Strings.Length(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => return 200;
  end Job30;

  function Job70 return integer32 is -- returns size of solution string

    use DoblDobl_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Link_to_Solution;
    fail : boolean;
    n : natural32;

  begin
    DoblDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 240;
    else
      n := DoblDobl_Solution_Strings.Length(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => return 240;
  end Job70;

  function Job110 return integer32 is -- returns size of solution string

    use QuadDobl_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Link_to_Solution;
    fail : boolean;
    n : natural32;

  begin
    QuadDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 280;
    else
      n := QuadDobl_Solution_Strings.Length(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => return 280;
  end Job110;

  function Job150 return integer32 is -- returns size of solution string

    use Multprec_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Link_to_Solution;
    fail : boolean;
    n : natural32;

  begin
    Multprec_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 280;
    else
      n := Multprec_Solution_Strings.Length(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => return 280;
  end Job150;

  function Job31 return integer32 is -- returns solution string

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 201;
    else
      declare
        s  : constant string := Standard_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      end;
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => return 201;
  end Job31;

  function Job71 return integer32 is -- returns solution string

    use DoblDobl_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    DoblDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 241;
    else
      declare
        s : constant string := DoblDobl_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      end;
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => return 241;
  end Job71;

  function Job111 return integer32 is -- returns solution string

    use QuadDobl_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    QuadDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 281;
    else
      declare
        s : constant string := QuadDobl_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      end;
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => return 281;
  end Job111;

  function Job151 return integer32 is -- returns solution string

    use Multprec_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    Multprec_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 201;
    else
      declare
        s : constant string := Multprec_Solution_Strings.Write(ls.all);
      begin
       -- put_line("The string in Job151 : " & s);
       -- put("  s'last = "); put(natural32(s'last),1);
       -- put("  n = "); put(n,1); new_line;
        sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      end;
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => return 281;
  end Job151;

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

  function Job38 return integer32 is -- append solution string to container

    use Standard_Complex_Solutions;
    va : constant C_Integer_Array 
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(va(va'first));
    nc : constant integer := integer(va(va'first+1));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : String(1..nc);   
    ind : integer := 1;
    sol : Solution(integer32(nv));
    fail : boolean;

  begin
   -- put_line("Inside the Job 38 ...");
    sv := C_Integer_Array_to_String(natural32(nc),vb);
   -- put_line("The string received in Job 38 : "); put_line(sv);
    Standard_Solution_Strings.Parse(sv,ind,nv,sol,fail);
   -- put_line("The parsed solution : ");
   -- Standard_Complex_Solutions_io.put(sol);
    if fail then
     -- put_line("Failure occurred !");
      return 208;
    else
     -- put_line("Appending the solution to the container...");
      Standard_Solutions_Container.Append(sol);
    end if;
    return 0;
  exception
    when others => 
 --     put_line("exception occurred in Job 38 of use_solcon, with nv :");
 --     put(nv,1); new_line;
 --     put_line("exception occurred in Job 38 of use_solcon, with nc :");
 --     put(natural32(nc),1); new_line;
 --     put_line("exception occurred in Job 38 of use_solcon, with string :");
 --     put_line(sv);
      return 208;
  end Job38;

  function Job78 return integer32 is -- append solution string to container

    use DoblDobl_Complex_Solutions;
    va : constant C_Integer_Array 
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(va(va'first));
    nc : constant integer := integer(va(va'first+1));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : String(1..nc);   
    ind : integer := 1;
    sol : Solution(integer32(nv));
    fail : boolean;

  begin
    sv := C_Integer_Array_to_String(natural32(nc),vb);
    DoblDobl_Solution_Strings.Parse(sv,ind,nv,sol,fail);
   -- put_line("the parsed solution : ");
   -- DoblDobl_Complex_Solutions_io.put(sol);
    if fail then
      return 208;
    else
      DoblDobl_Solutions_Container.Append(sol);
    end if;
    return 0;
  exception
    when others => return 248;
  end Job78;

  function Job118 return integer32 is -- append solution string to container

    use QuadDobl_Complex_Solutions;
    va : constant C_Integer_Array 
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(va(va'first));
    nc : constant integer := integer(va(va'first+1));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : String(1..nc);   
    ind : integer := 1;
    sol : Solution(integer32(nv));
    fail : boolean;

  begin
    sv := C_Integer_Array_to_String(natural32(nc),vb);
    QuadDobl_Solution_Strings.Parse(sv,ind,nv,sol,fail);
    if fail then
      return 288;
    else
      QuadDobl_Solutions_Container.Append(sol);
    end if;
    return 0;
  exception
    when others => return 288;
  end Job118;

  function Job158 return integer32 is -- append solution string to container

    use Multprec_Complex_Solutions;
    va : constant C_Integer_Array 
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(va(va'first));
    nc : constant integer := integer(va(va'first+1));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : String(1..nc);   
    ind : integer := 1;
    sol : Solution(integer32(nv));
    fail : boolean;

  begin
    sv := C_Integer_Array_to_String(natural32(nc),vb);
    Multprec_Solution_Strings.Parse(sv,ind,nv,sol,fail);
    if fail then
      return 288;
    else
      Multprec_Solutions_Container.Append(sol);
    end if;
    return 0;
  exception
    when others => return 288;
  end Job158;

  function Job39 return integer32 is -- replace solution string in container

    use Standard_Complex_Solutions;
    va : constant C_Integer_Array 
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    k : constant natural32 := natural32(va(va'first));
    nv : constant natural32 := natural32(va(va'first+1));
    nc : constant integer := integer(va(va'first+2));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : String(1..nc);   
    ind : integer := 1;
    sol : Solution(integer32(nv));
    fail : boolean;

  begin
    sv := C_Integer_Array_to_String(natural32(nc),vb);
    Standard_Solution_Strings.Parse(sv,ind,nv,sol,fail);
    if fail then
      return 209;
    else
      Standard_Solutions_Container.Replace(k,sol,fail);
      if fail
       then return 209;
       else return 0;
      end if;
    end if;
  exception
    when others => return 209;
  end Job39;

-- The jobs to drop a coordinate of a solution come in two flavors:
-- (1) by index: given the index of the variable in a[0];
-- (2) by name: given the number of characters of the symbol in a[0]
--     and the characters for the symbol name in b.

  function Job8 return integer32 is -- drop by index from solution list

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant natural32 := natural32(v_a(v_a'first));
    use Standard_Complex_Solutions;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    dropped : constant Solution_List := Solution_Drops.Drop(sols,ind);

  begin
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(dropped);
    return 0;
  end Job8;

  function Job9 return integer32 is -- drop by name from solution list

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant integer := integer(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use Standard_Complex_Solutions;
    sols,dropped : Solution_List;

  begin
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    sols := Standard_Solutions_Container.Retrieve;
    dropped := Solution_Drops.Drop(sols,ind);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(dropped);
    return 0;
  end Job9;

  function Job48 return integer32 is -- drop by index from solution list

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant natural32 := natural32(v_a(v_a'first));
    use DoblDobl_Complex_Solutions;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    dropped : constant Solution_List := Solution_Drops.Drop(sols,ind);

  begin
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(dropped);
    return 0;
  end Job48;

  function Job49 return integer32 is -- drop by name from solution list

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant integer := integer(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use DoblDobl_Complex_Solutions;
    sols,dropped : Solution_List;

  begin
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    sols := DoblDobl_Solutions_Container.Retrieve;
    dropped := Solution_Drops.Drop(sols,ind);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(dropped);
    return 0;
  end Job49;

  function Job88 return integer32 is -- drop by index from solution list

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant natural32 := natural32(v_a(v_a'first));
    use QuadDobl_Complex_Solutions;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    dropped : constant Solution_List := Solution_Drops.Drop(sols,ind);

  begin
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(dropped);
    return 0;
  end Job88;

  function Job89 return integer32 is -- drop by name from solution list

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant integer := integer(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use QuadDobl_Complex_Solutions;
    sols,dropped : Solution_List;

  begin
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    sols := QuadDobl_Solutions_Container.Retrieve;
    dropped := Solution_Drops.Drop(sols,ind);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(dropped);
    return 0;
  end Job89;

  function Job276 return integer32 is -- retrieve next standard solution

    use Standard_Complex_Solutions;
    ls : Link_to_Solution;
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    idx : natural32;

  begin
    if k = 0 then
      Standard_Solutions_Container.Retrieve_Next_Initialize;
    else
      Standard_Solutions_Container.Retrieve_Next(ls,idx);
      Assign(integer32(idx),a);
      if idx = 0
       then return 276;
       else Assign_Solution(ls,b,c);
      end if;
    end if;
    return 0;
  end Job276;

  function Job277 return integer32 is -- retrieve next dobldobl solution

    use DoblDobl_Complex_Solutions;
    ls : Link_to_Solution;
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    idx : natural32;

  begin
    if k = 0 then
      DoblDobl_Solutions_Container.Retrieve_Next_Initialize;
    else
      DoblDobl_Solutions_Container.Retrieve_Next(ls,idx);
      Assign(integer32(idx),a);
      if idx = 0
       then return 276;
       else Assign_Solution(ls,b,c);
      end if;
    end if;
    return 0;
  end Job277;

  function Job278 return integer32 is -- retrieve next quaddobl solution

    use QuadDobl_Complex_Solutions;
    ls : Link_to_Solution;
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    idx : natural32;

  begin
    if k = 0 then
      QuadDobl_Solutions_Container.Retrieve_Next_Initialize;
    else
      QuadDobl_Solutions_Container.Retrieve_Next(ls,idx);
      Assign(integer32(idx),a);
      if idx = 0
       then return 276;
       else Assign_Solution(ls,b,c);
      end if;
    end if;
    return 0;
  end Job278;

  function Job279 return integer32 is -- retrieve next multprec initialize
  begin
    Multprec_Solutions_Container.Retrieve_Next_Initialize;
    return 0;
  end Job279;

  function Job300 return integer32 is -- set pointer to next standard solution

    ind : natural32;

  begin
    Standard_Solutions_Container.Move_Current(ind);
    Assign(integer32(ind),a);
    return 0;
  exception
    when others => return 300;
  end Job300;

  function Job301 return integer32 is -- set pointer to next dobldobl solution

    ind : natural32;

  begin
    DoblDobl_Solutions_Container.Move_Current(ind);
    Assign(integer32(ind),a);
    return 0;
  exception
    when others => return 301;
  end Job301;

  function Job302 return integer32 is -- set pointer to next quaddobl solution

    ind : natural32;

  begin
    QuadDobl_Solutions_Container.Move_Current(ind);
    Assign(integer32(ind),a);
    return 0;
  exception
    when others => return 302;
  end Job302;

  function Job303 return integer32 is -- set pointer to next multprec solution

    ind : natural32;

  begin
    Multprec_Solutions_Container.Move_Current(ind);
    Assign(integer32(ind),a);
    return 0;
  exception
    when others => return 303;
  end Job303;

  function Job304 return integer32 is -- length of current standard solution

    use Standard_Complex_Solutions;
    ind,len : natural32;
    ls : Link_to_Solution;

  begin
    Standard_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      len := Standard_Solution_Strings.Length(ls.all);
      Assign(integer32(len),b);
    end if;
    return 0;
  exception
    when others => return 304;
  end Job304;

  function Job305 return integer32 is -- length of current dobldobl solution

    use DoblDobl_Complex_Solutions;
    ind,len : natural32;
    ls : Link_to_Solution;

  begin
    DoblDobl_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      len := DoblDobl_Solution_Strings.Length(ls.all);
      Assign(integer32(len),b);
    end if;
    return 0;
  exception
    when others => return 305;
  end Job305;

  function Job306 return integer32 is -- length of current quaddobl solution

    use QuadDobl_Complex_Solutions;
    ind,len : natural32;
    ls : Link_to_Solution;

  begin
    QuadDobl_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      len := QuadDobl_Solution_Strings.Length(ls.all);
      Assign(integer32(len),b);
    end if;
    return 0;
  exception
    when others => return 306;
  end Job306;

  function Job307 return integer32 is -- length of current multprec solution

    use Multprec_Complex_Solutions;
    ind,len : natural32;
    ls : Link_to_Solution;

  begin
    Multprec_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      len := Multprec_Solution_Strings.Length(ls.all);
      Assign(integer32(len),b);
    end if;
    return 0;
  exception
    when others => return 307;
  end Job307;

  function Job308 return integer32 is -- current standard solution string

    use Standard_Complex_Solutions;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    len : constant natural32 := natural32(v_a(v_a'first));
    ind : natural32;
    ls : Link_to_Solution;
    sv : Standard_Integer_Vectors.Vector(1..integer32(len));

  begin
    Standard_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      declare
        s : constant string := Standard_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(len,s));
      end;
      Assign(sv,b);
    end if;
    return 0;
  exception
    when others => return 308;
  end Job308;

  function Job309 return integer32 is -- current dobldobl solution string

    use DoblDobl_Complex_Solutions;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    len : constant natural32 := natural32(v_a(v_a'first));
    ind : natural32;
    ls : Link_to_Solution;
    sv : Standard_Integer_Vectors.Vector(1..integer32(len));

  begin
    DoblDobl_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      declare
        s : constant string := DoblDobl_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(len,s));
      end;
      Assign(sv,b);
    end if;
    return 0;
  exception
    when others => return 309;
  end Job309;

  function Job310 return integer32 is -- current quaddobl solution string

    use QuadDobl_Complex_Solutions;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    len : constant natural32 := natural32(v_a(v_a'first));
    ind : natural32;
    ls : Link_to_Solution;
    sv : Standard_Integer_Vectors.Vector(1..integer32(len));

  begin
    QuadDobl_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      declare
        s : constant string := QuadDobl_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(len,s));
      end;
      Assign(sv,b);
    end if;
    return 0;
  exception
    when others => return 310;
  end Job310;

  function Job311 return integer32 is -- current multprec solution string

    use Multprec_Complex_Solutions;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    len : constant natural32 := natural32(v_a(v_a'first));
    ind : natural32;
    ls : Link_to_Solution;
    sv : Standard_Integer_Vectors.Vector(1..integer32(len));

  begin
    Multprec_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      declare
        s : constant string := Multprec_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(len,s));
      end;
      Assign(sv,b);
    end if;
    return 0;
  exception
    when others => return 311;
  end Job311;

  function Job544 return integer32 is -- read standard sys+sols from file

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);

  begin
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
      sols : Standard_Complex_Solutions.Solution_List;
    begin
      Open(file,in_file,sv);
      Standard_System_and_Solutions_io.get(file,p,sols);
      Standard_PolySys_Container.Clear;
      Standard_PolySys_Container.Initialize(p.all);
      Standard_Solutions_Container.Clear;
      Standard_Solutions_Container.Initialize(sols);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 544;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 544;
    end;
    return 0;
  end Job544;

  function Job545 return integer32 is -- read double double sys+sols from file

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);

  begin
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
      sols : DoblDobl_Complex_Solutions.Solution_List;
    begin
      Open(file,in_file,sv);
      DoblDobl_System_and_Solutions_io.get(file,p,sols);
      DoblDobl_PolySys_Container.Clear;
      DoblDobl_PolySys_Container.Initialize(p.all);
      DoblDobl_Solutions_Container.Clear;
      DoblDobl_Solutions_Container.Initialize(sols);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 545;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 545;
    end;
    return 0;
  end Job545;

  function Job546 return integer32 is -- read quad double sys+sols from file

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);

  begin
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
      sols : QuadDobl_Complex_Solutions.Solution_List;
    begin
      Open(file,in_file,sv);
      QuadDobl_System_and_Solutions_io.get(file,p,sols);
      QuadDobl_PolySys_Container.Clear;
      QuadDobl_PolySys_Container.Initialize(p.all);
      QuadDobl_Solutions_Container.Clear;
      QuadDobl_Solutions_Container.Initialize(sols);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 546;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 546;
    end;
    return 0;
  end Job546;

  function Job547 return integer32 is -- read multiprecision sys+sols from file

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant natural := natural(v_a(v_a'first));
    use Interfaces.C;
    nbdeci : constant natural32 := natural32(v_a(v_a'first+1));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(nbdeci);

  begin
    Multprec_Complex_Polynomials_io.Set_Working_Precision(size);
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
      sols : Multprec_Complex_Solutions.Solution_List;
    begin
      Open(file,in_file,sv);
      Multprec_System_and_Solutions_io.get(file,p,sols);
      Multprec_PolySys_Container.Clear;
      Multprec_PolySys_Container.Initialize(p.all);
      Multprec_Solutions_Container.Clear;
      Multprec_Solutions_Container.Initialize(sols);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 547;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 547;
    end;
    return 0;
  end Job547;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- read from file into container
      when 1 => return Job1; -- write container to file
      when 2 =>
        Assign(integer32(Standard_Solutions_Container.Length),b); return 0;
      when 3 =>
        Assign(integer32(Standard_Solutions_Container.Dimension),b); return 0;
      when 4 => return Job4; -- return standard solution in container
      when 5 => return Job5; -- change standard solution in container
      when 6 => return Job6; -- append standard solution to container
      when 7 => Standard_Solutions_Container.Clear; return 0;
      when 8 => return Job8;   -- drop coordinate by index from solutions
      when 9 => return Job9;   -- drop coordinate by name from solutions
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
      when 30 => return Job30; -- returns size of solution string
      when 31 => return Job31; -- returns solution string
      when 32 => return Job32; -- returns size of solution intro
      when 33 => return Job33; -- returns solution intro string
      when 34 => return Job34; -- returns size of solution vector
      when 35 => return Job35; -- returns solution vector string
      when 36 => return Job36; -- returns size of solution diagnostics
      when 37 => return Job37; -- returns solution diagnostics string
      when 38 => return Job38; -- appends solution string to container
      when 39 => return Job39; -- replaces solution string in container
     -- corresponding operations for double double solutions
      when 40 => return Job40; -- read from file into container
      when 41 => return Job41; -- write container to file
      when 42 =>
        Assign(integer32(DoblDobl_Solutions_Container.Length),b); return 0;
      when 43 =>
        Assign(integer32(DoblDobl_Solutions_Container.Dimension),b); return 0;
      when 44 => return Job44; -- return dobldobl solution in container
      when 45 => return Job45; -- change dobldobl solution in container
      when 46 => return Job46; -- append dobldobl solution to container
      when 47 => DoblDobl_Solutions_Container.Clear; return 0;
      when 48 => return Job48; -- drop coordinate by index from solutions
      when 49 => return Job49; -- drop coordinate by name from solutions
      when 70 => return Job70; -- returns size of solution string
      when 71 => return Job71; -- returns solution string
      when 78 => return Job78; -- appends solution string to container
     -- corresponding operations for quad double solutions
      when 80 => return Job80; -- read from file into container
      when 81 => return Job81; -- write container to file
      when 82 =>
        Assign(integer32(QuadDobl_Solutions_Container.Length),b); return 0;
      when 83 =>
        Assign(integer32(QuadDobl_Solutions_Container.Dimension),b); return 0;
      when 84 => return Job84; -- return quaddobl solution in container
      when 85 => return Job85; -- change quaddobl solution in container
      when 86 => return Job86; -- append quaddobl solution to container
      when 87 => QuadDobl_Solutions_Container.Clear; return 0;
      when 88 => return Job88; -- drop coordinate by index from solutions
      when 89 => return Job89; -- drop coordinate by name from solutions
      when 110 => return Job110; -- returns size of solution string
      when 111 => return Job111; -- returns solution string
      when 118 => return Job118; -- appends solution string to container
     -- corresponding operations for multiprecision solutions
      when 120 => return Job120; -- read from file into container
      when 121 => return Job121; -- write container to file
      when 122 =>
        Assign(integer32(Multprec_Solutions_Container.Length),b); return 0;
      when 123 => 
        Assign(integer32(Multprec_Solutions_Container.Dimension),b); return 0;
      when 127 => Multprec_Solutions_Container.Clear; return 0;
      when 150 => return Job150; -- returns size of solution string
      when 151 => return Job151; -- returns solution string
      when 158 => return Job158; -- append solution string to container
     -- retrieve next solution
      when 276 => return Job276; -- retrieve next standard solution
      when 277 => return Job277; -- retrieve next double double solution
      when 278 => return Job278; -- retrieve next quad double solution
      when 279 => return Job279; -- retrieve next multprec initialize
     -- move pointer from current to next solution
      when 300 => return Job300; -- move to next standard solution
      when 301 => return Job301; -- move to next dobldobl solution
      when 302 => return Job302; -- move to next quaddobl solution
      when 303 => return Job303; -- move to next multprec solution
     -- return length of current solution string
      when 304 => return Job304; -- length of current standard solution string
      when 305 => return Job305; -- length of current dobldobl solution string
      when 306 => return Job306; -- length of current quaddobl solution string
      when 307 => return Job307; -- length of current multprec solution string
     -- return current solution string
      when 308 => return Job308; -- current standard solution string
      when 309 => return Job309; -- current dobldobl solution string
      when 310 => return Job310; -- current quaddobl solution string
      when 311 => return Job311; -- current multprec solution string
     -- reading system and solutions from given file name
      when 544 => return Job544; -- read standard system and solutions
      when 545 => return Job545; -- read double double system and solutions
      when 546 => return Job546; -- read quad double system and solutions
      when 547 => return Job547; -- read multiprecision system and solutions
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_solcon;
