with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;

package body Standard_Solutions_Container is

-- INTERNAL DATA :

  first,last : Solution_List;
  input_file,output_file,wfile1,wfile2 : file_type;
  ls1,ls2,solsym : Symbol_Table.Link_to_Array_of_Symbols;

-- OPERATIONS :

  procedure Initialize ( sols : in Solution_List ) is

    tmp : Solution_List := sols;
 
  begin
    for i in 1..Length_Of(sols) loop
      Append(first,last,Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize;

  function Length return natural32 is
  begin
    return Length_Of(first);
  end Length;

  function Dimension return natural32 is
  begin
    if Is_Null(first)
     then return 0;
     else return natural32(Head_Of(first).n);
    end if;
  end Dimension;

  function Retrieve return Solution_List is
  begin
    return first;
  end Retrieve;

  procedure Retrieve ( k : in natural32; s : out Solution;
                       fail : out boolean ) is

    ls : Link_to_Solution;

  begin
    Retrieve(k,ls,fail);
    if not fail
     then s := ls.all;
    end if;
  end Retrieve;

  procedure Retrieve ( k : in natural32; s : out Link_to_Solution;
                       fail : out boolean ) is

    tmp : Solution_List := first;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k then
        fail := false;
        s := Head_Of(tmp);
        return;
      else
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    fail := true;
  end Retrieve;

  procedure Replace ( k : in natural32; s : in Solution;
                      fail : out boolean ) is
	  
    tmp : Solution_List := first;
    ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k then
        fail := false;
        ls := Head_Of(tmp);
        ls.t := s.t;
        ls.m := s.m;
        ls.v := s.v;
        ls.err := s.err;
        ls.rco := s.rco;
        ls.res := s.res;
        Set_Head(tmp,ls);
        return;
      else
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    fail := true;
  end Replace;

  procedure Replace ( k : in natural32; s : in Link_to_Solution;
                      fail : out boolean ) is
	  
    tmp : Solution_List := first;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k then
        fail := false;
        Set_Head(tmp,s);
        return;
      else
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    fail := true;
  end Replace;

  procedure Append ( s : in Solution ) is
  begin
    Append(first,last,s);
  end Append;

  procedure Append ( s : in Link_to_Solution ) is
  begin
    Append(first,last,s);
  end Append;

  procedure Clear is
  begin
    Clear(first);
    last := first;
  end Clear;

-- handing files for incremental read and write of solution lists :

  procedure Silent_Open_Input_File ( filename : in string ) is
  begin
    Open(input_file,in_file,filename);
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File
              ( k : in natural32; filename : in string ) is
  begin
    if k = 1 then
      Open(wfile1,in_file,filename);
    elsif k = 2 then
      Open(wfile2,in_file,filename);
    end if;
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File is
  begin
    Read_Name_and_Open_File(input_file);
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File ( k : in natural32 ) is
  begin
    if k = 1 then
      Read_Name_and_Open_File(wfile1);
    elsif k = 2 then
      Read_Name_and_Open_File(wfile2);
    end if;
  end Silent_Open_Input_File;

  procedure Open_Input_File is
  begin
    put_line("Reading the name of the input file for solutions.");
    Read_Name_and_Open_File(input_file);
  end Open_Input_File;

  procedure Open_Input_File ( k : in natural32 ) is
  begin
    put("Reading the name of the input file for witness set ");
    put(k,1); put_line(".");
    if k = 1 then
      Read_Name_and_Open_File(wfile1);
    elsif k = 2 then
      Read_Name_and_Open_File(wfile2);
    end if;
  end Open_Input_File;

  procedure Create_Output_File is
  begin
    put_line("Reading the name of the output file for solutions.");
    Read_Name_and_Create_File(output_file);
  end Create_Output_File;

  function Solution_Input_File return file_type is
  begin
    return input_file;
  end Solution_Input_File;

  function Solution_Input_File ( k : natural32 ) return file_type is
  begin
    if k = 1 then
      return wfile1;
    elsif k = 2 then
      return wfile2;
    else
      return input_file;
    end if;
  end Solution_Input_File;

  function Solution_Output_File return file_type is
  begin
    return output_file;
  end Solution_Output_File;

  procedure Reset_Input_File ( k : in natural32 ) is
  begin
    if k = 1 then
      Reset(wfile1);
    elsif k = 2 then
      Reset(wfile2);
    end if;
  end Reset_Input_File;

  procedure Close_Input_File is
  begin
    close(input_file);
  end Close_Input_File;

  procedure Close_Input_File ( k : in natural32 ) is
  begin
    if k = 0 then
      close(input_file);
    elsif k = 1 then
      close(wfile1);
    elsif k = 2 then 
      close(wfile2);
    end if;
  end Close_Input_File;

  procedure Close_Output_File is
  begin
    close(output_file);
  end Close_Output_File;

-- Management of two symbol tables for diagonal homotopies :

  procedure Store_Symbol_Table
              ( k : in natural32; sbs : Symbol_Table.Array_of_Symbols ) is
  begin
    if k = 0 then
      solsym := new Symbol_Table.Array_of_Symbols'(sbs);
    elsif k = 1 then
      ls1 := new Symbol_Table.Array_of_Symbols'(sbs);
    elsif k = 2 then
      ls2 := new Symbol_Table.Array_of_Symbols'(sbs);
    end if;
  end Store_Symbol_Table;

  function Retrieve_Symbol_Table 
              ( k : in natural32 )
              return Symbol_Table.Link_to_Array_of_Symbols is

    res : Symbol_Table.Link_to_Array_of_Symbols := null;

  begin
    if k = 0 then
      return solsym;
    elsif k = 1 then
      return ls1;
    elsif k = 2 then
      return ls2;
    else 
      return res;
    end if;
  end Retrieve_Symbol_Table;

  procedure Clear_Symbol_Table ( k : in natural32 ) is
  begin
    if k = 0 then
      Symbol_Table.Clear(solsym);
    elsif k = 1 then
      Symbol_Table.Clear(ls1);
    elsif k = 2 then
      Symbol_Table.Clear(ls2);
    end if;
  end Clear_Symbol_Table;

end Standard_Solutions_Container;
