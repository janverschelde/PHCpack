with System_Call,File_Scanning;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Parse_Numbers;
with Characters_and_Numbers;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs_io;        use Standard_Integer_VecVecs_io;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;

package body DEMiCs_Command_Line is

  function Random_Name ( prefix : string ) return string is

    nbr : constant integer32
        := Standard_Random_Numbers.Random(10_000_000,99_999_999);

  begin
    return prefix & Characters_and_Numbers.Convert(nbr);
  end Random_Name;

  procedure Prepare_Input
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                iptname : in string;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

    file : file_type;

  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    if verbose then
      put_line("The supports : "); put(supports);
      put_line("Writing to file " & iptname & " ...");
    end if;
    Create_Output_File(file,iptname);
    put(file,"Dim = "); put(file,supports'last,1); new_line(file);
    put(file,"Support = "); put(file,supports'last,1); new_line(file);
    new_line(file);
    put(file,"Elem =");
    for i in supports'range loop
      put(file," ");
      put(file,Lists_of_Integer_Vectors.Length_Of(supports(i)),1);
    end loop;
    new_line(file);
    put(file,"Type =");
    for i in supports'range loop -- just assume fully mixed ...
      put(file," 1");
    end loop;
    new_line(file);
    new_line(file);
    put(file,supports);
    close(file);
  end Prepare_Input;

  procedure Call_DEMiCs
              ( infile,outfile : in string;
                execfile : in string := "/tmp/demics";
                verbose : in boolean := true ) is

    cmd : constant string
        := execfile & " -c " & infile & " > " & outfile;

  begin
    if verbose then
      put_line("Calling " & cmd);
    end if;
    System_Call.Call(cmd);
  exception
    when System_Call.SYSTEM_ERROR =>
      put_line("Executable not found at " & execfile & "!?");
      raise;
  end Call_DEMiCS;

  function Extract_Lifting_Values
             ( vals : string ) return Standard_Floating_Vectors.Vector is

    dlm : constant character := ' ';
    cnt : constant natural := String_Splitters.Count_Delimiters(vals,dlm);
    nbr : Array_of_Strings(1..integer(cnt));
    res : Standard_Floating_Vectors.Vector(1..integer32(cnt));
    pos : integer;

  begin
    nbr := String_Splitters.Split(cnt,vals,dlm);
    res := (res'range => 0.0);
    for i in nbr'range loop
      pos := nbr(i)'first;
      Standard_Parse_Numbers.Parse(nbr(i).all,pos,res(integer32(i)));
    end loop;
    return res;
  end Extract_Lifting_Values;

  procedure Parse_Lifting
              ( file : in file_type; dim : in integer32;
                lifting : out Standard_Floating_VecVecs.VecVec;
                verbose : in boolean := true ) is

    found : boolean;
    lifting_banner : constant string := "Lifting values";
    support_banner : constant string := "S";
    cnt : natural32 := 0;

  begin
    File_Scanning.Scan_and_Skip(file,lifting_banner,found);
    if not found then
      put_line("Error: no lifting values found!?");
    else
      if verbose then
        put_line("Found the lifting values banner.");
      end if;
      for i in 1..dim loop
        File_Scanning.Scan(file,support_banner,found);
        exit when not found;
        cnt := cnt + 1;
        File_Scanning.Scan(file,": ",found);
        declare
          strlifvals : constant string := text_io.get_line(file);
          lifvals : constant Standard_Floating_Vectors.Vector
                  := Extract_Lifting_Values(strlifvals);
        begin
          if verbose then
            put("The lifting for support "); put(cnt,1);
            put_line(" : " & strlifvals);
            put_line(lifvals);
          end if;
          lifting(integer32(i))
            := new Standard_Floating_Vectors.Vector'(lifvals);
        end;
      end loop;
      if verbose then
        put("Number of supports : "); put(cnt,1); new_line;
      end if;
    end if;
  end Parse_Lifting;

  function Extract_Cell_Indices
              ( dim : integer32; vals : string;
                verbose : boolean := true )
              return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..2*dim); -- fully mixed
    pos : integer := vals'first;
    idx : integer32;
    nbr : natural32;
    sign : character;
    idxfirst,idxsecond : integer32;

  begin
    for i in 1..dim loop
      while pos <= vals'last loop         -- skip spaces
        exit when (vals(pos) /= ' ');
        pos := pos + 1;
      end loop;
      Standard_Parse_Numbers.Parse(vals,pos,idx,nbr,sign);
      if verbose then
        put("Indices for component "); put(idx,1);
        put(" : ");
      end if;
      while pos <= vals'last loop          -- read till (
        exit when (vals(pos) = '(');
        pos := pos + 1;
      end loop;
      pos := pos + 1; -- skip the (
      while pos <= vals'last loop          -- skip spaces
        exit when (vals(pos) /= ' ');
        pos := pos + 1;
      end loop;
      Standard_Parse_Numbers.Parse(vals,pos,idxfirst,nbr,sign);
      while pos <= vals'last loop          -- skip spaces
        exit when (vals(pos) /= ' ');
        pos := pos + 1;
      end loop;
      Standard_Parse_Numbers.Parse(vals,pos,idxsecond,nbr,sign);
      while pos <= vals'last loop          -- read till )
        exit when (vals(pos) = ')');
        pos := pos + 1;
      end loop;
      pos := pos + 1; -- skip the )
      if verbose then
        put(" "); put(idxfirst,1);
        put(" "); put(idxsecond,1);
        new_line;
      end if;
      res(2*(idx-1) + 1) := idxfirst;
      res(2*(idx-1) + 2) := idxsecond;
    end loop;
    return res;
  end Extract_Cell_Indices;

  procedure Parse_Cells
              ( file : in file_type; dim : in integer32;
                cells : out Lists_of_Integer_Vectors.List;
                verbose : in boolean := true ) is

    found : boolean;
    ch : character;
    cnt : natural32 := 0;
    cells_last : Lists_of_Integer_Vectors.List := cells;

  begin
    loop
      File_Scanning.Scan(file,"# ",found);
      exit when not found;
      get(file,ch);
      exit when (ch = 'M'); -- final count of mixed cells
      cnt := cnt + 1;
      File_Scanning.Scan(file,":",found);
      declare
        strcell : constant string := text_io.get_line(file);
        idxcell : Standard_Integer_Vectors.Vector(1..2*dim);
      begin
        if verbose then
          put("Cell "); put(cnt,1);
          put_line(" :" & strcell);
        end if;
        idxcell := Extract_Cell_Indices(dim,strcell,verbose);
        if verbose then
          put("The indices of the mixed cell :");
          put(idxcell); new_line;
        end if;
        Lists_of_Integer_Vectors.Append(cells,cells_last,idxcell);
      end;
    end loop;
    if verbose
     then put("Number of mixed cells : "); put(cnt,1); put_line(".");
    end if;
  end Parse_Cells;

  procedure Process_Output
              ( dim : in integer32; filename : in string;
                mv : out natural32;
                lif : out Standard_Floating_VecVecs.VecVec;
                cells : out Lists_of_Integer_Vectors.List;
                verbose : in boolean := true ) is

    file : file_type;
    banner : constant string := "Mixed Volume:";
    found : boolean;
    ch,sign : character;
    nb : natural32;
    imv : integer32 := 0;

  begin
    if verbose then
      put_line("Opening " & filename & " ...");
    end if;
    Open_Input_File(file,filename);
    Parse_Lifting(file,dim,lif,verbose);
    Parse_Cells(file,dim,cells,verbose);
    if verbose
     then put_line("The lifting of the supports :"); put(lif);
    end if;
    File_Scanning.Scan(file,banner,found);
    if not found then
      put_line("No mixed volume computed?");
    else
      get(file,ch);
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,ch);
      Standard_Parse_Numbers.Parse(file,ch,imv,nb,sign);
      if verbose then
        new_line;
        put("The mixed volume : "); put(imv,1); new_line;
      end if;
      mv := natural32(imv);
    end if;
  end Process_Output;

end DEMiCs_Command_Line;
