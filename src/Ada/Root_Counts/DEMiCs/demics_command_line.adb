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
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;

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
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

    file : file_type;
    idx : integer32;
  
  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
      put_line("Writing to file " & iptname & " ...");
    end if;
    Create_Output_File(file,iptname);
    put(file,"Dim = "); -- the ambient dimension
    put(file,supports'last,1); new_line(file);
    put(file,"Support = "); -- number of different supports
    put(file,mix'last,1); new_line(file);
    new_line(file);
    put(file,"Elem ="); -- length of each support list
    idx := supports'first;
    for i in mix'range loop
      put(file," ");
      put(file,Lists_of_Integer_Vectors.Length_Of(supports(idx)),1);
      idx := idx + mix(i);
    end loop;
    new_line(file);
    put(file,"Type =");
    for i in mix'range loop
      put(file," "); put(file,mix(i),1);
    end loop;
    new_line(file);
    new_line(file);
    idx := supports'first;
    for i in mix'range loop
      put(file,supports(idx));
      idx := idx + mix(i);
    end loop;
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
              ( file : in file_type;
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
      for i in lifting'range loop
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

  function Offset_for_Index
              ( mix : Standard_Integer_Vectors.Vector;
                idx : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the offset for the index of the component idx,
  --   relative to the type of mixture in mix.

    res : integer32 := 0;

  begin
    for i in mix'first..(idx-1) loop
      res := res + mix(i)+1;
    end loop;
    return res;
  end Offset_for_Index;

  function Offset_for_Index
              ( mix : Standard_Integer_Vectors.Link_to_Vector;
                idx : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the offset for the index of the component idx,
  --   relative to the type of mixture in mix.

    res : integer32 := 0;

  begin
    for i in mix'first..(idx-1) loop
      res := res + mix(i)+1;
    end loop;
    return res;
  end Offset_for_Index;

  function Extract_Cell_Indices
              ( nbrpts : integer32;
                mix : Standard_Integer_Vectors.Vector;
                vals : string; verbose : boolean := true )
              return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..nbrpts);
    residx : integer32;
    pos : integer := vals'first;
    idx,idx2pnt : integer32;
    nbr : natural32;
    sign : character;

  begin
    for i in mix'range loop
      while pos <= vals'last loop         -- skip spaces
        exit when (vals(pos) /= ' ');
        pos := pos + 1;
      end loop;
      Standard_Parse_Numbers.Parse(vals,pos,idx,nbr,sign);
      residx := Offset_for_Index(mix,idx);
      if verbose then
        put("Indices for component "); put(idx,1);
        put(" : ");
      end if;
      while pos <= vals'last loop          -- read till (
        exit when (vals(pos) = '(');
        pos := pos + 1;
      end loop;
      pos := pos + 1; -- skip the (
      for k in 1..mix(idx)+1 loop          -- by greedy visit, idx /= i
        while pos <= vals'last loop        -- skip spaces
          exit when (vals(pos) /= ' ');
          pos := pos + 1;
        end loop;
        Standard_Parse_Numbers.Parse(vals,pos,idx2pnt,nbr,sign);
        if verbose then
          put(" "); put(idx2pnt);
        end if;
        residx := residx + 1;
        res(residx) := idx2pnt;
      end loop;
      if verbose
       then new_line;
      end if;
      while pos <= vals'last loop          -- read till )
        exit when (vals(pos) = ')');
        pos := pos + 1;
      end loop;
      pos := pos + 1; -- skip the )
    end loop;
    return res;
  end Extract_Cell_Indices;

  procedure Line2Cell_Indices
              ( line : in string; nbrpts : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                indices : out Standard_Integer_Vectors.Vector;
                verbose : boolean := true ) is

    residx : integer32;
    pos : integer := line'first;
    idx,idx2pnt : integer32;
    nbr : natural32;
    sign : character;

  begin
    if verbose then
      put("Inside Line2Cell_Indices, nbrpts = "); put(nbrpts,1); new_line;
      put("mix = "); put(mix.all); new_line;
      put("line = "); put_line(line);
    end if;
    while line(pos) /= ':' loop -- skip the number of the cell
      pos := pos + 1;
    end loop;
    pos := pos + 1; -- must skip the ':'
    for i in mix'range loop
      while pos <= line'last loop         -- skip spaces
        exit when (line(pos) /= ' ');
        pos := pos + 1;
      end loop;
      Standard_Parse_Numbers.Parse(line,pos,idx,nbr,sign);
      residx := Offset_for_Index(mix,idx);
      if verbose then
        put("Indices for component "); put(idx,1);
        put(" : ");
      end if;
      while pos <= line'last loop          -- read till (
        exit when (line(pos) = '(');
        pos := pos + 1;
      end loop;
      pos := pos + 1; -- skip the (
      for k in 1..mix(idx)+1 loop          -- by greedy visit, idx /= i
        while pos <= line'last loop        -- skip spaces
          exit when (line(pos) /= ' ');
          pos := pos + 1;
        end loop;
        Standard_Parse_Numbers.Parse(line,pos,idx2pnt,nbr,sign);
        if verbose then
          put(" "); put(idx2pnt);
        end if;
        residx := residx + 1;
        indices(residx) := idx2pnt;
      end loop;
      if verbose
       then new_line;
      end if;
      while pos <= line'last loop          -- read till )
        exit when (line(pos) = ')');
        pos := pos + 1;
      end loop;
      pos := pos + 1; -- skip the )
    end loop;
  end Line2Cell_Indices;

  function Number_of_Points_in_Cell
              ( mix : Standard_Integer_Vectors.Vector ) 
              return integer32 is

    res : integer32 := 0;

  begin
    for i in mix'range loop
      res := res + mix(i) + 1;
    end loop;
    return res;
  end Number_of_Points_in_Cell;

  procedure Parse_Cells
              ( file : in file_type; dim : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                cells : out Lists_of_Integer_Vectors.List;
                verbose : in boolean := true ) is

    found : boolean;
    ch : character;
    cnt : natural32 := 0;
    cells_last : Lists_of_Integer_Vectors.List := cells;
    nbr : constant integer32 := Number_of_Points_in_Cell(mix);

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
        idxcell : Standard_Integer_Vectors.Vector(1..nbr);
      begin
        if verbose then
          put("Cell "); put(cnt,1);
          put_line(" :" & strcell);
        end if;
        idxcell := Extract_Cell_Indices(nbr,mix,strcell,verbose);
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
                mix : in Standard_Integer_Vectors.Vector;
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
    Parse_Lifting(file,lif,verbose);
    Parse_Cells(file,dim,mix,cells,verbose);
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
