with File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with Standard_Select_Solutions;

package body DoblDobl_Select_Solutions is

  procedure Scan_Solutions
              ( file : in file_type; len,dim : in natural32;
                sv : in Vector; sa : out Solution_Array ) is

    s : Solution(integer32(dim));
    freq : natural32 := 1024; -- frequency updater
    ind : integer32 := sv'first;

  begin
    put("Scanning solutions ... ");
    for i in 1..len loop
      Read_Next(file,s);
      if i mod freq = 0
       then put(i,1); put(" ... "); freq := 2*freq;
      end if;
      if i = sv(ind)
       then sa(ind) := new Solution'(s); ind := ind + 1;
      end if;
      exit when (ind > sv'last);
    end loop;
    new_line;
  exception when others =>
    put_line("An exception was raised before reached end...");
    raise;
  end Scan_Solutions;

  procedure Scan_Solutions
              ( file : in file_type; len,dim : in natural32;
                sv : in Vector; sols : out Solution_List ) is

    last : Solution_List := sols;
    s : Solution(integer32(dim));
    freq : natural32 := 1024; -- frequency updater
    ind : integer32 := sv'first;

  begin
    put("Scanning solutions ... ");
    for i in 1..len loop
      Read_Next(file,s);
      if i mod freq = 0
       then put(i,1); put(" ... "); freq := 2*freq;
      end if;
      if i = sv(ind)
       then Append(sols,last,s); ind := ind + 1;
      end if;
      exit when (ind > sv'last);
    end loop;
    new_line;
  exception when others =>
    put_line("An exception was raised before reached end...");
    raise;
  end Scan_Solutions;

  function Sort ( v : Vector ) return Vector is

    res : Vector(v'range) := v;
    tmp : natural32;

  begin
    for i in v'first..v'last-1 loop
      for j in i+1..v'last loop
        if res(j) < res(i) then
          tmp := res(i);
          res(i) := res(j);
          res(j) := tmp;
        end if;
      end loop;
    end loop;
    return res;
  end Sort;

  function Find ( a : natural32; v : Vector ) return integer32 is
  begin
    for i in v'range loop
      if v(i) = a
       then return i;
      end if;
    end loop;
    return 0;
  end Find;

  procedure Write_Selection
               ( file : in file_type; dim : in natural32;
                 rv,sv : in Vector; sa : in Solution_Array ) is

    cnt : natural32 := 0;
    ind : integer32;

  begin
    Write_First(file,natural32(rv'last),dim);
    for i in sa'range loop
      cnt := rv(i) - 1;
      ind := Find(rv(i),sv);
      if ind > 0
       then Write_Next(file,cnt,sa(ind));
      end if;
    end loop;
  end Write_Selection;

  procedure Store_Selection
               ( sols : out Solution_List;
                 rv,sv : in Vector; sa : in Solution_Array ) is

    last : Solution_List := sols;
    ind : integer32;

  begin
    for i in sa'range loop
      ind := Find(rv(i),sv);
      if ind > 0
       then Append(sols,last,sa(ind).all);
      end if;
    end loop;
  end Store_Selection;

  procedure Select_Solutions
              ( sols : in Solution_List;
                rv : in Standard_Natural_Vectors.Vector;
                sv : out Standard_Natural_Vectors.Vector;
                sa : out Solution_Array ) is

    tmp : Solution_List := sols;
    ind : integer32 := rv'first;

  begin
    sv := Sort(rv);
    for i in 1..Length_Of(sols) loop
      if i = sv(ind) then
        sa(ind) := new Solution'(Head_Of(tmp).all);
        ind := ind + 1;
      end if;
      exit when (ind > sv'last);
      tmp := Tail_Of(tmp);
    end loop;
  end Select_Solutions;

  procedure Select_from_File
              ( file : in file_type;
                rv : in Standard_Natural_Vectors.Vector;
                sa : out Solution_Array; fail : out boolean ) is

   length,dim : natural32;
   sv : Standard_Natural_Vectors.Vector(rv'range);
   bannered : boolean;

  begin
    Standard_Select_Solutions.Scan_Banner_Dimensions
      (file,length,dim,bannered,fail);
    if fail then
      put_line("failed to scan the solution file again...");
    else
      put("ready to scan "); put(length,1);
      put(" solutions of dimension "); put(dim,1); put_line(" ...");
      sv := Sort(rv);
      Scan_Solutions(file,length,dim,sv,sa);
    end if;
  end Select_from_File;

  procedure Select_from_File
              ( file : in file_type;
                rv : in Standard_Natural_Vectors.Vector;
                sv : out Standard_Natural_Vectors.Vector;
                sa : out Solution_Array; fail : out boolean ) is

   length,dim : natural32;
   bannered : boolean;

  begin
    Standard_Select_Solutions.Scan_Banner_Dimensions
      (file,length,dim,bannered,fail);
    if fail then
      put_line("failed to scan the solution file again...");
    else
      put("ready to scan "); put(length,1);
      put(" solutions of dimension "); put(dim,1); put_line(" ...");
      sv := Sort(rv);
      Scan_Solutions(file,length,dim,sv,sa);
    end if;
  end Select_from_File;

  procedure Select_from_File
              ( file : in file_type;
                rv : in Standard_Natural_Vectors.Vector;
                sols : out Solution_List ) is

   fail : boolean;
   sv : Standard_Natural_Vectors.Vector(rv'range);
   sa : Solution_Array(rv'range);

  begin
    Select_from_File(file,rv,sv,sa,fail);
    if not fail then
      Store_Selection(sols,rv,sv,sa);
      Clear(sa);
    end if;
  end Select_from_File;

  procedure Select_from_File
              ( file : in file_type; bannered : in boolean;
                rv : in Standard_Natural_Vectors.Vector;
                sv : out Standard_Natural_Vectors.Vector;
                sa : out Solution_Array; fail : out boolean ) is

   length,dim : natural32;
   found : boolean;

  begin
    if bannered then
      File_Scanning.Scan_and_Skip(file,"THE SOLUTIONS",found);
      fail := not found;
    end if;
    if not fail
     then Standard_Select_Solutions.Read_Dimensions(file,length,dim,fail);
    end if;
    if fail then
      put_line("failed to scan the solution file again...");
    else
      put("ready to scan "); put(length,1);
      put(" solutions of dimension "); put(dim,1); put_line(" ...");
      sv := Standard_Select_Solutions.Sort(rv);
      Scan_Solutions(file,length,dim,sv,sa);
    end if;
  end Select_from_File;

  procedure Select_from_File
              ( file : in file_type; bannered : in boolean;
                rv : in Standard_Natural_Vectors.Vector;
                sols : out Solution_List ) is

   fail : boolean;
   sv : Standard_Natural_Vectors.Vector(rv'range);
   sa : Solution_Array(sv'range);

  begin
    Select_from_File(file,bannered,rv,sv,sa,fail);
    if not fail then
      Store_Selection(sols,rv,sv,sa);
      Clear(sa);
    end if;
  end Select_from_File;

end DoblDobl_Select_Solutions;
