package body File_Scanning is

  function Scan_Line_for_Number ( file : in file_type ) return natural is

    m : natural := 0;
    ch : character;

  begin
    while not END_OF_LINE(file) loop
      get(file,ch);
      case ch is
        when '0' => m := 10*m;
        when '1' => m := 10*m + 1;
        when '2' => m := 10*m + 2;
        when '3' => m := 10*m + 3;
        when '4' => m := 10*m + 4;
        when '5' => m := 10*m + 5;
        when '6' => m := 10*m + 6;
        when '7' => m := 10*m + 7;
        when '8' => m := 10*m + 8;
        when '9' => m := 10*m + 9;
        when others => null;
      end case;
    end loop;
    return m;
  end Scan_Line_for_Number;

  procedure Scan ( file : in file_type; ch : in character;
                   found : out boolean ) is

    c : character;

  begin
    while not End_of_File(file) loop
      get(file,c);
      if c = ch
       then found := true;
            return;
      end if;
    end loop;
    found := false;
  end Scan;

  procedure Scan ( file : in file_type; banner : in string;
                   found : out boolean ) is

    index : natural := banner'first-1;
    ch : character;

  begin
    while not End_of_File(file) loop
      get(file,ch);
      if index < banner'first
       then if ch = banner(banner'first)
             then index := banner'first+1;
            end if;
       else if ch = banner(index)
             then index := index + 1;
             else index := banner'first-1;
            end if;
      end if;
      exit when index > banner'last;
    end loop;
    if index > banner'last
     then found := true;
     else found := false;
    end if;
  exception
    when others => found := false; return;
  end Scan;

  procedure Scan_Line ( file : in file_type; banner : in string;
                        found : out boolean ) is

    index : natural := banner'first-1;
    ch : character;

  begin
    while not End_of_Line(file) loop
      get(file,ch);
      if index < banner'first
       then if ch = banner(banner'first)
             then index := banner'first+1;
            end if;
       else if ch = banner(index)
             then index := index + 1;
             else index := banner'first-1;
            end if;
      end if;
      exit when index > banner'last;
    end loop;
    if index > banner'last
     then found := true;
     else found := false;
    end if;
  exception
    when others => found := false; return;
  end Scan_Line;

  procedure Scan_and_Skip ( file : in file_type; banner : in string;
                            found : out boolean ) is

    fnd : boolean;

  begin
    Scan(file,banner,fnd);
    if fnd
     then skip_line(file);
    end if;
    found := fnd;
  end Scan_and_Skip;

end File_Scanning;
