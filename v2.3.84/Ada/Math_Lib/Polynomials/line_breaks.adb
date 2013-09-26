package body Line_Breaks is

  column : natural32 := 0;

  procedure Init_Line is
  begin
    column := 0;
  end Init_Line;

  procedure Line ( file : in file_type; n : natural32 ) is
  begin
    if column + n >= right then
      new_line(file);
      column := 0;
    else
      column := column +  n;
    end if;
  end Line;

  function Length ( d,i : integer32; std : boolean; pw : power )
                  return natural32 is

    l : natural32 := 0;
    sb : symbol;

  begin
    if std then
      if i < 10
       then l := l + 2;
       else l := l + 3;
      end if;
    else
      sb := Symbol_Table.get(natural32(i));
      l := l + Length(sb);
    end if;
    if d > 1 then
      if pw = '^'
       then l := l + 1;
       else l := l + 2;
      end if;
      if d < 10
       then l := l + 1;
       else l := l + 2;
      end if;
    end if;
    return l;
  end Length;

end Line_Breaks;
