with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;

package body Filtered_Points is

  function Create ( n,L,d : integer32 ) return Vector is

    res : Vector(0..n) := (0..n => 0);

  begin
    res(L) := d;
    res(n) := L;
    return res;
  end Create;

  procedure Update ( fp : in out List; d,i,L : in integer32 ) is

    tmp : List := fp;
    lpt : Link_to_Vector;
    cnt : integer32 := 0;
    found : boolean := false;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if lpt(lpt'last) = d then
        cnt := cnt + 1;
        if cnt = i then
          lpt(L) := lpt(L) + 1;
          Set_Head(tmp,lpt);
          found := true;
        end if;
      end if;
      exit when found;
      tmp := Tail_Of(tmp);
    end loop;
  end Update;

  procedure Write ( file : in file_type; fp : in List; i : in integer32 ) is

  -- DESCRIPTION :
  --   Writes the i-th component of every element in the list.

    tmp : List := fp;
    lpt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      put(file,lpt(i),3);
      tmp := Tail_Of(tmp);
    end loop;
    new_line(file);
  end Write;

  procedure Write ( file : in file_type; fp : in List ) is

    n : integer32;

  begin
    if not Is_Null(fp) then
      n := Head_Of(fp)'last;
      put(file,"dimension "); Write(file,fp,n);
      if n > 0 then
        put(file,"level "); put(file,n-1,1);
        put(file,"   "); Write(file,fp,n-1);
        for i in reverse 0..n-2 loop
          put(file,"      "); put(file,i,1);
          put(file,"   "); Write(file,fp,i);
        end loop;
      end if;
    end if;
  end Write;

end Filtered_Points;
