package body Checker_Boards is

  procedure Initialize ( b : in out Board ) is
  begin
    for i in b'range(1) loop
      for j in b'range(2) loop
        b(i,j) := ' ';
      end loop;
    end loop;
  end Initialize;

  function Configuration ( p : Vector ) return Board is

    res : Board(p'range,p'range);
 
  begin
    Initialize(res);
    for i in p'range loop
      res(integer32(p(i)),p'last+1-i) := '*';
    end loop;
    return res;
  end Configuration;

  function Column ( p : Vector; r : natural32 ) return integer32 is
  begin
    for i in p'range loop
      if p(i) = r
       then return p'last+1-i;
      end if;
    end loop;
    return 0;
  end Column;

  function Permutation ( b : Board ) return Vector is

    res : Vector(b'range(1));

  begin
    for j in reverse b'range(2) loop
      for i in b'range(1) loop
        if b(i,j) /= ' '
         then res(b'last(1)+1-j) := natural32(i);
        end if;
      end loop;
    end loop;
    return res;
  end Permutation;

  procedure Place_White ( b : in out Board; r,c : in Vector ) is

    col : integer32;

  begin
    for i in r'range loop
      col := integer32(c(c'last-i+1));
      if b(integer32(r(i)),col) = ' '
       then b(integer32(r(i)),col) := 'o';
       else b(integer32(r(i)),col) := '%';
      end if;
    end loop;
  end Place_White;

  function Configuration ( p,r,c : Vector ) return Board is

    res : Board(p'range,p'range) := Configuration(p);

  begin
    Place_White(res,r,c);
    return res;
  end Configuration;

end Checker_Boards;
