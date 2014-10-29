with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;

package body Checker_Moves is

-- PART I : Moves of the Black Checkers

  function Identity_Permutation ( n : natural32 ) return Vector is

    res : Vector(1..integer32(n));
 
  begin
    for i in 1..n loop
      res(integer32(i)) := i;
    end loop;
    return res;
  end Identity_Permutation;

  function Reverse_Permutation ( n : natural32 ) return Vector is

    res : Vector(1..integer32(n));
 
  begin
    for i in 1..n loop
      res(integer32(i)) := n-i+1;
    end loop;
    return res;
  end Reverse_Permutation;

  function Descending_Checker ( p : Vector ) return integer32 is
  begin
    for i in reverse p'first..p'last-1 loop
      if p(i+1) > p(i)
       then return i;
      end if;
    end loop;
    return 0;
  end Descending_Checker;

  function Falling_Checker ( p : Vector ) return integer32 is
  begin
    for i in p'range loop 
      if integer32(p(i)) > i
       then return i;
      end if;
    end loop;
    return 0;
  end Falling_Checker;

  function Rising_Checker ( p : Vector; i : integer32 ) return integer32 is
  begin
    if i /= 0 then
      for j in i+1..p'last loop
        if p(j) = p(i) + 1
         then return j;
        end if;
      end loop;
    end if;
    return 0;
  end Rising_Checker;

  function Ascending_Checker ( p : Vector; i : integer32 ) return integer32 is
  begin
    if i /= 0 then
      for j in i+1..p'last loop
        if p(j) = p(i) - 1
         then return j;
        end if;
      end loop;
    end if;
    return 0;
  end Ascending_Checker;

  procedure Move ( p : in out Vector; down,up : in integer32 ) is
  begin
    p(down) := p(down) + 1;
    p(up) := p(up) - 1;
  end Move;

  function Number_of_Moves ( n : natural32 ) return natural32 is

    res : constant natural32 := (n*n - n)/2 + 1;

  begin
    return res;
  end Number_of_Moves;

-- PART II : Moves of the White Checkers

  function Critical_Row 
              ( r,c : integer32; rows,cols : Vector ) return integer32 is
  begin
    for i in rows'range loop
      if (integer32(rows(i)) = r) and (integer32(cols(cols'last-i+1)) >= c)
       then return i;
      end if; 
    end loop;
    return 0;
  end Critical_Row;

  function In_Critical_Diagonal 
              ( dr,dc,n,row,column : integer32 ) return boolean is

    k : integer32 := 0;

  begin
    while ((dr + k <= n) and (dc + k <= n)) loop
      if ((dr + k = row) and (dc + k = column))
       then return true;
       else k := k+1;
      end if;
    end loop;
    return false;
  end In_Critical_Diagonal;

  function Top_White_Checker
              ( dr,dc,n : integer32; rows,cols : Vector ) return integer32 is

    res,c : integer32;

  begin
    res := 0;
    for i in rows'range loop
      c := integer32(cols(cols'last-i+1));
      if In_Critical_Diagonal(dr,dc,n,integer32(rows(i)),c) then
        if res = 0 then
          res := i;
        elsif rows(i) < rows(res) then
          res := i;
        end if;
      end if;
    end loop;
    return res;
  end Top_White_Checker;

  function Between ( a,b,x : natural32 ) return boolean is
  begin
    if a <= b
     then return (x >= a) and (x <= b);
     else return (x >= b) and (x <= a);
    end if;
  end Between;

  function Blocker ( rows,cols : Vector; t,r : integer32 ) return boolean is
  begin
    for i in rows'range loop
       if i /= r and i /= t then
         if Between(rows(t),rows(r),rows(i)) 
              and Between(cols(cols'last-t+1),
                          cols(cols'last-r+1),cols(cols'last-i+1))
          then return true;
         end if;
      end if;
    end loop;
    return false;
  end Blocker;

  function Central_Choice
              ( file : in file_type; p : Vector; d,r : integer32;
                rows,cols : Vector; cr,cd : integer32; debug : boolean ) 
              return integer32 is
  begin
    if debug then
      put(file,"Descending black checker at (");
      put(file,p(d),1); put(file,","); put(file,p'last-d+1,1);
      put_line(file,").");
      put(file,"White checker in critical row at (");
      put(file,rows(cr),1); put(file,",");
      put(file,cols(cols'last-cr+1),1); put_line(file,").");
    end if;
   -- is white checker in critical row in descending checker's square?
    if p(d) = rows(cr) and p'last-d+1 = integer32(cols(cols'last-cr+1)) then
      if debug
       then put_line(file,"No choice because in descending checker's square.");
      end if;
      return 1;   -- we must swap
    else
      if debug 
       then put_line(file,"White checker not in descending checker's square.");
      end if;
    end if;
    if debug then
      put(file,"Rising black checker at (");
      put(file,p(r),1); put(file,","); put(file,p'last-r+1,1);
      put_line(file,").");
      put(file,"Top white checker in critical diagonal at (");
      put(file,rows(cd),1); put(file,",");
      put(file,cols(cols'last-cd+1),1); put_line(file,").");
    end if;
   -- is top white checker in critical row in rising checker's square?
    if p(r) = rows(cd) and p'last-r+1 = integer32(cols(cols'last-cd+1)) then
      if debug
       then put_line(file,"No choice because in rising checker's square.");
      end if;
      return 1;   -- we must swap
    else
      if debug
       then put_line(file,"White checker not in rising checker's square.");
      end if;
    end if;
   -- is there a blocker?
    if Blocker(rows,cols,cr,cd) then
      if debug
       then put_line(file,"No choice because of a blocker.");
      end if;
      return 0;   -- we must stay
    else
      if debug
       then put_line(file,"There is not a blocker, so there is choice.");
      end if;
      return 2;
    end if;
  end Central_Choice;

  function Central_Choice
              ( p : Vector; d,r : integer32;
                rows,cols : Vector; cr,cd : integer32; debug : boolean ) 
              return integer32 is
  begin
   -- is white checker in critical row in descending checker's square?
    if p(d) = rows(cr) and p'last-d+1 = integer32(cols(cols'last-cr+1))
     then return 1;   -- we must swap
    end if;
   -- is top white checker in critical row in rising checker's square?
    if p(r) = rows(cd) and p'last-r+1 = integer32(cols(cols'last-cd+1))
     then return 1;   -- we must swap
    end if;
   -- is there a blocker?
    if Blocker(rows,cols,cr,cd)
     then return 0;   -- we must stay
     else return 2;
    end if;
  end Central_Choice;

  procedure Swap ( rows : in out Vector; i,j : in integer32 ) is

    tmp : natural32;

  begin
    if i /= 0 then
      tmp := rows(j);
      rows(j) := rows(i);
      rows(i) := tmp;
    end if;
  end Swap;

  function Happy_in_Row ( p : Vector; r,c : natural32 ) return boolean is
  begin
    for k in p'range loop 
      if p(k) = r then                   -- black checker in row r ?
        if p'last+1-k <= integer32(c)    -- and column <= c ?
         then return true;
        end if;
      end if;
    end loop;
    return false;
  end Happy_in_Row;

  function Happy_in_Column ( p : Vector; r,c : natural32 ) return boolean is

    k : constant integer32 := p'last+1-integer32(c);

  begin
    return (p(k) <= r);
  end Happy_in_Column;

  procedure Make_Happy ( p : in Vector; rows,cols : in out Vector ) is

    r,c : natural32;

  begin
    for i in rows'range loop
      r := rows(i);
      c := cols(cols'last+1-i);
     -- put("White checker at (");
     -- put(r,1); put(",");
     -- put(c,1); put(")");
      if not Happy_in_Row(p,r,c) then
        -- put_line(" is not happy in its row.");
        while r > 1 loop
          r := r - 1;
          exit when Happy_in_Row(p,r,c);
        end loop;
        rows(i) := r;
       -- else put_line(" is happy in its row.");
      end if;
     -- put("White checker at (");
     -- put(r,1); put(",");
     -- put(c,1); put(")");
      if not Happy_in_Column(p,r,c) then
        -- put_line(" is not happy in its column.");
        while c > 1 loop
          c := c - 1;
          exit when Happy_in_Column(p,r,c);
        end loop;
        cols(cols'last+1-i) := c;
       -- else put_line(" is happy in its column.");
      end if;
    end loop;
  end Make_Happy;

  procedure Write_Position ( r,c : in natural32 ) is

  -- DESCRIPTION :
  --   Auxiliary writing routine when debugging Make_Happy.

  begin
    put("White checker at (");
    put(r,1); put(","); put(c,1); put(")");
  end Write_Position;

  procedure Make_Happy ( p : in Vector; rows,cols : in out Vector;
                         debug : in boolean ) is

    r,c : natural32;

  begin
    if debug
     then put_line("entering Make_Happy ...");
    end if;
    for i in rows'range loop
      r := rows(i);
      c := cols(cols'last+1-i);
      if debug
       then Write_Position(r,c);
      end if;
      if not Happy_in_Row(p,r,c) then
        if debug
         then put_line(" is not happy in its row.");
        end if;
        while (r > 1) loop
          r := r - 1;
          exit when Happy_in_Row(p,r,c);
        end loop;
        rows(i) := r;
      else
        if debug
         then put_line(" is happy in its row.");
        end if;
      end if;
      if debug
       then Write_Position(r,c);
      end if;
      if not Happy_in_Column(p,r,c) then
        if debug
         then put_line(" is not happy in its column.");
        end if;
        while (c > 1) loop
          c := c - 1;
          exit when Happy_in_Column(p,r,c);
        end loop;
        cols(cols'last+1-i) := c;
      else
        if debug
         then put_line(" is happy in its column.");
        end if;
      end if;
    end loop;
    if debug then 
      Write_Position(r,c);
      if Happy_in_Row(p,r,c)
       then put_line(" is happy in its row");
       else put_line(" is NOT happy in its row");
      end if;
      if Happy_in_Column(p,r,c)
       then put_line(" and happy in its column.");
       else put_line(" and NOT happy in its column.");
      end if;
    end if;
    if debug
     then put_line(" ... leaving Make_Happy.");
    end if;
  end Make_Happy;

  procedure Check_Happiness
              ( p,rows,cols : in Vector; happy : out boolean ) is

    r,c : natural32;

  begin
    happy := true;
    for i in rows'range loop
      r := rows(i); c := cols(cols'last+1-i);
      put("White checker at (");
      put(r,1); put(","); put(c,1); put(")");
      if Happy_in_Row(p,r,c)
       then put(" is happy in its row");
       else put(" is NOT happy in its row"); happy := false;
      end if;
      if Happy_in_Column(p,r,c)
       then put_line(" and happy in its column.");
       else put_line(" and NOT happy in its column."); happy := false;
      end if;
    end loop;
  end Check_Happiness;

  function Happy_Checkers ( p,rows,cols : Vector ) return boolean is

    r,c : natural32;

  begin
    for i in rows'range loop
      r := rows(i); 
      c := cols(cols'last+1-i);
      if not Happy_in_Row(p,r,c)
       then return false;
      end if;
      if not Happy_in_Column(p,r,c)
       then return false;
      end if;
    end loop;
    return true;
  end Happy_Checkers;

end Checker_Moves;
