with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;

package body Checker_Boards_io is

-- AUXILIARY ROUTINES :

  procedure Write_Line ( file : in file_type;
                         b : in Board; i : in integer32 ) is

  -- DESCRIPTION :
  --   Auxiliary routine to write the i-th line of a board b to file.

  begin
    put(file,i,2);
    for j in b'range(2) loop
      if b'last(1) < 10 
       then put(file," ");
       else put(file,"  ");
      end if;
      put(file,b(i,j));
    end loop;
  end Write_Line;

  procedure Write_Line ( file : in file_type;
                         m : in Matrix; i : in integer32 ) is

  -- DESCRIPTION :
  --   Writes the i-th line of a pattern in m to file.

  begin
    for j in m'range(2) loop
      if m(i,j) = 2
       then put(file," *");
       else put(file,m(i,j),2);
      end if;
    end loop;
  end Write_Line;

-- TARGET ROUTINES :
 
  procedure Write ( p : in Vector ) is
  begin
    Write(standard_output,p);
  end Write;

  procedure Write ( file : in file_type; p : in Vector ) is
  begin
    if p'last < 10 then
      put(file,p);
    else
      for i in p'range loop
        if p(i) < 10
         then put(file,"  ");
         else put(file," ");
        end if;
        put(file,p(i),1);
      end loop;
    end if;
  end Write;

  procedure Write ( b : in Board ) is
  begin
    Write(standard_output,b);
  end Write;

  procedure Write ( file : in file_type; b : in Board ) is
  begin
    for i in b'range(1) loop
      Write_Line(file,b,i);
      new_line(file);
    end loop;
  end Write;

  procedure Write ( a,b : in Board ) is
  begin
    Write(standard_output,a,b);
  end Write;

  procedure Write ( file : in file_type; a,b : in Board ) is
  begin
    for i in b'range(1) loop
      Write_Line(file,a,i);
      Write_Line(file,b,i);
      new_line(file);
    end loop;
  end Write;

  procedure Write ( b : in Board; f : in Matrix ) is
  begin
    Write(standard_output,b,f);
  end Write;

  procedure Write ( file : in file_type; b : in Board; f : in Matrix ) is
  begin
    for i in b'range(1) loop
      Write_Line(file,b,i);
      put(file," |");
      for j in f'range(2) loop
        put(file," "); put(file,f(i,j),1);
      end loop;
      new_line(file);
    end loop;
  end Write;

  procedure Write ( b : in Board; f,t : in Matrix ) is
  begin
    Write(standard_output,b,f,t);
  end Write;

  procedure Write ( file : in file_type; b : in Board; f,t : in Matrix ) is

    ft : constant Matrix(f'range(1),t'range(2)) := f*t;

  begin
    for i in b'range(1) loop
      Write_Line(file,b,i);
      put(file," |");
      for j in f'range(2) loop
        put(file," "); put(file,f(i,j),1);
      end loop;
      put(file," T");
      for j in t'range(2) loop
        put(file," "); put(file,t(i,j),1);
      end loop;
      put(file," !");
      for j in ft'range(2) loop
        put(file," "); put(file,ft(i,j),1);
      end loop;
      new_line(file);
    end loop;
  end Write;

  procedure Write ( b : in Board; f,t,x : in Matrix ) is
  begin
    Write(standard_output,b,f,t,x);
  end Write;

  procedure Write ( file : in file_type; b : in Board; f,t,x : in Matrix ) is

    ft : constant Matrix(f'range(1),t'range(2)) := f*t;

  begin
    for i in b'range(1) loop
      Write_Line(file,b,i);
      put(file," |");
      for j in x'range(2) loop
        put(file," "); put(file,x(i,j),1);
      end loop;
      put(file," |");
      for j in f'range(2) loop
        put(file," "); put(file,f(i,j),1);
      end loop;
      put(file," T");
      for j in t'range(2) loop
        put(file," "); put(file,t(i,j),1);
      end loop;
      put(file," !");
      for j in ft'range(2) loop
        put(file," "); put(file,ft(i,j),1);
      end loop;
      new_line(file);
    end loop;
  end Write;
  
  procedure Read_Permutation ( p : out Vector ) is
  begin
    put("Give "); put(p'last,1);
    put(" natural numbers : ");
    -- get(p);  -- does not work with all validity checks !!!
    -- ts_flagcond option #6, crashes 2nd time when reading cols
    for i in p'range loop
      p(i) := 0; get(p(i));
    end loop;
  end Read_Permutation;

  procedure Write_Permutation ( p : in Vector ) is
  begin
    Write_Permutation(standard_output,p);
  end Write_Permutation;

  procedure Write_Permutation ( file : in file_type; p : in Vector ) is

    b : constant Board(p'range,p'range) := Configuration(p);

  begin
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b);
    put(file,"  "); Write(file,Permutation(b)); new_line(file);
  end Write_Permutation;

  procedure Write_Permutation ( p : in Vector; f : in Matrix ) is
  begin
    Write_Permutation(standard_output,p,f);
  end Write_Permutation;

  procedure Write_Permutation ( file : in file_type;
                                p : in Vector; f : in Matrix ) is

    b : constant Board(p'range,p'range) := Configuration(p);

  begin
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b,f);
    put(file,"  "); Write(file,Permutation(b)); new_line(file);
  end Write_Permutation;

  procedure Write_Permutation ( p : in Vector; f,t : in Matrix ) is
  begin
    Write_Permutation(standard_output,p,f,t);
  end Write_Permutation;

  procedure Write_Permutation ( file : in file_type;
                                p : in Vector; f,t : in Matrix ) is

    b : constant Board(p'range,p'range) := Configuration(p);

  begin
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b,f,t);
    put(file,"  "); Write(file,Permutation(b)); new_line(file);
  end Write_Permutation;

  procedure Write_Permutation ( p : in Vector; f,t,x : in Matrix ) is
  begin
    Write_Permutation(standard_output,p,f,t,x);
  end Write_Permutation;

  procedure Write_Permutation ( file : in file_type;
                                p : in Vector; f,t,x : in Matrix ) is

    b : constant Board(p'range,p'range) := Configuration(p);

  begin
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b,f,t,x);
    put(file,"  "); Write(file,Permutation(b)); new_line(file);
  end Write_Permutation;

  procedure Write_Permutation ( p,r,c : in Vector; f : in Matrix ) is
  begin
    Write_Permutation(standard_output,p,r,c,f);
  end Write_Permutation;

  procedure Write_Permutation ( file : in file_type;
                                p,r,c : in Vector; f : in Matrix ) is

    b : Board(p'range,p'range) := Configuration(p);

  begin
    Place_White(b,r,c);
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b,f);
    put(file,"  "); Write_Coordinates(file,r,c);
  end Write_Permutation;

  procedure Write_Permutation ( p,r,c : in Vector; f,t : in Matrix ) is
  begin
    Write_Permutation(standard_output,p,r,c,f,t);
  end Write_Permutation;

  procedure Write_Permutation ( file : in file_type;
                                p,r,c : in Vector; f,t : in Matrix ) is

    b : Board(p'range,p'range) := Configuration(p);

  begin
    Place_White(b,r,c);
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b,f,t);
    put(file,"  "); Write_Coordinates(file,r,c);
  end Write_Permutation;

  procedure Write_Permutation ( p,r,c : in Vector; f,t,x : in Matrix ) is
  begin
    Write_Permutation(standard_output,p,r,c,f,t,x);
  end Write_Permutation;

  procedure Write_Permutation ( file : in file_type;
                                p,r,c : in Vector; f,t,x : in Matrix ) is

    b : Board(p'range,p'range) := Configuration(p);

  begin
    Place_White(b,r,c);
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b,f,t,x);
    put(file,"  "); Write_Coordinates(file,r,c);
  end Write_Permutation;

  procedure Write_Bracket ( v : in Vector ) is
  begin
    Write_Bracket(standard_output,v);
  end Write_Bracket;

  procedure Write_Bracket ( file : in file_type; v : in Vector ) is
  begin
    put(file,"[");
    put(file,v(v'first),1);
    for i in v'first+1..v'last loop
      put(file," "); put(file,v(i),1);
    end loop;
    put(file,"]");
  end Write_Bracket;

  procedure Write_Coordinates ( rows,cols : in Vector ) is
  begin
    Write_Coordinates(standard_output,rows,cols);
  end Write_Coordinates;

  procedure Write_Coordinates ( file : in file_type; rows,cols : in Vector ) is
  begin
    Write_Bracket(file,rows);
    Write_Bracket(file,cols);
    new_line(file);
  end Write_Coordinates;

  procedure Write ( b : in Board; p,rows,cols : in Vector ) is
  begin
    Write(standard_output,b,p,rows,cols);
  end Write;

  procedure Write ( file : in file_type;
                    b : in Board; p,rows,cols : in Vector ) is
  begin
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b);
    put(file,"  "); Write_Coordinates(file,rows,cols);
  end Write;

  procedure Write ( b : in Board; f : in Matrix; p,rows,cols : in Vector ) is
  begin
    Write(standard_output,b,f,p,rows,cols);
  end Write;

  procedure Write ( file : in file_type;
                    b : in Board; f : in Matrix; p,rows,cols : in Vector ) is
  begin
    put(file,"  "); Write(file,p); new_line(file);
    Write(file,b,f);
    put(file,"  "); Write_Coordinates(file,rows,cols);
  end Write;

  procedure Write_Permutations ( p,q,pr,pc,qr,qc : in Vector ) is
  begin
    Write_Permutations(standard_output,p,q,pr,pc,qr,qc);
  end Write_Permutations;

  procedure Write_Permutations
              ( file : in file_type; p,q,pr,pc,qr,qc : in Vector ) is

    a : Board(p'range,p'range) := Configuration(p);
    b : Board(q'range,q'range) := Configuration(q);

  begin
    Place_White(a,pr,pc);
    Place_White(b,qr,qc);
    put(file,"  "); Write(file,p);
    put(file,"  "); Write(file,q); new_line(file);
    Write(file,a,b);
    Write_Bracket(file,pr); Write_Bracket(file,pc);
    put(file,"  ");
    Write_Bracket(file,qr); Write_Bracket(file,qc);
    new_line(file);
  end Write_Permutations;

  procedure Write_Patterns ( fp,fq,xp,xq : in Matrix ) is
  begin
    Write_Patterns(standard_output,fp,fq,xp,xq);
  end Write_Patterns;

  procedure Write_Patterns ( file : in file_type; fp,fq,xp,xq : in Matrix ) is
  begin
    for i in fp'range loop
      Write_Line(file,fp,i); put(file," !");
      Write_Line(file,xp,i); put(file," |");
      Write_Line(file,fq,i); put(file," !");
      Write_Line(file,xq,i); new_line(file);
    end loop;
  end Write_Patterns;

  procedure Write_Permutations_and_Patterns
              ( p,q,pr,pc,qr,qc : in Vector; fp,fq,xp,xq : in Matrix ) is
  begin
    Write_Permutations_and_Patterns
      (standard_output,p,q,pr,pc,qr,qc,fp,fq,xp,xq);
  end Write_Permutations_and_Patterns;

  procedure Write_Permutations_and_Patterns
              ( file : in file_type;
                p,q,pr,pc,qr,qc : in Vector; fp,fq,xp,xq : in Matrix ) is

    a : Board(p'range,p'range) := Configuration(p);
    b : Board(q'range,q'range) := Configuration(q);

  begin
    Place_White(a,pr,pc);
    Place_White(b,qr,qc);
    put(file,"  "); Write(file,p);
    put(file,"  "); Write(file,q); new_line(file);
    for i in p'range loop
      Write_Line(file,a,i); Write_Line(file,b,i); put(file," |");
      Write_Line(file,fp,i); put(file," !");
      Write_Line(file,xp,i); put(file," |");
      Write_Line(file,fq,i); put(file," !");
      Write_Line(file,xq,i); new_line(file);
    end loop;
    Write_Bracket(file,pr); Write_Bracket(file,pc);
    put(file,"  ");
    Write_Bracket(file,qr); Write_Bracket(file,qc);
    new_line(file);
  end Write_Permutations_and_Patterns;

end Checker_Boards_io;
