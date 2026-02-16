with String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Complex_Numbers;
with Characters_and_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;
with Standard_Complex_Poly_Strings;
with Standard_Complex_Laur_Readers;
with Standard_Complex_Laur_Strings;
with Standard_Complex_Laurentials_IO;
with Real_Powered_Series_IO;

package body Real_Powered_Homotopy_IO is

-- WRITE OUTPUT :

  function to_string
             ( d : Standard_Complex_Laurentials.Degrees ) return string is

  -- DESCRIPTION :
  --   Returns the string representation of the terms defined
  --   by the degrees d.

    function write ( k : integer32; accu : string ) return string is

    -- DESCRIPTION :
    --   Recursive write where k is the position in the degrees d.

    begin
      if k > d'last then
        return accu;
      elsif d(k) = 0 then
        return write(k+1,accu);
      else
        declare
          sb : constant Symbol_Table.Symbol := Symbol_Table.Get(natural32(k));
        begin
          if d(k) = 1 then
            declare
              new_accu : constant string := accu & '*'
                & Standard_Complex_Poly_Strings.Concat_Symbol0("",sb);
            begin
              return write(k+1,new_accu);
            end;
          else
            declare
              new_accu : constant string := accu & '*'
                & Standard_Complex_Poly_Strings.Concat_Symbol0("",sb)
                & '^' & Characters_and_Numbers.Convert(d(k));
            begin
              return write(k+1,new_accu);
            end;
          end if;
        end;
      end if;
    end write;

  begin
    return write(d'first,"");
  end to_string;

  function to_string ( q : Standard_Complex_Laurentials.Poly;
                       c : Standard_Complex_VecVecs.VecVec;
                       p : Standard_Floating_VecVecs.VecVec;
                       t : character := 't'; vrblvl : integer32 := 0 )
                     return string is

    terms : array(c'range) of Standard_Complex_Laurentials.Term;
    idx : integer32 := 0;

    procedure Get_Term ( trm : in Standard_Complex_Laurentials.Term;
                         continue : out boolean ) is

    -- DESCRIPTION :
    --   Stores the current term trm in the array terms at position idx
    --   and continues if continue is set to true.

    begin
      idx := idx + 1;
      terms(idx) := trm;
      continue := true;
    end Get_Term;

    procedure Get_Terms is new 
      Standard_Complex_Laurentials.Visiting_Iterator(Get_Term);

    function Write ( k : integer32; accu : string ) return string is

    -- DESCRIPTION :
    --   Recursive function to write the terms, stops when k > c'last,
    --   returning the accumulator accu.

    begin
      if k > c'last then
        return accu;
      else
        declare
          sk : constant string
             := Real_Powered_Series_IO.to_string(c(k).all,p(k).all,t,vrblvl-1);
          str_trm : constant string := to_string(terms(k).dg);
        begin
          if k = 1 then -- do not write the first +
            declare
              new_accu : constant string := accu & "(" & sk & ")";
            begin
              return Write(k+1,new_accu & str_trm);
            end;
          else
            declare
              new_accu : constant string := accu & " + (" & sk & ")";
            begin
              return Write(k+1,new_accu & str_trm);
            end;
          end if;
        end;
      end if;
    end Write;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.to_string ...");
    end if;
    Get_Terms(q);
    return Write(1,"");
  end to_string;

  procedure put ( q : in Standard_Complex_Laurentials.Poly;
                  c : in Standard_Complex_VecVecs.VecVec;
                  p : in Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put 0 ...");
    end if;
    put(standard_output,q,c,p,t,vrblvl);
  end put;

  procedure put ( file : in file_type;
                  q : in Standard_Complex_Laurentials.Poly;
                  c : in Standard_Complex_VecVecs.VecVec;
                  p : in Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;

    procedure Write_Term ( trm : in Standard_Complex_Laurentials.Term;
                           continue : out boolean ) is

    -- DESCRIPTION :
    --   Writes the term trm, and continues if continue is set to true.

      zerodeg : boolean := true;

    begin
      if idx = 0
       then put(file,"(");
       else put(file," + (");
      end if;
      idx := idx + 1;
      Real_Powered_Series_io.put(file,c(idx).all,p(idx).all,t);
      for i in trm.dg'range loop
        if trm.dg(i) /= 0
         then zerodeg := false; exit;
        end if;
      end loop;
      if zerodeg then
        put(file,")");
      else
        put(file,")*");
        Standard_Complex_Laurentials_IO.put
          (file,trm.dg,std => true, pow => '^');
      end if;
      continue := true;
    end Write_Term;

    procedure Write_Terms is
      new Standard_Complex_Laurentials.Visiting_Iterator(Write_Term);

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put 1 ...");
    end if;
    Write_Terms(q); put(file,";");
  end put;

  procedure put_line ( q : in Standard_Complex_Laurentials.Poly;
                       c : in Standard_Complex_VecVecs.VecVec;
                       p : in Standard_Floating_VecVecs.VecVec;
                       t : in character := 't';
                       vrblvl : in integer32 := 0 ) is

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put_line 0 ...");
    end if;
    put_line(standard_output,q,c,p,t,vrblvl);
  end put_line;

  procedure put_line ( file : in file_type;
                       q : in Standard_Complex_Laurentials.Poly;
                       c : in Standard_Complex_VecVecs.VecVec;
                       p : in Standard_Floating_VecVecs.VecVec;
                       t : in character := 't';
                       vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;

    procedure Write_Term ( trm : in Standard_Complex_Laurentials.Term;
                           continue : out boolean ) is

    -- DESCRIPTION :
    --   Writes the term trm, and continues if continue is set to true.

      zerodeg : boolean := true;

    begin
      if idx = 0
       then put(file,"(");
       else put(file," + (");
      end if;
      idx := idx + 1;
      Real_Powered_Series_io.put_line(file,c(idx).all,p(idx).all,t);
      for i in trm.dg'range loop
        if trm.dg(i) /= 0
         then zerodeg := false; exit;
        end if;
      end loop;
      if zerodeg then
        put(file,")");
      else
        put(file,")*");
        Standard_Complex_Laurentials_IO.put
          (file,trm.dg,std => true, pow => '^');
      end if;
      continue := true;
    end Write_Term;

    procedure Write_Terms is
      new Standard_Complex_Laurentials.Visiting_Iterator(Write_Term);

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put_line 1 ...");
    end if;
    Write_Terms(q);
    put_line(file,";");
  end put_line;

-- WRITE SYSTEMS :

  procedure put ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                  c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                  p : in Standard_Floating_VecVecs.Array_of_VecVecs;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put 2 ...");
    end if;
    put(standard_output,q,c,p,t,vrblvl);
  end put;

  procedure put ( file : in file_type;
                  q : in Standard_Complex_Laur_Systems.Laur_Sys;
                  c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                  p : in Standard_Floating_VecVecs.Array_of_VecVecs;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put 3 ...");
    end if;
    for i in q'range loop
      put(file,q(i),c(i).all,p(i).all,t,vrblvl-1);
    end loop;
  end put;

  procedure put_line ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                       c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                       p : in Standard_Floating_VecVecs.Array_of_VecVecs;
                       t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put_line 2 ...");
    end if;
    put_line(standard_output,q,c,p,t,vrblvl);
  end put_line;

  procedure put_line ( file : in file_type;
                       q : in Standard_Complex_Laur_Systems.Laur_Sys;
                       c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                       p : in Standard_Floating_VecVecs.Array_of_VecVecs;
                       t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put_line 3 ...");
    end if;
    for i in q'range loop
      put_line(file,q(i),c(i).all,p(i).all,t,vrblvl-1);
    end loop;
  end put_line;

  procedure put ( npol,nvar,size : in integer32;
                  q : in Standard_Complex_Laur_Systems.Laur_Sys;
                  c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                  p : in Standard_Floating_VecVecs.Array_of_VecVecs;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put 4 ...");
    end if;
    put(standard_output,npol,nvar,size,q,c,p,t,vrblvl);
  end put;

  procedure put ( file : in file_type; npol,nvar,size : in integer32;
                  q : in Standard_Complex_Laur_Systems.Laur_Sys;
                  c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                  p : in Standard_Floating_VecVecs.Array_of_VecVecs;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put 5 ...");
    end if;
    put(npol,1); put(" ");
    put(nvar,1); put(" ");
    put(size,1); new_line;
    put_line(file,q,c,p,t,vrblvl);
  end put;

  procedure put_line ( npol,nvar,size : in integer32;
                       q : in Standard_Complex_Laur_Systems.Laur_Sys;
                       c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                       p : in Standard_Floating_VecVecs.Array_of_VecVecs;
                       t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put_line 4 ...");
    end if;
    put_line(standard_output,npol,nvar,size,q,c,p,t,vrblvl);
  end put_line;

  procedure put_line ( file : in file_type; npol,nvar,size : in integer32;
                       q : in Standard_Complex_Laur_Systems.Laur_Sys;
                       c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                       p : in Standard_Floating_VecVecs.Array_of_VecVecs;
                       t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put_line 5 ...");
    end if;
    put(file,npol,1); put(file," ");
    put(file,nvar,1); put(file," ");
    put(file,size,1); new_line(file);
    put_line(file,q,c,p,t,vrblvl);
  end put_line;

-- PARSE INPUT :

  function number_of_terms
             ( s : string; vrblvl : integer32 := 0 ) return integer32 is

    cnt : integer32 := 0;
    res : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.parse_string ...");
    end if;
    for i in s'range loop
      if s(i) = '(' then
        cnt := cnt + 1;
      elsif s(i) = ')' then
        cnt := cnt - 1;
        if cnt = 0
         then res := res + 1;
        end if;
        if cnt < 0 then
          if vrblvl > 0
           then put_line("Bracket count mismatch.  Error, returning -1 ...");
          end if;
          return -1;
        end if;
      end if;
    end loop;
    return res;
  end number_of_terms;

  procedure parse_series
              ( s : in string;
                c : out Standard_Complex_VecVecs.VecVec;
                p : out Standard_Floating_VecVecs.VecVec;
                t : in character := 't'; vrblvl : in integer32 := 0 ) is

    cnt : integer32 := 0;
    i1,i2 : integer := 0;
    idx : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.parse_series ...");
    end if;
    for i in s'range loop
      if s(i) = '(' then
        cnt := cnt + 1;
        if cnt = 1
         then i1 := integer(i) +1;
        end if;
      elsif s(i) = ')' then
        cnt := cnt - 1;
        if cnt = 0 then
          i2 := integer(i) - 1;
          idx := idx + 1;
          if vrblvl > 0 then
            put("parsing " & s(i1..i2));
            put(", assigning as series "); put(idx,1); new_line;
          end if;
          Real_Powered_Series_IO.parse_string
            (s(i1..i2),c(idx),p(idx),t,vrblvl-1);
        end if;
      end if;
    end loop;
  end parse_series;

  function extract_polynomial_string
             ( s : in string; m : in integer32;
               vrblvl : in integer32 := 0 ) return string is

    cnt : integer32 := 0;
    idx : integer := s'first;

    function extract ( k : integer32; accu : string ) return string is

    -- DESCRIPTION :
    --   Recursive removal of the series coefficients,
    --   where k runs over all terms.
    --   The accumulator accu is returned when k > m.

      newterm : boolean;
      startidx,endidx : integer;

    begin
      if k > m then
        return accu & '1' & s(idx..s'last);
      else
        newterm := false;
        startidx := idx;
        loop
          if s(idx) = '(' then
            cnt := cnt + 1;
            if cnt = 1
             then endidx := idx - 1; -- do not include '('
            end if;
          elsif s(idx) = ')' then
            cnt := cnt - 1;
            if cnt = 0
             then newterm := true;
            end if;
          end if;
          idx := idx + 1; -- skip the last ')'
          exit when newterm;
        end loop;
        if k = 1 then
          declare
            new_accu : constant string := accu & s(startidx..endidx);
          begin
            return extract(k+1,new_accu);
          end;
        else
          declare
            new_accu : constant string := accu & '1' & s(startidx..endidx);
          begin
            return extract(k+1,new_accu);
          end;
        end if;
      end if;
    end extract;

  begin
    if vrblvl > 0 then
      put_line("-> in Real_Powered_Homotopy_IO.extract_polynomial_string ...");
    end if;
    return extract(1,"");
  end extract_polynomial_string;

  procedure parse_string
              ( s : in string; n : in integer32;
                q : out Standard_Complex_Laurentials.Poly;
                c : out Standard_Complex_VecVecs.VecVec;
                p : out Standard_Floating_VecVecs.VecVec;
                t : in character := 't'; vrblvl : in integer32 := 0 ) is

    m : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.parse_string ...");
    end if;
    m := number_of_terms(s,vrblvl-1);
    if vrblvl > 0
     then put("number of terms : "); put(m,1); new_line;
    end if;
    parse_series(s,c,p,t,vrblvl);
    declare
      sq : constant string
         := extract_polynomial_string(s,m,vrblvl) & ';';
    begin
      if vrblvl > 0 then
        put_line("sq : " & sq);
      end if;
      q := Standard_Complex_Laur_Strings.parse(natural32(n),sq);
    end;
  end parse_string;

  procedure get ( n,size : in integer32;
                  q : out Standard_Complex_Laurentials.Poly;
                  c : out Standard_Complex_VecVecs.VecVec;
                  p : out Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.get 0 ...");
    end if;
    get(standard_input,n,size,q,c,p,t,vrblvl);
  end get;

  procedure get ( file : in file_type; n,size : in integer32;
                  q : out Standard_Complex_Laurentials.Poly;
                  c : out Standard_Complex_VecVecs.VecVec;
                  p : out Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is

    ch : character;
    idx : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.get 1 ...");
    end if;
    while not end_of_file(file) loop
      loop
        get(file,ch);
        exit when end_of_file(file) or ch = '(';
      end loop;
      exit when end_of_file(file);
      idx := idx + 1;
      if vrblvl > 0 then
        put("reading series "); put(idx,1); put_line(" ...");
      end if;
      Real_Powered_Series_IO.get(file,size,c(idx),p(idx),t,vrblvl-1);
     -- the current character read is ')' the closing of the series
      ch := ' ';
      declare
        cnt : integer32 := 0;
        trm : Standard_Complex_Laurentials.Term;
        pol : Standard_Complex_Laurentials.Poly;
      begin
        trm.cf := Standard_Complex_Numbers.create(1.0);
        trm.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
        while not end_of_file(file) loop
          get(file,ch);
          exit when end_of_file(file) or ch = '*';
        end loop;
        exit when end_of_file(file);
        get(file,ch); -- skip '*'
        Standard_Complex_Laur_Readers.Read_Factor
          (file,cnt,ch,natural32(n),trm.dg,pol);
        if vrblvl > 0 then
          put("trm : "); Standard_Complex_Laurentials_IO.put(trm); new_line;
         -- put("pol : "); Standard_Complex_Laurentials_IO.put(pol); new_line;
        end if;
        Standard_Complex_Laurentials.Add(q,trm);
        Standard_Complex_Laurentials.Clear(trm);
      end;
    end loop;
  end get;

-- READING SYSTEMS :

  procedure get ( npol,nvar : out integer32;
                  q : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                  c : out Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                  p : out Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.get 2 ...");
    end if;
    get(standard_input,npol,nvar,q,c,p,t,vrblvl);
  end get;

  procedure get ( file : in file_type; npol,nvar : out integer32;
                  q : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                  c : out Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                  p : out Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is

    n,m : natural;
    sys : String_Splitters.Link_to_Array_of_Strings;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.get 3 ...");
    end if;
    String_Splitters.get(file,n,m,sys);
    if vrblvl > 0 then
      put("Read "); put(integer32(n),1); put(" strings, ");
      put(integer32(m),1); put_line(" .. as number of variables ...");
    end if;
    npol := integer32(n);
    nvar := integer32(m);
    if Symbol_Table.Empty
     then Symbol_Table.Init(natural32(nvar));
    end if;
    declare
      sq : Standard_Complex_Laur_Systems.Laur_Sys(1..npol);
      cf : Standard_Complex_VecVecs.Array_of_VecVecs(1..npol);
      pw : Standard_Floating_VecVecs.Array_of_VecVecs(1..npol);
    begin
      for i in 1..npol loop
        declare
          nbt : constant integer32
              := number_of_terms(sys(integer(i)).all,vrblvl-1);
          cfi : Standard_Complex_VecVecs.VecVec(1..nbt);
          pwi : Standard_Floating_VecVecs.VecVec(1..nbt);
        begin
          parse_string(sys(integer(i)).all,nvar,sq(i),cfi,pwi,t,vrblvl-1);
          cf(i) := new Standard_Complex_VecVecs.VecVec'(cfi);
          pw(i) := new Standard_Floating_VecVecs.VecVec'(pwi);
        end;
      end loop;
      q := new Standard_Complex_Laur_Systems.Laur_Sys'(sq);
      c := new Standard_Complex_VecVecs.Array_of_VecVecs'(cf);
      p := new Standard_Floating_VecVecs.Array_of_VecVecs'(pw);
    end;
  end get;

  procedure get ( q : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                  c : out Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                  p : out Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.get 4 ...");
    end if;
    get(standard_input,q,c,p,t,vrblvl); -- not so useful ...
  end get;

  procedure get ( file : in file_type;
                  q : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                  c : out Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                  p : out Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is

    npol,nvar : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.get 5 ...");
    end if;
    get(file,npol,nvar,q,c,p,t,vrblvl);
  end get;

end Real_Powered_Homotopy_IO;
