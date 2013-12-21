with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Symbol_Table,Symbol_Table_io;
with Characters_and_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Write_Numbers;
with Standard_Binomial_Varieties;

package body Standard_Binomial_Varieties_io is

-- OUTPUT ROUTINES :

  procedure Fill_Symbol_Table ( n : in natural32 ) is
  begin
    if Symbol_Table.Empty then
      Symbol_Table.Init(n);
      for i in 1..n loop
        declare
          s : Symbol_Table.Symbol;
          c : constant string := Characters_and_Numbers.nConvert(i);
        begin
          s := (s'range => ' ');
          s(1) := 'x';
          for j in c'range loop
            s(j+1) := c(j);
          end loop;
          Symbol_Table.Add(s);
        end;
      end loop;
    end if;
  end Fill_Symbol_Table;

  procedure Write_Header ( file : in file_type; n,d : in natural32 ) is
  begin
    put(file,n,1); put(file," "); put(file,n+d,1); new_line(file);
    for i in 1..(d-1) loop
      put(file,"t"); put(file,i,1); put(file," + ");
    end loop;
    put(file,"t"); put(file,d,1); put(file," - ");
    for i in 1..(d-1) loop
      put(file,"t"); put(file,i,1); put(file," - ");
    end loop;
    put(file,"t"); put(file,d,1); put(file," + "); new_line(file);
  end Write_Header;

  procedure Write_Solution
               ( file : in file_type; d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := M'last(1);
    tc : constant Standard_Complex_Vectors.Vector
       := Standard_Binomial_Varieties.Transform_Coefficients(integer32(d),M,c); 
    one : constant Complex_Number := Create(1.0);
    cnt : natural32 := 0;
    bracket : boolean;
 
  begin
    for j in 1..n loop
      declare
        s : constant Symbol_Table.Symbol := Symbol_Table.get(natural32(j)); 
        first : boolean := true;
      begin
        Symbol_Table_io.put(file,s); put(file," - ");
        if tc(j) /= one then
          bracket := Standard_Write_Numbers.Is_Real(tc(j))
                  or Standard_Write_Numbers.Is_Imag(tc(j));
          if bracket then put(file,"("); end if;
          Standard_Write_Numbers.Write_Number(file,tc(j),cnt);
          if bracket then put(file,")"); end if;
          first := false;
        end if;
        for i in 1..integer32(d) loop
          if M(i,j) /= 0 then
            if first
             then first := false;
             else put(file,"*");
            end if;
            put(file,"t"); put(file,i,1);
            if M(i,j) /= 1 
             then put(file,"^"); put(file,M(i,j),1);
            end if;
          end if;
        end loop;
        if first then -- we have the constant term
          if tc(j) = one
           then put(file,"1");
          end if;
        end if;
        put_line(file,";"); 
      end;
    end loop;
  end Write_Solution;

  procedure Write_System
               ( d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector;
                 p : out Standard_Complex_Laur_Systems.Laur_Sys ) is

    n : constant integer32 := integer32(M'last(1));
    tc : constant Standard_Complex_Vectors.Vector
       := Standard_Binomial_Varieties.Transform_Coefficients(integer32(d),M,c); 
    tm : Term;
 
  begin
    tm.dg := new Standard_Integer_Vectors.Vector(1..n+integer32(d));
    for j in 1..n loop
      for i in 1..n+integer32(d) loop  -- p(j) = x(j) - tc(j)*t^M
        tm.dg(i) := 0;
      end loop;
      tm.dg(integer32(d)+j) := 1;
      tm.cf := Create(1.0);
      p(j) := Create(tm);   -- initialize p(j) with x(j)
      tm.dg(integer32(d)+j) := 0;
      tm.cf := -tc(j);      -- subtract parameter product
      for i in 1..integer32(d) loop
        tm.dg(i) := M(i,j);
      end loop;
      Add(p(j),tm);
    end loop;
    Clear(tm);
  end Write_System;

  procedure Write_Free_Affine_Solution
               ( file : in file_type;
                 s,f : in Standard_Integer_Vectors.Vector ) is

    n : constant integer32 := s'last;
    cnt_free : natural32 := 0;

  begin
    for j in 1..n loop
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.get(natural32(j)); 
      begin
        Symbol_Table_io.put(file,sb);
        if s(j) = 1 then
          put_line(file,";");
        elsif f(j) = 1 then
          cnt_free := cnt_free + 1;
          put(file," - t");
          put(file,cnt_free,1); put_line(file,";");
        end if;
      end;
    end loop;
  end Write_Free_Affine_Solution;

  procedure Write_Free_Affine_System
               ( d : in natural32;
                 f : in Standard_Integer_Vectors.Vector;
                 p : out Standard_Complex_Laur_Systems.Laur_Sys ) is

    n : constant integer32 := f'last;
    cnt_free : natural32 := 0;
    tm : Term;

  begin
    tm.dg := new Standard_Integer_Vectors.Vector'(1..n+integer32(d) => 0);
    tm.cf := Create(1.0);
    for j in 1..n loop
      tm.dg(j+integer32(d)) := 1;
      p(j) := Create(tm);
      tm.dg(j+integer32(d)) := 0;
      if f(j) = 1 then
        cnt_free := cnt_free + 1;
        tm.dg(integer32(cnt_free)) := 1;
        Sub(p(j),tm);
        tm.dg(integer32(cnt_free)) := 0;
      end if;
    end loop;
    Clear(tm);
  end Write_Free_Affine_System;

  procedure Write_Affine_Solution
               ( file : in file_type; d : in natural32;
                 s,f : in Standard_Integer_Vectors.Vector; 
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := s'last;
    tc : constant Standard_Complex_Vectors.Vector
       := Standard_Binomial_Varieties.Transform_Coefficients(integer32(d),M,c); 
    one : constant Complex_Number := Create(1.0);
    ind : integer32 := 0;
    cnt : natural32 := 0;
    cnt_free : natural32 := 0;
    bracket,first : boolean;

  begin
    for j in 1..n loop
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.get(natural32(j)); 
      begin
        Symbol_Table_io.put(file,sb);
        if s(j) = 1 then
          put_line(file,";");
        elsif f(j) = 1 then
          cnt_free := cnt_free + 1;
          put(file," - t");
          put(file,d+cnt_free,1); put_line(file,";");
        else
          ind := ind + 1;
          first := true;
          put(file," - ");
          if tc(ind) /= one then
            bracket := Standard_Write_Numbers.Is_Real(tc(ind))
                    or Standard_Write_Numbers.Is_Imag(tc(ind));
            if bracket then put(file,"("); end if;
            Standard_Write_Numbers.Write_Number(file,tc(ind),cnt);
            if bracket then put(file,")"); end if;
            first := false;
          end if;
          for i in 1..integer32(d) loop
            if M(i,ind) /= 0 then
              if first
               then first := false;
               else put(file,"*");
              end if;
              put(file,"t"); put(file,i,1);
              if M(i,ind) /= 1
               then put(file,"^"); put(file,M(i,ind),1);
              end if;
            end if;
          end loop;
          if first then -- we have the constant term
            if tc(ind) = one
             then put(file,"1");
            end if;
          end if;
          put_line(file,";"); 
        end if;
      end;
    end loop;
  end Write_Affine_Solution;

  procedure Write_Affine_System
               ( d,e : in natural32;
                 s,f : in Standard_Integer_Vectors.Vector; 
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector;
                 p : out Standard_Complex_Laur_Systems.Laur_Sys ) is

    n : constant integer32 := s'last;
    tc : constant Standard_Complex_Vectors.Vector
       := Standard_Binomial_Varieties.Transform_Coefficients(integer32(d),M,c); 
    ind : integer32 := 0;
    cnt_free : natural32 := 0;
    tm : Term;

  begin
    tm.dg := new Standard_Integer_Vectors.Vector(1..n+integer32(d+e));
    for j in 1..n loop
      for i in 1..n+integer32(d+e) loop
        tm.dg(i) := 0;
      end loop;
      tm.dg(integer32(d)+j) := 1;
      tm.cf := Create(1.0);
      p(j) := Create(tm);
      tm.dg(integer32(d)+j) := 0;
      if f(j) = 1 then
        cnt_free := cnt_free + 1;
        tm.dg(integer32(cnt_free)) := 1;
        Sub(p(j),tm);
      elsif s(j) = 0 then
        ind := ind + 1;
        tm.cf := -tc(j);
        for i in 1..integer32(d) loop
          tm.dg(i) := M(i,ind);
        end loop;
      end if;
    end loop;
  end Write_Affine_System;

  procedure Write_Solution
               ( file : in file_type; d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 w : in Standard_Integer_Vectors.Vector;
                 c : in Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := M'last(1);
    tc : constant Standard_Complex_Vectors.Vector
       := Standard_Binomial_Varieties.Transform_Coefficients(integer32(d),M,c); 
    one : constant Complex_Number := Create(1.0);
    cnt : natural32 := 0;
    bracket : boolean;
 
  begin
    for j in 1..n loop
      declare
        s : constant Symbol_Table.Symbol := Symbol_Table.get(natural32(j)); 
        first : boolean := true;
      begin
        Symbol_Table_io.put(file,s); put(file," - ");
        if tc(j) /= one then
          bracket := Standard_Write_Numbers.Is_Real(tc(j))
                  or Standard_Write_Numbers.Is_Imag(tc(j));
          if bracket then put(file,"("); end if;
          Standard_Write_Numbers.Write_Number(file,tc(j),cnt);
          if bracket then put(file,")"); end if;
          first := false;
        end if;
        for i in 1..integer32(d) loop
          if M(i,j) /= 0 then
            if first
             then first := false;
             else put(file,"*");
            end if;
            put(file,"t"); put(file,i,1); put(file,"^");
            put(file,M(i,j),1);
            if w(i) /= 1 
             then put(file,"/"); put(file,w(i),1);
            end if;
          end if;
        end loop;
        if first then -- we have the constant term
          if tc(j) = one
           then put(file,"1");
          end if;
        end if;
        put_line(file,";"); 
      end;
    end loop;
  end Write_Solution;

-- INPUT ROUTINES :

  function Variable_Term ( t : Term; n,d : natural32 ) return boolean is
  begin
    for i in integer32(d+1)..integer32(d+n) loop
      if t.dg(i) /= 0
       then return true;
      end if;
    end loop;
    return false;
  end Variable_Term;

  procedure Extract_Binomial_Variety
               ( p : in Poly; n,d,i : in natural32;
                 T : out Standard_Integer_Matrices.Matrix;
                 c : out Standard_Complex_Vectors.Vector ) is

    procedure Visit ( m : in Term; continue : out boolean ) is
    begin
      if not Variable_Term(m,n,d) then
        c(integer32(i)) := -m.cf;
        for j in 1..integer32(d) loop
          T(integer32(i),j) := m.dg(j);
        end loop;
      end if;
      continue := true;
    end Visit;
    procedure Visit_Terms is new Visiting_Iterator(Visit);
 
  begin
    Visit_Terms(p);
  end Extract_Binomial_Variety;

  procedure Extract_Binomial_Variety
               ( s : in Laur_Sys; n,d : in natural32;
                 T : out Standard_Integer_Matrices.Matrix;
                 c : out Standard_Complex_Vectors.Vector ) is
  begin
    for i in 1..integer32(n) loop
      Extract_Binomial_Variety(s(i),n,d,natural32(i),T,c);
      if Number_of_Terms(s(i)) = 1 then  -- zero component
        c(i) := Create(0.0);
        for j in T'range(2) loop
          T(i,j) := 0;
        end loop;
      end if;
    end loop;
  end Extract_Binomial_Variety;

  procedure Parse_Binomial_Variety
              ( p : in Laur_Sys;
                T : out Standard_Integer_Matrices.Link_to_Matrix;
                c : out Standard_Complex_Vectors.Link_to_Vector ) is

    n : constant integer32 := p'last;
    d : constant integer32 := integer32(Number_of_Unknowns(p(p'first))) - n;
    Tm : Standard_Integer_Matrices.Matrix(1..n,1..d);
    cv : Standard_Complex_Vectors.Vector(1..n);

  begin
    Extract_Binomial_Variety(p,natural32(n),natural32(d),Tm,cv);
    T := new Standard_Integer_Matrices.Matrix'(Tm);
    c := new Standard_Complex_Vectors.Vector'(cv);
  end Parse_Binomial_Variety;

end Standard_Binomial_Varieties_io;
