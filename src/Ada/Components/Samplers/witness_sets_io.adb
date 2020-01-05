with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Symbol_Table_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Permutations;                       use Permutations;
with Permute_Operations;                 use Permute_Operations;               
with Planes_and_Polynomials;             use Planes_and_Polynomials;
with Witness_Sets;                       use Witness_Sets;

package body Witness_Sets_io is

  function Maximum ( a,b : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximum of a and b.

  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Maximum;

-- AUXILIARIES TO DETERMINE THE ELIMINATION ORDER :

  procedure Write_Symbol_Table ( file : in file_type ) is
 
  -- DESCRPITION :
  --   Writes the content of the symbol table in the current order.

  begin
    for i in 1..Symbol_Table.Number loop
      Symbol_Table_io.put(file,Symbol_Table.Get(i)); put(file," ");
    end loop;
    new_line;
  end Write_Symbol_Table;
 
  procedure Write_Symbol_Order ( perm : in Permutation ) is
 
  -- DESCRIPTION :
  --   Writes the symbols according to the given permutation.
 
  begin
    for i in 1..integer32(Symbol_Table.Number) loop
      Symbol_Table_io.put(Symbol_Table.Get(natural32(perm(i)))); put(" ");
    end loop;
  end Write_Symbol_Order;
 
  procedure Permute_Symbol_Table ( perm : in Permutation ) is
 
  -- DESCRIPTION :
  --   Permutes the symbols in the symbol table.
 
    nb : constant natural32 := Symbol_Table.Number;
    syms : array(1..integer32(nb)) of Symbol;
 
  begin
    for i in syms'range loop
      syms(i) := Symbol_Table.Get(natural32(i));
    end loop;
    for i in syms'range loop
      Symbol_Table.Replace(natural32(i),syms(perm(i)));
    end loop;
  end Permute_Symbol_Table;

  procedure Read_Permutation ( perm : out Permutation ) is

    nb : constant integer32 := integer32(Symbol_Table.Number);
    sb : Symbol;
    error : boolean;

  begin
    for i in 1..nb loop
      loop
        put("Give symbol at place "); put(i,1); put(" : ");
        Symbol_Table_io.Get(sb);
        put("  Your symbol : "); Symbol_Table_io.put(sb);
        perm(i) := integer32(Symbol_Table.Get(sb));
        error := false;
        if perm(i) = 0 then
          put_line(" did not occur before!  Please try again...");
          error := true;
        else
          for j in 1..i-1 loop
            if perm(j) = perm(i) then
              error := true;
              put(" occurs already at new place ");
              put(j,1); put_line("!  Please try again...");
            end if;
            exit when error;
          end loop;
          if not error
           then put(" was at position "); put(perm(i),1); new_line;
          end if;
        end if;
        exit when not error;
      end loop;
    end loop;
  end Read_Permutation;

  function Interactive_Read_Permutation
             ( nbunk : in natural32 ) return Standard_Integer_Vectors.Vector is
 
  -- DESCRIPTION :
  --   Returns a permutation, interactively determined by the user.
  --   The loop allows to adjust and to restart.
 
    perm : Standard_Integer_Vectors.Vector(1..integer32(nbunk));
    ans : character;
 
  begin
    loop
      new_line;
      put_line("Reading a permutation of the unknowns.");
      Read_Permutation(Permutation(perm));
      put("The permutation : "); put(perm); new_line;
      put("The new order of symbols : ");
      Write_Symbol_Order(Permutation(perm));
      new_line;
      put("Happy with this new order ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans = 'y');
    end loop;
    return perm;
  end Interactive_Read_Permutation;          

  procedure Permute_Solutions
               ( perm : in Permutation;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is
 
  -- DESCRIPTION :
  --   Permutes the solutions in the list.
 
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ls.v := perm*ls.v;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Permute_Solutions;

  procedure Permute_Solutions
               ( perm : in Permutation;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List ) is
 
  -- DESCRIPTION :
  --   Permutes the solutions in the list.
 
    use DoblDobl_Complex_Solutions;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ls.v := perm*ls.v;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Permute_Solutions;

  procedure Permute_Solutions
               ( perm : in Permutation;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List ) is
 
  -- DESCRIPTION :
  --   Permutes the solutions in the list.
 
    use QuadDobl_Complex_Solutions;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      ls.v := perm*ls.v;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Permute_Solutions;

  procedure Permute_Polynomial
              ( perm : in Permutation;
                f : in out Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      pt : Term := perm*t;

    begin
      Add(res,pt);
      Clear(pt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(f);
    Clear(f);
    f := res;
  end Permute_Polynomial;

  procedure Permute_Polynomial
              ( perm : in Permutation;
                f : in out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      pt : Term := perm*t;

    begin
      Add(res,pt);
      Clear(pt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(f);
    Clear(f);
    f := res;
  end Permute_Polynomial;

  procedure Permute_Polynomial
              ( perm : in Permutation;
                f : in out DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      pt : Term := perm*t;

    begin
      Add(res,pt);
      Clear(pt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(f);
    Clear(f);
    f := res;
  end Permute_Polynomial;

  procedure Permute_Polynomial
              ( perm : in Permutation;
                f : in out DoblDobl_Complex_Laurentials.Poly ) is

    use DoblDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      pt : Term := perm*t;

    begin
      Add(res,pt);
      Clear(pt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(f);
    Clear(f);
    f := res;
  end Permute_Polynomial;

  procedure Permute_Polynomial
              ( perm : in Permutation;
                f : in out QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      pt : Term := perm*t;

    begin
      Add(res,pt);
      Clear(pt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(f);
    Clear(f);
    f := res;
  end Permute_Polynomial;

  procedure Permute_Polynomial
              ( perm : in Permutation;
                f : in out QuadDobl_Complex_Laurentials.Poly ) is

    use QuadDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      pt : Term := perm*t;

    begin
      Add(res,pt);
      Clear(pt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(f);
    Clear(f);
    f := res;
  end Permute_Polynomial;

  procedure Permute_Polynomial_System
              ( perm : in Permutation;
                f : in out Standard_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in f'range loop
      Permute_Polynomial(perm,f(i));
    end loop;
  end Permute_Polynomial_System;

  procedure Permute_Polynomial_System
              ( perm : in Permutation;
                f : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in f'range loop
      Permute_Polynomial(perm,f(i));
    end loop;
  end Permute_Polynomial_System;

  procedure Permute_Polynomial_System
              ( perm : in Permutation;
                f : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in f'range loop
      Permute_Polynomial(perm,f(i));
    end loop;
  end Permute_Polynomial_System;

  procedure Permute_Polynomial_System
              ( perm : in Permutation;
                f : in out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in f'range loop
      Permute_Polynomial(perm,f(i));
    end loop;
  end Permute_Polynomial_System;

  procedure Permute_Polynomial_System
              ( perm : in Permutation;
                f : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in f'range loop
      Permute_Polynomial(perm,f(i));
    end loop;
  end Permute_Polynomial_System;

  procedure Permute_Polynomial_System
              ( perm : in Permutation;
                f : in out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in f'range loop
      Permute_Polynomial(perm,f(i));
    end loop;
  end Permute_Polynomial_System;

-- MANIPULATION OF SYMBOLS TO SHIFT ADDED VARIABLES TO END :

  function Is_Match ( s1,s2 : string ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all initial characters in both strings match.

  begin
    if s1'last < s2'last then
      return false;
    else
      for i in s1'range loop
        exit when (i > s2'last);
        if s1(i) /= s2(i)
         then return false;
        end if;
      end loop;
    end if;
    return true;
  end Is_Match;

  function First_Embed_Symbol ( n : natural32; s : string ) return natural32 is

  -- DESCRIPTION :
  --   Returns the place in the symbol table of the first symbol whose
  --   initial characters match the given string.
  --   This is extremely important because the program assumes that
  --   these added variables come last and this is not necessarily so
  --   when the embedded system is read from file.
  --   Returns n+1 if there is not a matching symbol.

  begin
    for i in 1..n loop
      declare 
        sb : constant Symbol := Symbol_Table.Get(i);
      begin
        if Is_Match(sb,s)
         then return i;
        end if;
      end;
    end loop;
    return n+1;
  end First_Embed_Symbol;

  function Count_Embed_Symbols
             ( sbt : Array_of_Symbols; s : string ) return natural32 is

    cnt : natural32 := 0;

  begin
    for i in sbt'range loop
      if Is_Match(sbt(i),s)
       then cnt := cnt+1;
      end if;
    end loop;
    return cnt;
  end Count_Embed_Symbols;

  function Count_Embed_Symbols
             ( n : natural32; s : string ) return natural32 is

    cnt : natural32 := 0;

  begin
    for i in 1..n loop
      declare
        sb : constant Symbol := Symbol_Table.Get(i);
      begin
        if Is_Match(sb,s)
         then cnt := cnt+1;
        end if;
      end;
    end loop;
    return cnt;
  end Count_Embed_Symbols;

  function Swap ( t : Standard_Complex_Polynomials.Term; i,j : natural32 ) 
                return Standard_Complex_Polynomials.Term is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the term.

    use Standard_Complex_Polynomials;
    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
    res.dg(integer32(i)) := t.dg(integer32(j));
    res.dg(integer32(j)) := t.dg(integer32(i));
    return res;
  end Swap;

  function Swap ( t : Standard_Complex_Laurentials.Term; i,j : natural32 ) 
                return Standard_Complex_Laurentials.Term is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the term.

    use Standard_Complex_Laurentials;
    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
    res.dg(integer32(i)) := t.dg(integer32(j));
    res.dg(integer32(j)) := t.dg(integer32(i));
    return res;
  end Swap;

  function Swap ( t : DoblDobl_Complex_Polynomials.Term; i,j : natural32 ) 
                return DoblDobl_Complex_Polynomials.Term is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the term.

    use DoblDobl_Complex_Polynomials;
    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
    res.dg(integer32(i)) := t.dg(integer32(j));
    res.dg(integer32(j)) := t.dg(integer32(i));
    return res;
  end Swap;

  function Swap ( t : DoblDobl_Complex_Laurentials.Term; i,j : natural32 ) 
                return DoblDobl_Complex_Laurentials.Term is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the term.

    use DoblDobl_Complex_Laurentials;
    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
    res.dg(integer32(i)) := t.dg(integer32(j));
    res.dg(integer32(j)) := t.dg(integer32(i));
    return res;
  end Swap;

  function Swap ( t : QuadDobl_Complex_Polynomials.Term; i,j : natural32 ) 
                return QuadDobl_Complex_Polynomials.Term is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the term.

    use QuadDobl_Complex_Polynomials;
    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
    res.dg(integer32(i)) := t.dg(integer32(j));
    res.dg(integer32(j)) := t.dg(integer32(i));
    return res;
  end Swap;

  function Swap ( t : QuadDobl_Complex_Laurentials.Term; i,j : natural32 ) 
                return QuadDobl_Complex_Laurentials.Term is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the term.

    use QuadDobl_Complex_Laurentials;
    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
    res.dg(integer32(i)) := t.dg(integer32(j));
    res.dg(integer32(j)) := t.dg(integer32(i));
    return res;
  end Swap;

  procedure Swap ( p : in out Standard_Complex_Polynomials.Poly;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the polynomial p.

    use Standard_Complex_Polynomials;
    swp : Poly := Null_Poly;

    procedure Swap_in_Term ( t : in Term; cont : out boolean ) is

      nt : constant Term := Swap(t,i,j);

    begin
      Add(swp,nt);
      cont := true;
    end Swap_in_Term;
    procedure Swap_in_Terms is new Visiting_Iterator(Swap_in_Term);

  begin
    Swap_in_Terms(p);
    Clear(p);
    p := swp;
  end Swap;

  procedure Swap ( p : in out Standard_Complex_Laurentials.Poly;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the polynomial p.

    use Standard_Complex_Laurentials;
    swp : Poly := Null_Poly;

    procedure Swap_in_Term ( t : in Term; cont : out boolean ) is

      nt : constant Term := Swap(t,i,j);

    begin
      Add(swp,nt);
      cont := true;
    end Swap_in_Term;
    procedure Swap_in_Terms is new Visiting_Iterator(Swap_in_Term);

  begin
    Swap_in_Terms(p);
    Clear(p);
    p := swp;
  end Swap;

  procedure Swap ( p : in out DoblDobl_Complex_Polynomials.Poly;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the polynomial p.

    use DoblDobl_Complex_Polynomials;
    swp : Poly := Null_Poly;

    procedure Swap_in_Term ( t : in Term; cont : out boolean ) is

      nt : constant Term := Swap(t,i,j);

    begin
      Add(swp,nt);
      cont := true;
    end Swap_in_Term;
    procedure Swap_in_Terms is new Visiting_Iterator(Swap_in_Term);

  begin
    Swap_in_Terms(p);
    Clear(p);
    p := swp;
  end Swap;

  procedure Swap ( p : in out DoblDobl_Complex_Laurentials.Poly;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the polynomial p.

    use DoblDobl_Complex_Laurentials;
    swp : Poly := Null_Poly;

    procedure Swap_in_Term ( t : in Term; cont : out boolean ) is

      nt : constant Term := Swap(t,i,j);

    begin
      Add(swp,nt);
      cont := true;
    end Swap_in_Term;
    procedure Swap_in_Terms is new Visiting_Iterator(Swap_in_Term);

  begin
    Swap_in_Terms(p);
    Clear(p);
    p := swp;
  end Swap;

  procedure Swap ( p : in out QuadDobl_Complex_Polynomials.Poly;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the polynomial p.

    use QuadDobl_Complex_Polynomials;
    swp : Poly := Null_Poly;

    procedure Swap_in_Term ( t : in Term; cont : out boolean ) is

      nt : constant Term := Swap(t,i,j);

    begin
      Add(swp,nt);
      cont := true;
    end Swap_in_Term;
    procedure Swap_in_Terms is new Visiting_Iterator(Swap_in_Term);

  begin
    Swap_in_Terms(p);
    Clear(p);
    p := swp;
  end Swap;

  procedure Swap ( p : in out QuadDobl_Complex_Laurentials.Poly;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in the polynomial p.

    use QuadDobl_Complex_Laurentials;
    swp : Poly := Null_Poly;

    procedure Swap_in_Term ( t : in Term; cont : out boolean ) is

      nt : constant Term := Swap(t,i,j);

    begin
      Add(swp,nt);
      cont := true;
    end Swap_in_Term;
    procedure Swap_in_Terms is new Visiting_Iterator(Swap_in_Term);

  begin
    Swap_in_Terms(p);
    Clear(p);
    p := swp;
  end Swap;

  procedure Swap ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in every polynomial in the system p.

  begin
    for k in p'range loop
      Swap(p(k),i,j);
    end loop;
  end Swap;

  procedure Swap ( p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in every polynomial in the system p.

  begin
    for k in p'range loop
      Swap(p(k),i,j);
    end loop;
  end Swap;

  procedure Swap ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in every polynomial in the system p.

  begin
    for k in p'range loop
      Swap(p(k),i,j);
    end loop;
  end Swap;

  procedure Swap ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in every polynomial in the system p.

  begin
    for k in p'range loop
      Swap(p(k),i,j);
    end loop;
  end Swap;

  procedure Swap ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in every polynomial in the system p.

  begin
    for k in p'range loop
      Swap(p(k),i,j);
    end loop;
  end Swap;

  procedure Swap ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the variables i and j in every polynomial in the system p.

  begin
    for k in p'range loop
      Swap(p(k),i,j);
    end loop;
  end Swap;

  procedure Swap ( s : in out Standard_Complex_Solutions.Solution;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the values of i and j in the solution vector.

    use Standard_Complex_Numbers;
    tmp : constant Complex_Number := s.v(integer32(i));

  begin
    s.v(integer32(i)) := s.v(integer32(j)); 
    s.v(integer32(j)) := tmp;
  end Swap;

  procedure Swap ( s : in out DoblDobl_Complex_Solutions.Solution;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the values of i and j in the solution vector.

    use DoblDobl_Complex_Numbers;
    tmp : constant Complex_Number := s.v(integer32(i));

  begin
    s.v(integer32(i)) := s.v(integer32(j)); 
    s.v(integer32(j)) := tmp;
  end Swap;

  procedure Swap ( s : in out QuadDobl_Complex_Solutions.Solution;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the values of i and j in the solution vector.

    use QuadDobl_Complex_Numbers;
    tmp : constant Complex_Number := s.v(integer32(i));

  begin
    s.v(integer32(i)) := s.v(integer32(j)); 
    s.v(integer32(j)) := tmp;
  end Swap;

  procedure Swap ( sols : in out Standard_Complex_Solutions.Solution_List;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the values of i and j in every solution vector of the list.

    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Swap(ls.all,i,j);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Swap;

  procedure Swap ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the values of i and j in every solution vector of the list.

    use DoblDobl_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Swap(ls.all,i,j);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Swap;

  procedure Swap ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                   i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the values of i and j in every solution vector of the list.

    use QuadDobl_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Swap(ls.all,i,j);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Swap;

  procedure Swap_Symbols ( i,j : in natural32 ) is

  -- DESCRIPTION :
  --   Swaps the symbols for i and j in the symbol table.

    sbi : constant Symbol := Symbol_Table.Get(i);
    sbj : constant Symbol := Symbol_Table.Get(j);

  begin
    Symbol_Table.Replace(i,sbj);
    Symbol_Table.Replace(j,sbi);
  end Swap_Symbols;

  function Is_Added_Term
                ( t : Standard_Complex_Polynomials.Term; k : natural32 )
                return boolean is

  -- DESCRIPTION :
  --   Returns true all exponents in t.dg are zero, except the k last ones.

    use Standard_Complex_Polynomials;

  begin
    for i in t.dg'first..t.dg'last-integer32(k) loop
      if t.dg(i) /= 0
       then return false;                  -- we have an original term
      end if;
    end loop;
    for i in t.dg'last-integer32(k)+1..t.dg'last loop
      if t.dg(i) /= 0
       then return true;                   -- the term was of form c*z
      end if;
    end loop;
    return false;                          -- we found the constant term
  end Is_Added_Term;

  function Remove_Last_Variables
              ( p : Standard_Complex_Polynomials.Poly; k : natural32 )
              return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Removes the last k variables from the polynomial p.
  --   Also removes the terms c*z, for every added variable z.

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Remove_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      if not Is_Added_Term(t,k) then
        rt.cf := t.cf;
        rt.dg := new Standard_Natural_Vectors.Vector'
                       (t.dg(t.dg'first..t.dg'last-integer32(k)));
        Add(res,rt);
        Clear(rt);
      end if;
      cont := true;
    end Remove_in_Term;
    procedure Remove_in_Terms is new Visiting_Iterator(Remove_in_Term);

  begin
    Remove_in_Terms(p);
    return res;
  end Remove_Last_Variables;

  function Remove_Slices
              ( p : Standard_Complex_Poly_Systems.Poly_Sys; k : natural32 )
              return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns a system with the k slices and k last variables removed.

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last-integer32(k));

  begin
    for i in res'range loop
      res(i) := Remove_Last_Variables(p(i),k);
    end loop;
    return res;
  end Remove_Slices;

-- TARGET ROUTINES :

  procedure Add_Slack_Symbols ( k : in natural32 ) is

  -- DESCRIPTION :
  --   Adds k new symbols of the form "ssk", with k as index.

    sb : Symbol;

  begin
    Symbol_Table.Enlarge(k);
    for i in 1..k loop
      sb := (sb'range => ' ');
      sb(1..2) := "ss";
      declare
        stri : constant string := Convert(integer32(i));
      begin
        for j in stri'range loop
          sb(j+2) := stri(j);
        end loop;
      end;
     -- sb(3) := Convert_Decimal(i);
      Symbol_Table.Add(sb);
    end loop;
  end Add_Slack_Symbols;

  procedure Add_Embed_Symbols ( k : in natural32 ) is

  -- DESCRIPTION :
  --   Adds k new symbols of the form "zzk", with k as index.

    sb : Symbol;

  begin
    Symbol_Table.Enlarge(k);
    for i in 1..k loop
      sb := (sb'range => ' ');
      sb(1..2) := "zz";
      declare
        stri : constant string := Convert(integer32(i));
      begin
        for j in stri'range loop
          sb(j+2) := stri(j);
        end loop;
      end;
     -- sb(3) := Convert_Decimal(i);
      Symbol_Table.Add(sb);
    end loop;
  end Add_Embed_Symbols;

  procedure Read_and_Add_Symbol ( i : in natural32 ) is

  -- DESCRIPTION :
  --   Reads a symbol to denote the i-th variable.

    sb : Symbol;

  begin
    put("Give symbol for variable "); put(i,1); put(" : ");
    sb := (sb'range => ' ');
    Symbol_Table_io.get(sb);
    Symbol_Table.add(sb);
  end Read_and_Add_Symbol;

  procedure Add_Extra_Symbols ( k : in natural32 ) is

    ns : constant natural32 := Symbol_Table.Number;

  begin
    put("The current symbols are : ");
    Write_Symbol_Table(Standard_Output);
    Symbol_Table.Enlarge(k);
    if k = 1 then
      Read_and_Add_Symbol(ns+1);
    else
      put("Reading "); put(k,1);
      put_line(" extra symbols ...");
      for i in 1..k loop
        Read_and_Add_Symbol(ns+i);
      end loop;
    end if;
  end Add_Extra_Symbols;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out Standard_Complex_Laur_Systems.Laur_Sys ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap(sols,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap(sols,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap(sols,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap(sols,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap(sols,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  procedure Swap_Symbols_to_End
              ( n,k : in natural32; s : in string;
                p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    addsym : natural32 := First_Embed_Symbol(n,s);
    last : natural32 := n;

  begin
    while addsym /= n+1-k loop
      Swap(p,addsym,last);
      Swap(sols,addsym,last);
      Swap_Symbols(addsym,last);
      addsym := First_Embed_Symbol(last,s);
      last := last-1;
    end loop;
  end Swap_Symbols_to_End;

  function Permutation_for_Embed_Symbols
             ( nv,n,d : natural32 ) return Permutation is

  -- DESCRIPTION :
  --   Returns the permutation to sort the embed symbols.

    res : Permutation(1..integer32(nv));

  begin
    for i in 1..integer32(n) loop
      res(i) := i;
    end loop;
    for i in n+1..n+d loop
      declare
        s : constant Symbol := Symbol_Table.Get(i);
        k : constant natural32 := Convert(s(3..s'last));
      begin
        res(integer32(i)) := integer32(n+k);
      end;
    end loop;
    for i in integer32(n+d+1)..integer32(nv) loop
      res(i) := i;
    end loop;
    return res;
  end Permutation_for_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out Standard_Complex_Laur_Systems.Laur_Sys ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Solutions(p,sols);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Solutions(p,sols);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Solutions(p,sols);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Solutions(p,sols);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Solutions(p,sols);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Sort_Embed_Symbols
              ( nv,n,d : in natural32;
                f : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    p : constant Permutation(1..integer32(nv))
      := Permutation_for_Embed_Symbols(nv,n,d);

  begin
    Permute_Symbol_Table(p);
    Permute_Solutions(p,sols);
    Permute_Polynomial_System(p,f);
  end Sort_Embed_Symbols;

  procedure Standard_Read_Embedding
              ( file : in file_type;
                lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    use Standard_Complex_Poly_Systems;

    n : natural32 := 0;
    found : boolean;

  begin
    get(file,n); skip_line(file); -- to skip the #variables
    lp := new Poly_Sys(1..integer32(n));
    Symbol_Table.Init(n);
    get(file,lp.all);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    if found then
      get(file,sols);
    else 
      new_line;
      Read(sols);
    end if;
    dim := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,dim,"zz",lp.all,sols);
    if dim > 1
     then Sort_Embed_Symbols(n,n-dim,dim,lp.all,sols);
    end if;
  end Standard_Read_Embedding;

  procedure Standard_Read_Embedding
              ( file : in file_type;
                lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    use Standard_Complex_Laur_Systems;

    n : natural32 := 0;
    found : boolean;

  begin
    get(file,n); skip_line(file); -- to skip the #variables
    lp := new Laur_Sys(1..integer32(n));
    Symbol_Table.Init(n);
    get(file,lp.all);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    if found then
      get(file,sols);
    else 
      new_line;
      Read(sols);
    end if;
    dim := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,dim,"zz",lp.all,sols);
    if dim > 1
     then Sort_Embed_Symbols(n,n-dim,dim,lp.all,sols);
    end if;
  end Standard_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( file : in file_type;
                lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    use DoblDobl_Complex_Poly_Systems;

    n : natural32 := 0;
    found : boolean;

  begin
    get(file,n); skip_line(file); -- to skip the #variables
    lp := new Poly_Sys(1..integer32(n));
    Symbol_Table.Init(n);
    get(file,lp.all);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    if found then
      get(file,sols);
    else 
      new_line;
      Read(sols);
    end if;
    dim := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,dim,"zz",lp.all,sols);
    if dim > 1
     then Sort_Embed_Symbols(n,n-dim,dim,lp.all,sols);
    end if;
  end DoblDobl_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( file : in file_type;
                lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    use DoblDobl_Complex_Laur_Systems;

    n : natural32 := 0;
    found : boolean;

  begin
   -- get(file,n); skip_line(file); -- to skip the #variables
   -- lp := new Laur_Sys(1..integer32(n));
   -- Symbol_Table.Init(n);
    get(file,lp);
    n := natural32(lp'last);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    if found then
      get(file,sols);
    else 
      new_line;
      Read(sols);
    end if;
    dim := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,dim,"zz",lp.all,sols);
    if dim > 1
     then Sort_Embed_Symbols(n,n-dim,dim,lp.all,sols);
    end if;
  end DoblDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( file : in file_type;
                lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    use QuadDobl_Complex_Poly_Systems;

    n : natural32 := 0;
    found : boolean;

  begin
    get(file,n); skip_line(file); -- to skip the #variables
    lp := new Poly_Sys(1..integer32(n));
    Symbol_Table.Init(n);
    get(file,lp.all);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    if found then
      get(file,sols);
    else 
      new_line;
      Read(sols);
    end if;
    dim := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,dim,"zz",lp.all,sols);
    if dim > 1
     then Sort_Embed_Symbols(n,n-dim,dim,lp.all,sols);
    end if;
  end QuadDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( file : in file_type;
                lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    use QuadDobl_Complex_Laur_Systems;

    n : natural32 := 0;
    found : boolean;

  begin
   -- get(file,n); skip_line(file); -- to skip the #variables
   -- lp := new Laur_Sys(1..integer32(n));
   -- Symbol_Table.Init(n);
    get(file,lp);
    n := natural32(lp'last);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    if found then
      get(file,sols);
    else 
      new_line;
      Read(sols);
    end if;
    dim := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,dim,"zz",lp.all,sols);
    if dim > 1
     then Sort_Embed_Symbols(n,n-dim,dim,lp.all,sols);
    end if;
  end QuadDobl_Read_Embedding;

  procedure Standard_Read_Embedding
              ( file : in file_type;
                lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is

    n : natural32;

  begin
    Standard_Read_Embedding(file,lp,sols,dim);
    n := natural32(lp'last);
    nsl := Count_Embed_Symbols(n,"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(n-dim,nsl,"ss",lp.all,sols);
      if nsl > 1
       then Sort_Embed_Symbols(n,n-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end Standard_Read_Embedding;

  procedure Standard_Read_Embedding
              ( file : in file_type;
                lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is

    n : natural32;

  begin
    Standard_Read_Embedding(file,lp,sols,dim);
    n := natural32(lp'last);
    nsl := Count_Embed_Symbols(n,"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(n-dim,nsl,"ss",lp.all,sols);
      if nsl > 1
       then Sort_Embed_Symbols(n,n-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end Standard_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( file : in file_type;
                lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is

    n : natural32;

  begin
    DoblDobl_Read_Embedding(file,lp,sols,dim);
    n := natural32(lp'last);
    nsl := Count_Embed_Symbols(n,"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(n-dim,nsl,"ss",lp.all,sols);
      if nsl > 1
       then Sort_Embed_Symbols(n,n-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end DoblDobl_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( file : in file_type;
                lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is

    n : natural32;

  begin
    DoblDobl_Read_Embedding(file,lp,sols,dim);
    n := natural32(lp'last);
    nsl := Count_Embed_Symbols(n,"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(n-dim,nsl,"ss",lp.all,sols);
      if nsl > 1
       then Sort_Embed_Symbols(n,n-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end DoblDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( file : in file_type;
                lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is

    n : natural32;

  begin
    QuadDobl_Read_Embedding(file,lp,sols,dim);
    n := natural32(lp'last);
    nsl := Count_Embed_Symbols(n,"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(n-dim,nsl,"ss",lp.all,sols);
      if nsl > 1
       then Sort_Embed_Symbols(n,n-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end QuadDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( file : in file_type;
                lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is

    n : natural32;

  begin
    QuadDobl_Read_Embedding(file,lp,sols,dim);
    n := natural32(lp'last);
    nsl := Count_Embed_Symbols(n,"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(n-dim,nsl,"ss",lp.all,sols);
      if nsl > 1
       then Sort_Embed_Symbols(n,n-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end QuadDobl_Read_Embedding;

  procedure Standard_Read_Embedding
              ( lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    infile : file_type;

  begin
    new_line;
    put_line("Reading the name of the file for the embedded system.");
    Read_Name_and_Open_File(infile);
    Standard_Read_Embedding(infile,lp,sols,dim);
    close(infile);
  end Standard_Read_Embedding;

  procedure Standard_Read_Embedding
              ( lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    infile : file_type;

  begin
    new_line;
    put_line("Reading the name of the file for the embedded system.");
    Read_Name_and_Open_File(infile);
    Standard_Read_Embedding(infile,lp,sols,dim);
    close(infile);
  end Standard_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    infile : file_type;

  begin
    new_line;
    put_line("Reading the name of the file for the embedded system.");
    Read_Name_and_Open_File(infile);
    DoblDobl_Read_Embedding(infile,lp,sols,dim);
    close(infile);
  end DoblDobl_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    infile : file_type;

  begin
    new_line;
    put_line("Reading the name of the file for the embedded system.");
    Read_Name_and_Open_File(infile);
    DoblDobl_Read_Embedding(infile,lp,sols,dim);
    close(infile);
  end DoblDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    infile : file_type;

  begin
    new_line;
    put_line("Reading the name of the file for the embedded system.");
    Read_Name_and_Open_File(infile);
    QuadDobl_Read_Embedding(infile,lp,sols,dim);
    close(infile);
  end QuadDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    infile : file_type;

  begin
    new_line;
    put_line("Reading the name of the file for the embedded system.");
    Read_Name_and_Open_File(infile);
    QuadDobl_Read_Embedding(infile,lp,sols,dim);
    close(infile);
  end QuadDobl_Read_Embedding;

  procedure Standard_Read_Embedding
              ( lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is
  begin
    Standard_Read_Embedding(lp,sols,dim);
    nsl := Count_Embed_Symbols(natural32(lp'last),"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(natural32(lp'last)-dim,nsl,"ss",lp.all,sols);
      if nsl > 1 then
        Sort_Embed_Symbols
          (natural32(lp'last),natural32(lp'last)-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end Standard_Read_Embedding;

  procedure Standard_Read_Embedding
              ( lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is
  begin
    Standard_Read_Embedding(lp,sols,dim);
    nsl := Count_Embed_Symbols(natural32(lp'last),"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(natural32(lp'last)-dim,nsl,"ss",lp.all,sols);
      if nsl > 1 then
        Sort_Embed_Symbols
          (natural32(lp'last),natural32(lp'last)-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end Standard_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is
  begin
    DoblDobl_Read_Embedding(lp,sols,dim);
    nsl := Count_Embed_Symbols(natural32(lp'last),"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(natural32(lp'last)-dim,nsl,"ss",lp.all,sols);
      if nsl > 1 then
        Sort_Embed_Symbols
          (natural32(lp'last),natural32(lp'last)-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end DoblDobl_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is
  begin
    DoblDobl_Read_Embedding(lp,sols,dim);
    nsl := Count_Embed_Symbols(natural32(lp'last),"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(natural32(lp'last)-dim,nsl,"ss",lp.all,sols);
      if nsl > 1 then
        Sort_Embed_Symbols
          (natural32(lp'last),natural32(lp'last)-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end DoblDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is
  begin
    QuadDobl_Read_Embedding(lp,sols,dim);
    nsl := Count_Embed_Symbols(natural32(lp'last),"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(natural32(lp'last)-dim,nsl,"ss",lp.all,sols);
      if nsl > 1 then
        Sort_Embed_Symbols
          (natural32(lp'last),natural32(lp'last)-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end QuadDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim,nsl : out natural32 ) is
  begin
    QuadDobl_Read_Embedding(lp,sols,dim);
    nsl := Count_Embed_Symbols(natural32(lp'last),"ss");
    if nsl > 0 then
      Swap_Symbols_to_End(natural32(lp'last)-dim,nsl,"ss",lp.all,sols);
      if nsl > 1 then
        Sort_Embed_Symbols
          (natural32(lp'last),natural32(lp'last)-dim-nsl,nsl,lp.all,sols);
      end if;
    end if;
  end QuadDobl_Read_Embedding;

  procedure Standard_Read_Embedding
              ( name : in string;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    Open(file,in_file,name);
    Standard_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put(name); put_line(".");
      put("Or incorrect formats of system/solutions in ");
      put(name); put_line(".");
      Standard_Read_Embedding(lp,sols,dim);
  end Standard_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( name : in string;
                lp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    Open(file,in_file,name);
    DoblDobl_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put(name); put_line(".");
      put("Or incorrect formats of system/solutions in ");
      put(name); put_line(".");
      DoblDobl_Read_Embedding(lp,sols,dim);
  end DoblDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( name : in string;
                lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    Open(file,in_file,name);
    QuadDobl_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put(name); put_line(".");
      put("Or incorrect formats of system/solutions in ");
      put(name); put_line(".");
      QuadDobl_Read_Embedding(lp,sols,dim);
  end QuadDobl_Read_Embedding;

  procedure Standard_Read_Embedding
              ( name : in string;
                lp : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    Open(file,in_file,name);
    Standard_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put(name); put_line(".");
      put("Or incorrect formats of system/solutions in ");
      put(name); put_line(".");
      Standard_Read_Embedding(lp,sols,dim);
  end Standard_Read_Embedding;

  procedure DoblDobl_Read_Embedding
              ( name : in string;
                lp : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    Open(file,in_file,name);
    DoblDobl_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put(name); put_line(".");
      put("Or incorrect formats of system/solutions in ");
      put(name); put_line(".");
      DoblDobl_Read_Embedding(lp,sols,dim);
  end DoblDobl_Read_Embedding;

  procedure QuadDobl_Read_Embedding
              ( name : in string;
                lp : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    Open(file,in_file,name);
    QuadDobl_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put(name); put_line(".");
      put("Or incorrect formats of system/solutions in ");
      put(name); put_line(".");
      QuadDobl_Read_Embedding(lp,sols,dim);
  end QuadDobl_Read_Embedding;

  procedure Get_Multprec_System 
              ( stsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                mpsys : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                size,dim : in natural32 ) is

    res : Multprec_Complex_Poly_Systems.Poly_Sys
            (stsys'first..stsys'last-integer32(dim));
    op : Standard_Complex_Poly_Systems.Poly_Sys
            (stsys'first..stsys'last-integer32(dim));
    ans : character;
    infile : file_type;
    m : natural32 := 0;
    n : constant natural32 := natural32(stsys'last);
    newres : Multprec_Complex_Polynomials.Poly;

  begin
    put("Do you wish to read in the multi-precision coefficients? (y/n) "); 
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Reading the polynomial system"
               & " with multi-precision coefficients.");
      Read_Name_and_Open_File(infile);
      get(infile,m,res);
      for i in res'range loop     -- patch!!! #vars = size of symbol table
        for j in 1..dim loop
          newres := Remove_Variable(res(i),integer32(n+1-j));
          Multprec_Complex_Polynomials.Copy(newres,res(i));
          Multprec_Complex_Polynomials.Clear(newres);
        end loop;
      end loop;
    else
      op := Remove_Slices(stsys,dim);
      res := Convert(op);
    end if;
    Set_Size(res,size);
    mpsys := new Multprec_Complex_Poly_Systems.Poly_Sys'(res);
  end Get_Multprec_System;

  procedure Determine_Order
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
 
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    perm : Standard_Integer_Vectors.Vector(1..integer32(nbunk));
    pp : Poly_Sys(p'range);
    ans : character;
 
  begin
    put("There are "); put(nbunk,1); put(" unknowns, ordered as ");
    Write_Symbol_Table(Standard_Output);
    put("Do you wish to determine a new order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      perm := Interactive_Read_Permutation(nbunk);
      pp := p*Permutation(perm);
      Permute_Symbol_Table(Permutation(perm));
      Clear(p); p := pp;
    end if;
  end Determine_Order;

  procedure Determine_Order
               ( p : in out Standard_Complex_Laur_Systems.Laur_Sys ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;
 
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    perm : Standard_Integer_Vectors.Vector(1..integer32(nbunk));
    pp : Laur_Sys(p'range);
    ans : character;
 
  begin
    put("There are "); put(nbunk,1); put(" unknowns, ordered as ");
    Write_Symbol_Table(Standard_Output);
    put("Do you wish to determine a new order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      perm := Interactive_Read_Permutation(nbunk);
      pp := p*Permutation(perm);
      Permute_Symbol_Table(Permutation(perm));
      Clear(p); p := pp;
    end if;
  end Determine_Order;

  procedure Determine_Order
               ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
 
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    perm : Standard_Integer_Vectors.Vector(1..integer32(nbunk));
    pp : Poly_Sys(p'range);
    ans : character;
 
  begin
    put("There are "); put(nbunk,1); put(" unknowns, ordered as ");
    Write_Symbol_Table(Standard_Output);
    put("Do you wish to determine a new order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      perm := Interactive_Read_Permutation(nbunk);
      pp := p*Permutation(perm);
      Permute_Symbol_Table(Permutation(perm));
      Clear(p); p := pp;
    end if;
  end Determine_Order;

  procedure Determine_Order
               ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;
 
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    perm : Standard_Integer_Vectors.Vector(1..integer32(nbunk));
    pp : Laur_Sys(p'range);
    ans : character;
 
  begin
    put("There are "); put(nbunk,1); put(" unknowns, ordered as ");
    Write_Symbol_Table(Standard_Output);
    put("Do you wish to determine a new order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      perm := Interactive_Read_Permutation(nbunk);
      pp := p*Permutation(perm);
      Permute_Symbol_Table(Permutation(perm));
      Clear(p); p := pp;
    end if;
  end Determine_Order;

  procedure Determine_Order
               ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
 
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    perm : Standard_Integer_Vectors.Vector(1..integer32(nbunk));
    pp : Poly_Sys(p'range);
    ans : character;
 
  begin
    put("There are "); put(nbunk,1); put(" unknowns, ordered as ");
    Write_Symbol_Table(Standard_Output);
    put("Do you wish to determine a new order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      perm := Interactive_Read_Permutation(nbunk);
      pp := p*Permutation(perm);
      Permute_Symbol_Table(Permutation(perm));
      Clear(p); p := pp;
    end if;
  end Determine_Order;

  procedure Determine_Order
               ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;
 
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    perm : Standard_Integer_Vectors.Vector(1..integer32(nbunk));
    pp : Laur_Sys(p'range);
    ans : character;
 
  begin
    put("There are "); put(nbunk,1); put(" unknowns, ordered as ");
    Write_Symbol_Table(Standard_Output);
    put("Do you wish to determine a new order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      perm := Interactive_Read_Permutation(nbunk);
      pp := p*Permutation(perm);
      Permute_Symbol_Table(Permutation(perm));
      Clear(p); p := pp;
    end if;
  end Determine_Order;

  procedure Determine_Order
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
 
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    perm : Standard_Integer_Vectors.Vector(1..integer32(nbunk));
    pp : Poly_Sys(p'range);
    ans : character;
 
  begin
    put("There are "); put(nbunk,1); put(" unknowns, ordered as ");
    Write_Symbol_Table(Standard_Output);
    put("Do you wish to determine a new order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      perm := Interactive_Read_Permutation(nbunk);       
      pp := p*Permutation(perm);
      Permute_Symbol_Table(Permutation(perm));
      Clear(p); p := pp;
      Permute_Solutions(Permutation(perm),sols);
    end if;
  end Determine_Order;

  function Square_and_Embed
              ( p : Standard_Complex_Poly_Systems.Poly_Sys; k : natural32 )
              return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    nbequ : constant natural32 := natural32(p'length);
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));
    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : Poly_Sys(1..integer32(max)) := Square(p);
    ep : Poly_Sys(1..integer32(max+k));

  begin
    if nbequ /= nbunk then
      if nbequ > nbunk
       then Add_Slack_Symbols(nbequ-nbunk);
      end if;
    end if;
    if k = 0 then
      return sp;
    else
      Add_Embed_Symbols(k);
      if nbunk > nbequ then
        for i in nbunk-k+1..nbunk loop
          exit when (i = nbequ);  -- do not wipe out original eqs
          Clear(sp(integer32(i)));
        end loop;
      end if;
      ep := Slice_and_Embed(sp,k);
      return ep;
    end if;
  end Square_and_Embed;

end Witness_Sets_io;
