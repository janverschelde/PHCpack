with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Binomial_Varieties;
with Standard_Binomial_Varieties_io;
with Affine_Binomial_Iterator;
with Standard_Affine_Binomials;
with Standard_Monomial_Maps_io;         use Standard_Monomial_Maps_io;
with Standard_Monomial_Map_Solvers;

package body Standard_Permanent_Factors is

  procedure Solve_Affine_Subsystem 
               ( output : in boolean; p : in Laur_Sys;
                 eqcnt,s0cnt : in integer32;
                 eq,s0,free : in Standard_Integer_Vectors.Vector ) is

    q : constant Laur_Sys(1..eqcnt)
      := Standard_Affine_Binomials.Subsystem(p,eqcnt,eq);
    q1 : Laur_Sys(q'range)
       := Standard_Affine_Binomials.Eliminate_Variables(q,s0,s0cnt);
    nbfree : constant integer32 := Standard_Integer_Vectors.Sum(free);
    sf : constant Standard_Integer_Vectors.Vector
       := Standard_Affine_Binomials.Eliminate_Variables(free,s0,s0cnt);
    q2 : Laur_Sys(q1'range)
       := Standard_Affine_Binomials.Eliminate_Variables(q1,sf,nbfree);
    fail : boolean;
    d : integer32;
    M : Standard_Integer_Matrices.Link_to_Matrix;
    cff,tmp : Solution_List;
    use Standard_Binomial_Varieties;
 
  begin
    put_line("the subsystem (after elimination of variables) : "); put(q2);
    if output
     then Black_Box_Solver(standard_output,q2,fail,d,M,cff);
     else Black_Box_Solver(q2,fail,d,M,cff);
    end if;
    if not Is_Null(cff) then
      tmp := cff;
      if s0cnt = 0 then
        for i in 1..Length_Of(cff) loop
          put("solution "); put(i,1); put_line(" : ");
          Standard_Binomial_Varieties_io.Write_Header
            (standard_output,natural32(M'last(1)),natural32(d));
          Standard_Binomial_Varieties_io.Write_Solution
            (standard_output,natural32(d),M.all,Head_Of(tmp).v);
          tmp := Tail_of(tmp);
        end loop;
      else
        for i in 1..Length_Of(cff) loop
          put("solution "); put(i,1); put_line(" : ");
          Standard_Binomial_Varieties_io.Write_Header
            (standard_output,natural32(s0'last),natural32(d+nbfree));
          Standard_Binomial_Varieties_io.Write_Affine_Solution
            (standard_output,natural32(d),s0,free,M.all,Head_Of(tmp).v);
          tmp := Tail_of(tmp);
        end loop;
      end if;
    end if;
    Clear(q1); Clear(q2);
  end Solve_Affine_Subsystem;

  function Monomial_Map_Solution
             ( output : in boolean; p : in Laur_Sys;
               eqcnt,s0cnt : in integer32;
               eq,s0,free : in Standard_Integer_Vectors.Vector )
             return Link_to_Monomial_Map_Array is

    res : Standard_Monomial_Maps.Link_to_Monomial_Map_Array;
    q : constant Laur_Sys(1..eqcnt)
      := Standard_Affine_Binomials.Subsystem(p,eqcnt,eq);
    q1 : Laur_Sys(q'range)
       := Standard_Affine_Binomials.Eliminate_Variables(q,s0,s0cnt);
    nbfree : constant integer32 := Standard_Integer_Vectors.Sum(free);
    sf : constant Standard_Integer_Vectors.Vector
       := Standard_Affine_Binomials.Eliminate_Variables(free,s0,s0cnt);
    q2 : Laur_Sys(q1'range)
       := Standard_Affine_Binomials.Eliminate_Variables(q1,sf,nbfree);
 
  begin
    if output then
      put_line("the subsystem (after elimination of variables) : ");
      put(q2);
    end if;
    if s0cnt = 0 then
      res := Standard_Monomial_Map_Solvers.Toric_Solve(q2); 
    elsif nbfree = 0 then
      res := Standard_Monomial_Map_Solvers.Affine_Solve(q2,s0); 
    else
      res := Standard_Monomial_Map_Solvers.Affine_Solve(q2,nbfree,s0,free); 
    end if;
    if res /= null and output then
      put_line("the solutions to the subsystem : ");
      put(res.all);
    end if;
    Clear(q1); Clear(q2);
    return res;
  end Monomial_Map_Solution;

  procedure Show_Selection
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 s0,s1 : in Standard_Integer_Vectors.Vector;
                 cnt,s0cnt,eqcnt : in integer32; fail : out boolean ) is

    nq : constant integer32 := A'last(1)/2;
    eq : Standard_Integer_Vectors.Vector(1..nq);
    eqcnt2,nbfree : integer32;
    free : Standard_Integer_Vectors.Vector(s0'range);
    ok : boolean;

  begin
    put(cnt,3); put(" : ");
    put(" s[1] : "); put(s1); 
    put(" #eq : "); put(eqcnt,1);
    put(" s[0] : "); put(s0); 
    put(" # : "); put(s0cnt,1);
    put(" D : "); put(s0'last - s0cnt - eqcnt,1); new_line;
    Standard_Affine_Binomials.Nonzero_Binomials(A,s0,eq,eqcnt2,ok);
    put("   -> eqs : "); put(eq); put(" #eq : "); put(eqcnt2,1);
    if not ok then
      put_line("  invalid s[0], BUG!");
    else
      put("  okay");
      free := Standard_Affine_Binomials.Free_Variables(A,s0);
      put("  free : "); put(free); new_line;
      if eqcnt2 > 0 then
        Solve_Affine_Subsystem(false,p,eqcnt2,s0cnt,eq,s0,free);
      else
        nbfree := Standard_Integer_Vectors.Sum(free);
        put_line("free solution :");
        Standard_Binomial_Varieties_io.Write_Header
          (standard_output,natural32(s0'last),natural32(nbfree));
        Standard_Binomial_Varieties_io.Write_Free_Affine_Solution
          (standard_output,s0,free);
      end if;
    end if;
    fail := not ok;
  end Show_Selection;

  procedure Append_Solution_Maps
               ( output : in boolean;
                 p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 s0,s1 : in Standard_Integer_Vectors.Vector;
                 cnt,s0cnt,eqcnt : in integer32; fail : out boolean;
                 first,last : in out Monomial_Map_List ) is

    nq : constant integer32 := A'last(1)/2;
    eq : Standard_Integer_Vectors.Vector(1..nq);
    eqcnt2,nbfree : integer32;
    free : Standard_Integer_Vectors.Vector(s0'range);
    sols : Link_to_Monomial_Map_Array;
    ok : boolean;

  begin
    if output then
      put(cnt,3); put(" : ");
      put(" s[1] : "); put(s1); 
      put(" #eq : "); put(eqcnt,1);
      put(" s[0] : "); put(s0); 
      put(" # : "); put(s0cnt,1);
      put(" D : "); put(s0'last - s0cnt - eqcnt,1); new_line;
    end if;
    Standard_Affine_Binomials.Nonzero_Binomials(A,s0,eq,eqcnt2,ok);
    if output then
      put("   -> eqs : "); put(eq); put(" #eq : "); put(eqcnt2,1);
      if not ok
       then put_line("  invalid s[0], BUG!");
       else put("  okay");
      end if;
    end if;
    if ok then
      free := Standard_Affine_Binomials.Free_Variables(A,s0);
      if output
       then put("  free : "); put(free); new_line;
      end if;
      if eqcnt2 > 0 then
        sols := Monomial_Map_Solution(output,p,eqcnt2,s0cnt,eq,s0,free);
        if sols /= null
         then Concatenate(sols.all,first,last); Clear(sols);
        end if;
      else
        nbfree := Standard_Integer_Vectors.Sum(free);
        declare
          map : constant Monomial_Map(free'last)
              := Standard_Monomial_Map_Solvers.Free_Monomial_Map
                   (free'last,nbfree,free);
        begin
          Append(first,last,map);
          if output
           then put_line("free solution : "); put(map);
          end if;
        end;
      end if;
    end if;
    fail := not ok;
  end Append_Solution_Maps;

  procedure Selection_Iterator
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 max : in integer32; puredim : in boolean ) is

    n : constant integer32 := A'last(2);
    nq : constant integer32 := A'last(1)/2;
    cnt : integer32 := 0;
    s0 : Standard_Integer_Vectors.Vector(1..n);
    s1 : Standard_Integer_Vectors.Vector(1..n);
    eqcnt : integer32;
    s0cnt,dimension : integer32;
    bug : boolean := false;

  begin
    put_line("enumeration of selections with iterator : ");
    Affine_Binomial_Iterator.Initialize_Iterator(A,max);
    loop
      Affine_Binomial_Iterator.Next_Selection(s0,s0cnt,s1,eqcnt);
      exit when (s0cnt < 0);
      if puredim then
        dimension := s0'last - s0cnt - eqcnt;
        if dimension = n - nq then
          cnt := cnt + 1;
          Show_Selection(p,A,s0,s1,cnt,s0cnt,eqcnt,bug);
        end if;
      else
        cnt := cnt + 1;
        Show_Selection(p,A,s0,s1,cnt,s0cnt,eqcnt,bug);
      end if;
      exit when bug;
    end loop;
    Affine_Binomial_Iterator.Clear_Iterator;
  end Selection_Iterator;

  procedure Selection_Iterator
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 max : in integer32; puredim,output : in boolean;
                 sols : out Monomial_Map_List ) is

    n : constant integer32 := A'last(2);
    nq : constant integer32 := A'last(1)/2;
    cnt : integer32 := 0;
    s0 : Standard_Integer_Vectors.Vector(1..n);
    s1 : Standard_Integer_Vectors.Vector(1..n);
    eqcnt : integer32;
    s0cnt,dimension : integer32;
    bug : boolean := false;
    last : Monomial_Map_List := sols;

  begin
    if output
     then put_line("enumeration of selections with iterator : ");
    end if;
    Affine_Binomial_Iterator.Initialize_Iterator(A,max);
    loop
      Affine_Binomial_Iterator.Next_Selection(s0,s0cnt,s1,eqcnt);
      exit when (s0cnt < 0);
      if puredim then
        dimension := s0'last - s0cnt - eqcnt;
        if dimension = n - nq then
          cnt := cnt + 1;
          Append_Solution_Maps(output,p,A,s0,s1,cnt,s0cnt,eqcnt,bug,sols,last);
        end if;
      else
        cnt := cnt + 1;
        Append_Solution_Maps(output,p,A,s0,s1,cnt,s0cnt,eqcnt,bug,sols,last);
      end if;
      exit when bug;
    end loop;
    Affine_Binomial_Iterator.Clear_Iterator;
  end Selection_Iterator;

  procedure Recursive_Enumerator
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 max : in integer32; puredim : in boolean ) is

    n : constant integer32 := A'last(2);
    nq : constant integer32 := A'last(1)/2;
    cnt : integer32 := 0;

    procedure Write ( s0 : in Standard_Integer_Vectors.Vector;
                      s0cnt : in integer32;
                      s1 : in Standard_Integer_Vectors.Vector;
                      eqcnt : in integer32;
                      continue : out boolean ) is

      dimension : constant integer32 := s0'last - s0cnt - eqcnt;
      bug : boolean := false;

    begin
      if not puredim then
        cnt := cnt + 1;
        Show_Selection(p,A,s0,s1,cnt,s0cnt,eqcnt,bug);
      else
        if dimension = n - nq then
          cnt := cnt + 1;
          Show_Selection(p,A,s0,s1,cnt,s0cnt,eqcnt,bug);
        end if;
      end if;
      continue := not bug;
    end Write;
    procedure Enum is 
      new Affine_Binomial_Iterator.Enumerate(Report=>Write);

  begin
    Enum(A,max);
  end Recursive_Enumerator;

  procedure Recursive_Enumerator
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 max : in integer32; puredim,output : in boolean;
                 sols : out Monomial_Map_List ) is

    n : constant integer32 := A'last(2);
    nq : constant integer32 := A'last(1)/2;
    cnt : integer32 := 0;
    last : Monomial_Map_List := sols;

    procedure Write ( s0 : in Standard_Integer_Vectors.Vector;
                      s0cnt : in integer32;
                      s1 : in Standard_Integer_Vectors.Vector;
                      eqcnt : in integer32;
                      continue : out boolean ) is

      dimension : constant integer32 := s0'last - s0cnt - eqcnt;
      bug : boolean := false;

    begin
      if not puredim then
        cnt := cnt + 1;
        Append_Solution_Maps
          (output,p,A,s0,s1,cnt,s0cnt,eqcnt,bug,sols,last);
      else
        if dimension = n - nq then
          cnt := cnt + 1;
          Append_Solution_Maps
            (output,p,A,s0,s1,cnt,s0cnt,eqcnt,bug,sols,last);
        end if;
      end if;
      continue := not bug;
    end Write;
    procedure Enum is 
      new Affine_Binomial_Iterator.Enumerate(Report=>Write);

  begin
    Enum(A,max);
  end Recursive_Enumerator;

  procedure Pruning_Maximum ( nv,nq : in integer32; max : out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for the new maximum on variables to be set
  --   to zero, after making a suggestion based on the number nv of variables 
  --   and the number nq of equations.

  begin
    max := nv;
    put("current max is "); put(max,1); 
    put_line(" (#variables)");
    if nv > nq then
      put("expected dimension : "); put(nv-nq,1);
      put_line(" (#variables - #equations)");
    end if;
    put("give new max : "); get(max); 
  end Pruning_Maximum;

  procedure Silent_Affine_Solutions_with_Recursion
              ( p : in Laur_Sys; sols : out Monomial_Map_List;
                fail : out boolean ) is

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    nq : constant integer32 := p'last;
    A : Standard_Integer_Matrices.Matrix(1..2*nq,1..n);
  
  begin
    Standard_Affine_Binomials.Incidence_Matrix(p,A,fail);
    if not fail
     then Recursive_Enumerator(p,A,n,false,false,sols);
    end if;
  end Silent_Affine_Solutions_with_Recursion;

  procedure Interactive_Affine_Solutions_with_Recursion
              ( p : in Laur_Sys; sols : out Monomial_Map_List;
                fail : out boolean ) is

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    nq : constant integer32 := p'last;
    A : Standard_Integer_Matrices.Matrix(1..2*p'last,1..n);
    puredim,output : boolean;
    s0_max : integer32 := n;
    ans : character;
  
  begin
    Standard_Affine_Binomials.Incidence_Matrix(p,A,fail);
   -- if fail then
   --   put_line("The system is not binomial!");
   -- else
    if not fail then
      put_line("The incidence matrix : "); put(A); new_line;
      Pruning_Maximum(n,nq,s0_max);
      put("pure dimensional ? (y/n) "); Ask_Yes_or_No(ans);
      puredim := (ans = 'y');
     -- put("write solution maps to file ? (y/n) "); Ask_Yes_or_No(ans);
     -- if ans = 'n' then
     --   Recursive_Enumerator(p,A,s0_max,puredim);
     -- else
        put("write solution maps to screen ? (y/n) "); Ask_Yes_or_No(ans);
        output := (ans = 'y');
        Recursive_Enumerator(p,A,s0_max,puredim,output,sols);
     -- end if;
    end if;
  end Interactive_Affine_Solutions_with_Recursion;

  procedure Silent_Affine_Solutions_with_Iterator
              ( p : in Laur_Sys; puretopdim : in boolean;
                sols : out Monomial_Map_List; fail : out boolean ) is

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    nq : constant integer32 := p'last;
    A : Standard_Integer_Matrices.Matrix(1..2*nq,1..n);
  
  begin
    Standard_Affine_Binomials.Incidence_Matrix(p,A,fail);
    if not fail
     then Selection_Iterator(p,A,n,puretopdim,false,sols);
    end if;
  end Silent_Affine_Solutions_with_Iterator;

  procedure Interactive_Affine_Solutions_with_Iterator
              ( p : in Laur_Sys; sols : out Monomial_Map_List;
                fail : out boolean ) is

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    nq : constant integer32 := p'last;
    A : Standard_Integer_Matrices.Matrix(1..2*p'last,1..n);
    puredim,output : boolean;
    s0_max : integer32 := n;
    ans : character;
  
  begin
    Standard_Affine_Binomials.Incidence_Matrix(p,A,fail);
   -- if fail then
   --   put_line("The system is not binomial!");
   -- else
    if not fail then
      put_line("The incidence matrix : "); put(A); new_line;
      Pruning_Maximum(n,nq,s0_max);
      put("pure dimensional ? (y/n) "); Ask_Yes_or_No(ans);
      puredim := (ans = 'y');
     -- put("write solution maps to file ? (y/n) "); Ask_Yes_or_No(ans);
     -- if ans = 'n' then
     --   Selection_Iterator(p,A,s0_max,puredim);
     -- else
        put("write solution maps to screen ? (y/n) "); Ask_Yes_or_No(ans);
        output := (ans = 'y');
        Selection_Iterator(p,A,s0_max,puredim,output,sols);
     -- end if;
    end if;
  end Interactive_Affine_Solutions_with_Iterator;

end Standard_Permanent_Factors;
