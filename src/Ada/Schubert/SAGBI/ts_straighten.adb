with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Polynomials_io;    use Standard_Bracket_Polynomials_io;
with Straightening_Syzygies;             use Straightening_Syzygies;

procedure ts_straighten is

-- DESCRIPTION :
--   Enumerates all brackets in the Pluecker embedding and generates
--   the van der Waerden syzygies.

  function Number_of_Brackets ( n,d : natural32 ) return natural32 is

  -- DESCIPTION :
  --   Returns the number of brackets of d entries chosen from n numbers.

    a,b : natural32;

  begin
    a := 1;
    for i in d+1..n loop
      a := a*i;
    end loop;
    b := 1;
    for i in 1..n-d loop
      b := b*i;
    end loop;
    return a/b;
  end Number_of_Brackets;

  function Number_of_Linear_Equations ( n,d : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of linear equations you need to determine a
  --   d-plane in C^n completely.

  begin
    return (n-d)*d;
  end Number_of_Linear_Equations;

  function Number_of_Zero_Brackets ( n,d : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximal number of brackets that can be set to zero.

  begin
    return (Number_of_Brackets(n,d) - Number_of_Linear_Equations(n,d) - 1);
  end Number_of_Zero_Brackets;

  function Create ( b1,b2 : Bracket ) return Bracket_Term is

    bm : Bracket_Monomial := Create(b1);
    bt : Bracket_Term;

  begin
    Multiply(bm,b2);
    bt.coeff := Create(1.0);
    bt.monom := bm;
    return bt;   
  end Create;

  procedure Enumerate ( n,d,k,start : in integer32; accu : in out Bracket;
                        cnt : in out integer32 ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of all brackets.

  begin
    if k > d then
      -- put(accu); new_line;
      cnt := cnt + 1;
    else
      for l in start..n loop
        accu(k) := natural32(l);
        Enumerate(n,d,k+1,l+1,accu,cnt);
      end loop;
    end if;
  end Enumerate;

  procedure Enumerate2 ( n,d,k,start : in integer32; b : in Bracket;
                         accu : in out Bracket;
                         cntstd,cntnonstd : in out integer32;
                         nonstd : in out Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of all brackets, with b < accu and with
  --   a test whether the pair b*accu forms a Standard monomial.

  -- ON ENTRY :
  --   n             total number of indices to choose from;
  --   d             number of indices in the brackets;
  --   k             current entry in accu that is to be filled;
  --   start         indices will be taken in between start..n;
  --   b             previously enumerated bracket;
  --   accu          accumulating parameter, filled in upto (k-1)th entry;
  --   cntstd        current number of quadratic standard monomials;
  --   cntnonstd     current number of quadratic nonstandard monomials;
  --   nonstd        current polynomial of quadratic nonstandard monomials.

  -- ON RETURN :
  --   accu          updated accumulating parameter, accu(k) is filled in;
  --   cnstd         updated number of quadratic standard monomials;
  --   cntnonstd     updated number of quadratic nonstandard monomials;
  --   nonstd        updated polynomial of quadratic nonstandard monomials.

    s : natural32;

  begin
    if k > d then
      if b < accu then
        -- put(b); put(accu); 
        s := Brackets.Is_Standard(b,accu);
        if s = 0 then
          -- put_line(" is a Standard monomial.");
          cntstd := cntstd + 1;
        else
          -- put(" is not a Standard monomial with s = ");
          -- put(s,1); new_line;
          cntnonstd := cntnonstd + 1;
          Add(nonstd,Create(b,accu));
        end if;
      end if;
    else
      for ell in start..n loop
        accu(k) := natural32(ell);
        Enumerate2(n,d,k+1,ell+1,b,accu,cntstd,cntnonstd,nonstd);
      end loop;
    end if;
  end Enumerate2;

  procedure Enumerate1 ( n,d,k,start : integer32; acc1,acc2 : in out Bracket;
                         cntstd,cntnonstd : in out integer32;
                         nonstd : in out Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of all brackets, with acc1 < acc2 and with
  --   a test whether the pair acc1*acc2 forms a Standard monomial.
  --   Counts the standard and nonstandard monomials and adds every
  --   nonStandard monomial to the polynomial nonstd.

  begin
    if k > d then
      Enumerate2(n,d,1,integer32(acc1(1)),acc1,acc2,cntstd,cntnonstd,nonstd);
    else
      for ell in start..n loop
        acc1(k) := natural32(ell);
        Enumerate1(n,d,k+1,ell+1,acc1,acc2,cntstd,cntnonstd,nonstd);
      end loop;
    end if;
  end Enumerate1;

  procedure Enumerate_Pairs ( n,d : in integer32;
                              nonstd : in out Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Enumerates all pairs (b1,b2), with b1 < b2 and checks whether
  --   the corresponding bracket monomial b1*b2 is standard or not.

    b1,b2 : Bracket(1..d);
    cntstd,cntnonstd : integer32 := 0;

  begin
    Enumerate1(n,d,1,1,b1,b2,cntstd,cntnonstd,nonstd);
    put("#Standard bi-monomials    : "); put(cntstd,1);    new_line;
    put("#nonStandard bi-monomials : "); put(cntnonstd,1); new_line;
  end Enumerate_Pairs;

  procedure Enumerate_Brackets ( n,d : in integer32 ) is

    acc : Bracket(1..d);
    cnt : integer32 := 0;
    nonstd : Bracket_Polynomial;

  begin
    Enumerate(n,d,1,1,acc,cnt);
    put("#brackets : "); put(cnt,1); new_line;
    Enumerate_Pairs(n,d,nonstd);
    put_line("The polynomial of all quadratic nonstandard monomials :");
    put(nonstd);
  end Enumerate_Brackets;

  procedure Read_Bracket ( b : in out Bracket ) is

    d : constant integer32 := b'last;
    sig : integer32;

  begin
    put("Give "); put(d,1); put(" row indices : "); get(b,sig);
    put("The sorted bracket : ");
    if sig > 0
     then put("+");
     else put("-");
    end if;
    put(b);
    if Is_Zero(b)
     then put_line(" is zero.");
     else put_line(" is different from zero.");
    end if;
  end Read_Bracket;

  procedure Test_Straighten ( d : in integer32 ) is

    b,bb : Bracket(1..d);
    s : natural32;
    ans : character;
    bp : Bracket_Polynomial;

  begin
    loop
      Read_Bracket(b);
      Read_Bracket(bb);
      if Is_Equal(b,bb)
       then put_line("Both brackets are equal.");
      end if;
      if b < bb then
        put(b); put(" < "); put(bb);
        s := Brackets.Is_Standard(b,bb);
        put(" and "); put(b); put(bb);
        bp := Straightening_Syzygy(b,bb);
      else
        put(b); put(" >= "); put(bb); 
        s := Brackets.Is_Standard(bb,b);
        put(" and "); put(bb); put(b);
        bp := Straightening_Syzygy(bb,b);
      end if;
      if s = 0
       then put_line(" is a Standard monomial.");
       else put(" is not a Standard monomial, with s = "); put(s,1); new_line;
      end if;
      put_line("The straightening syzygy : "); put(bp);
      put("Do you want to test more pairs of brackets ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Straighten;

  procedure Test_Laplace ( n,d : integer32 ) is

    p : constant Bracket_Polynomial
      := Laplace_Expansion(natural32(n),natural32(d));

  begin
    put("The Laplace expansion of "); put(n,1); put("*"); put(n,1);
    put("-determinant as product of "); put(d,1); put("- and ");
    put(n-d,1); put_line("-blocks : ");
    put(p);
  end Test_Laplace;

  function Decompose ( n,d : integer32; p : Bracket_Polynomial )
                     return Bracket_Polynomial is

  -- DESCRIPTION :
  --   Decomposes the ideal of square-free nonStandard monomials.

    res : Bracket_Polynomial;
    np : constant integer32 := integer32(Number_of_Monomials(p));
    brackmat : array(1..np,1..integer32(2)) of Link_to_Bracket;

    procedure Initialize is

    -- DESCRIPTION :
    --   Initializes the matrix of the brackets that occur in p.

      cnt_term : integer32 := 0;

      procedure List_Term ( t : in Bracket_Term; continue : out boolean ) is

        cnt_mon : integer32 := 0;

        procedure Store_Bracket ( b : in Bracket; continue : out boolean ) is
        begin
          cnt_mon := cnt_mon + 1;
          brackmat(cnt_term,cnt_mon) := new Bracket'(b);
          continue := true;
        end Store_Bracket;
        procedure Store_Brackets is
          new Bracket_Monomials.Enumerate_Brackets(Store_Bracket);

      begin
        cnt_term := cnt_term+1;
        Store_Brackets(t.monom);
        continue := true;
      end List_Term;
      procedure List_Terms is new Enumerate_Terms(List_Term);

    begin
      List_Terms(p);
    end Initialize;

   -- procedure Write is
   --
    -- DESCRIPTION :
    --   Writes the matrix of brackets.
   --
   -- begin
   --   for i in 1..np loop
   --     for j in 1..2 loop
   --       put("  ");
   --       put("b("); put(i,1); put(","); put(j,1); put(") :");
   --       put(brackmat(i,j).all);
   --     end loop;
   --     new_line;
   --   end loop;
   -- end Write;

   -- function Occured_Yet ( k : natural; b : Bracket ) return boolean is
   --
    -- DESCRIPTION :
    --   Returns true if the bracket b occurs in one of the rows <= k.
   --
   -- begin
   --   for i in 1..k loop
   --     for j in 1..2 loop
   --       if Is_Equal(brackmat(i,j).all,b)
   --        then return true;
   --       end if;
   --     end loop;
   --   end loop;
   --   return false;
   -- end Occured_Yet;

    procedure Solve is

    -- DESCRIPTION :
    --   Solves the initial ideal of quadratic nonStandard monomials.

      accu : Bracket_Monomial;
      nzrb : constant integer32
           := integer32(Number_of_Zero_Brackets(natural32(n),natural32(d)));

      procedure Solve ( k : in integer32; acc : in Bracket_Monomial ) is
      begin
        if k > np then
          Add(res,Create(acc));
        else
          if Divisible(acc,brackmat(k,1).all)
              or else Divisible(acc,brackmat(k,2).all) then
            Solve(k+1,acc);
          else 
            if Number_of_Brackets(acc) < natural32(nzrb) then
              for i in 1..integer32(2) loop
                declare
                  newacc : Bracket_Monomial;
                begin
                  Copy_Multiply(acc,newacc);
                  Multiply(newacc,brackmat(k,i).all);
                  Solve(k+1,newacc);
                end;
              end loop;
            end if;
          end if;
        end if;
      end Solve;

    begin
      Solve(1,accu);
    end Solve;

  begin
    Initialize;
   -- Write;
    Solve;
    return res;
  end Decompose;

  procedure Enumerate_Straightening_Syzygies ( n,d : in integer32) is

  -- DESCRIPTION :
  --   Sets up of the straightening syzygies for all nonStandard
  --   quadratic monomials.

    nonstd,decom : Bracket_Polynomial;

    procedure Write_Syzygy ( s : in Bracket_Polynomial; cont : out boolean ) is
    begin
      put(s); cont := true;
    end Write_Syzygy;
    procedure Write_Syzygies is new Enumerate_Syzygies(Write_Syzygy);

  begin
    put("(n,d) = "); put("("); put(n,1); put(","); put(d,1); put(")");
    put("  #Brackets : ");
    put(Number_of_Brackets(natural32(n),natural32(d)),1); 
    put("  #Linear Equations : ");
    put(Number_of_Linear_Equations(natural32(n),natural32(d)),1); new_line;
    put_line("Type of linear equation :");
    put(Laplace_Expansion(natural32(n),natural32(n-d)));
    nonstd := nonStandard_Monomials(natural32(n),natural32(d));
    put_line("The polynomial with all quadratic nonStandard monomials :");
    put(nonstd);
    put_line("All quadratic straightening syzygies : ");
    Write_Syzygies(nonstd);
    put("#Zero-Brackets : ");
    put(Number_of_Zero_Brackets(natural32(n),natural32(d)),1); new_line;
    decom := Decompose(n,d,nonstd);
    put_line("Decomposition of the ideal of nonStandard monomials :");
    put(decom);
    put("Number of solutions : "); put(Number_of_Monomials(decom),1);
  end Enumerate_Straightening_Syzygies;

  procedure Main is

    n,d : integer32 := 0; 
    ans : character;

  begin
    new_line;
    put_line("Interactive testing of straightening algorithm.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. Exit this program.");
      put_line("  1. Enumerate all brackets.");
      put_line("  2. Test the straightening algorithm.");
      put_line("  3. Test Laplace expansion.");
      put_line("  4. Enumerate straightening syzygies.");
      put("Make your choice (0,1,2,3,or 4) : "); get(ans);
      exit when ans = '0';
      new_line;
      put("Give the number of entries in bracket : "); get(d);
      put("Give the number of elements to choose from : "); get(n);
      case ans is
        when '1' => Enumerate_Brackets(n,d);
        when '2' => Test_Straighten(d);
        when '3' => Test_Laplace(n,d);
        when '4' => Enumerate_Straightening_Syzygies(n,d);
        when others => put_line("Bad answer.");
      end case;
    end loop;
  end Main;
begin
  Main;
end ts_straighten;
