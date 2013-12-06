with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Brackets,Brackets_io;             use Brackets,Brackets_io;

procedure ts_brackets is

-- DESCRIPTION :
--   Enumerates all brackets in the Pluecker embedding.

  procedure Enumerate ( n,d,k,start : in integer32; accu : in out Bracket;
                        cnt : in out integer32 ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of all brackets.

  begin
    if k > d then
      -- put(accu); new_line;
      cnt := cnt + 1;
    else
      for ell in start..n loop
        accu(k) := natural32(ell);
        Enumerate(n,d,k+1,ell+1,accu,cnt);
      end loop;
    end if;
  end Enumerate;

  procedure Enumerate2 ( n,d,k,start : in integer32; b,accu : in out Bracket;
                         cntstd,cntnonstd : in out integer32 ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of all brackets, with b < accu and with
  --   a test whether the pair b*accu forms a Standard monomial.

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
        end if;
      end if;
    else
      for ell in start..n loop
        accu(k) := natural32(ell);
        Enumerate2(n,d,k+1,ell+1,b,accu,cntstd,cntnonstd);
      end loop;
    end if;
  end Enumerate2;

  procedure Enumerate1 ( n,d,k,start : integer32; acc1,acc2 : in out Bracket;
                         cntstd,cntnonstd : in out integer32 ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of all brackets, with acc1 < acc2 and with
  --   a test whether the pair acc1*acc2 forms a Standard monomial.
  --   Counts the standard and nonstandard monomials.

  begin
    if k > d then
      Enumerate2(n,d,1,integer32(acc1(1)),acc1,acc2,cntstd,cntnonstd);
    else
      for ell in start..n loop
        acc1(k) := natural32(ell);
        Enumerate1(n,d,k+1,ell+1,acc1,acc2,cntstd,cntnonstd);
      end loop;
    end if;
  end Enumerate1;

  procedure Enumerate_Pairs ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Enumerates all pairs (b1,b2), with b1 < b2 and checks whether
  --   the corresponding bracket monomial b1*b2 is standard or not.

    b1,b2 : Bracket(1..d);
    cntstd,cntnonstd : integer32 := 0;

  begin
    Enumerate1(n,d,1,1,b1,b2,cntstd,cntnonstd);
    put("#Standard bi-monomials    : "); put(cntstd,1);    new_line;
    put("#nonStandard bi-monomials : "); put(cntnonstd,1); new_line;
  end Enumerate_Pairs;

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

  procedure Test_Sort ( d : in integer32 ) is

    b,bb : Bracket(1..d);
    s : natural32;
    ans : character;

  begin
    loop
      Read_Bracket(b);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
      Read_Bracket(bb);
      if Is_Equal(b,bb)
       then put_line("Both brackets are equal.");
      end if;
      if b < bb then
        put(b); put(" < "); put(bb);
        s := Brackets.Is_Standard(b,bb);
        put(" and "); put(b); put(bb);
      else
        put(b); put(" >= "); put(bb); 
        s := Brackets.Is_Standard(bb,b);
        put(" and "); put(bb); put(b);
      end if;
      if s = 0
       then put_line(" is a Standard monomial.");
       else put(" is not a Standard monomial, with s = "); put(s,1); new_line;
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Sort;

  procedure Main is

    n,d : integer32 := 0;
    ans : character;

  begin
    loop
      put("Give the number of entries in bracket : "); get(d);
      put("Give the number of elements to choose from : "); get(n);
      declare
        acc : Bracket(1..d);
        cnt : integer32 := 0;
      begin
        Enumerate(n,d,1,1,acc,cnt);
        put("#brackets : "); put(cnt,1); new_line;
      end;
      Enumerate_Pairs(n,d);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
    Test_Sort(d);
  end Main;

begin
  Main;
end ts_brackets;
