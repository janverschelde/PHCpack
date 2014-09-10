--with text_io; use text_io;
--with Multprec_Integer_Numbers_io; use Multprec_Integer_Numbers_io;

package body Multprec_Common_Divisors is

  function pos_gcd ( a,b : Integer_Number ) return Integer_Number is

    res,r : Integer_Number;

  begin
    if Equal(b,0) then
      Copy(a,res);
    else
      r := Rmd(a,b);
      if Equal(r,0)
       then Copy(b,res);
       else res := pos_gcd(b,r);
      end if;
      Clear(r);
    end if;
    return res;
  end pos_gcd;

  function gcd ( a,b : Integer_Number ) return Integer_Number is

    res,mina,minb : Integer_Number;

  begin
    if a < 0 then
      if b < 0 then
        mina := -a; minb := -b;
        res := pos_gcd(mina,minb);
        Clear(mina); Clear(minb);
      else
        mina := -a;
        res := pos_gcd(mina,b);
        Clear(mina);
      end if;
    elsif b < 0 then
      minb := -b;
      res := pos_gcd(a,minb);
      Clear(minb);
    else
      res := pos_gcd(a,b);
    end if;
    return res;
  end gcd;

  function lcm ( a,b : Integer_Number ) return Integer_Number is

    res : Integer_Number;
    g : Integer_Number := gcd(a,b);

  begin
    if Equal(g,0)
     then res := Create(integer(0));
     else res := a/g; Mul(res,b);
    end if;
    Clear(g);
    return res;
  end lcm;

  procedure pos_gcd ( a,b : in Integer_Number; k,l,d : out Integer_Number ) is

    r,aa,bb,k1,l1,k0,l0,h,tmp : Integer_Number;

  -- REQUIRED :
  --   a >= b and a > 0 and b > 0

  begin
   -- put("Calling pos_gcd with a = "); put(a);
   -- put(" and b = "); put(b); new_line;
    r := Rmd(a,b);
    if Equal(r,0) then
      k := Create(integer(0)); l := Create(integer(1)); 
      Copy(b,d);
    else
      k0 := Create(integer(0)); l0 := Create(integer(1));
     -- bb := Create(0); 
      Copy(r,bb);
     -- aa := Create(0); 
      Copy(b,aa);
      k1 := Create(integer(1)); l1 := -a; Div(l1,b); --l1 := (-a)/b;
      Clear(r);
      loop
        r := Rmd(aa,bb);
        exit when Equal(r,0);
        h := aa / bb;
       -- nk1 := k0 - h*k1; Clear(k0); k0 := k1; k1 := nk1;
       -- nl1 := l0 - h*l1; Clear(l0); l0 := l1; l1 := nl1;
        Copy(k1,tmp); k1 := k0 - h*k1; Copy(tmp,k0);
        Copy(l1,tmp); l1 := l0 - h*l1; Copy(tmp,l0);
        Clear(h); Clear(tmp);
        Copy(bb,aa);
        Copy(r,bb);
        Clear(r);
      end loop;
      k := k1; l := l1; d := bb;
    end if;
 -- exception
 --   when others => put_line("exception in pos_gcd");
 --     put("a = "); put(a); put("  b = "); put(b); new_line; raise;
  end pos_gcd;

  procedure gcd ( a,b : in Integer_Number; k,l,d : out Integer_Number ) is

    mina,minb,kk,ll : Integer_Number;

  begin
   -- put("calling gcd with a = "); put(a); put(" and b = "); put(b); new_line;
    if Equal(a,0) then
      if b < 0
       then d := -b; k := create(integer(0)); l := create(integer(1)); Min(l); 
       else Copy(b,d); k := create(integer(0)); l := create(integer(1));
      end if;
      return;
    elsif Equal(b,0) then
      if a < 0
       then d := -a; k := create(integer(1)); Min(k); l := create(integer(0)); 
       else Copy(a,d); k := create(integer(1)); l := create(integer(0));
      end if;
      return;
    end if;
    if a < b then
      if b < 0 then
        minb := -b; mina := -a;
        pos_gcd(mina,minb,kk,ll,d);
        Clear(minb); Clear(mina);
        Min(kk); k := kk;
        Min(ll); l := ll;
      elsif a < 0 then
        mina := -a;
        if b > mina
         then pos_gcd(b,mina,l,kk,d);
         else pos_gcd(mina,b,kk,l,d);
        end if;
        Clear(mina);
        Min(kk); k := kk;
      else
        pos_gcd(a=>b,b=>a,k=>l,l=>k,d=>d);
      end if;
    else -- a >= b
      if a < 0 then
        mina := -a; minb := -b;
        pos_gcd(minb,mina,ll,kk,d);
        Clear(mina); Clear(minb);
        Min(kk); k := kk;
        Min(ll); l := ll;
      elsif b < 0 then
        minb := -b;
        if a > minb
         then pos_gcd(a,minb,k,ll,d);
         else pos_gcd(minb,a,ll,k,d);
        end if;
        Clear(minb);
        Min(ll); l := ll;
      else
        pos_gcd(a,b,k,l,d);
      end if;
    end if;
  end gcd;

end Multprec_Common_Divisors;
