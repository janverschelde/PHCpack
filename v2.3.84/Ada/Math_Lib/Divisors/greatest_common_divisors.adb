-- with text_io; use text_io;

package body Greatest_Common_Divisors is

  function pos_gcd ( a,b : number ) return number is

    res,r : number;

  begin
    res := Create(0);
    if Equal(b,zero) then
      Copy(a,res);
    else
      r := Rmd(a,b);
      if Equal(r,zero)
       then Copy(b,res);
       else res := pos_gcd(b,r);
      end if;
      Clear(r);
    end if;
    return res;
  end pos_gcd;

  function gcd ( a,b : number ) return number is

    res,mina,minb : number;

  begin
    if a < zero then
      if b < zero then
        mina := -a; minb := -b;
        res := pos_gcd(mina,minb);
        Clear(mina); Clear(minb);
      else
        mina := -a;
        res := pos_gcd(mina,b);
        Clear(mina);
      end if;
    elsif b < zero then
      minb := -b;
      res := pos_gcd(a,minb);
      Clear(minb);
    else
      res := pos_gcd(a,b);
    end if;
    return res;
  end gcd;

  function lcm ( a,b : number ) return number is

    res : number;
    g : number := gcd(a,b);

  begin
    if Equal(g,zero)
     then res := zero;
     else res := a/g; Mul(res,b);
    end if;
    Clear(g);
    return res;
  end lcm;

  procedure pos_gcd ( a,b : in number; k,l,d : out number ) is

    r,aa,bb,k1,l1,k0,l0,h,tmp : number;
   -- nk1,nl1 : number;

  -- REQUIRED :
  --   a >= b and a > 0 and b > 0

  begin
    r := Rmd(a,b);
    if Equal(r,zero) then
      k := Create(0); l := Create(1); Copy(b,d);
    else
      k0 := Create(0); l0 := Create(1);
      Copy(r,bb); Copy(b,aa);
      k1 := Create(1); l1 := -a; Div(l1,b); --l1 := (-a)/b;
      Clear(r);
      loop
        r := Rmd(aa,bb);
        exit when Equal(r,zero);
        h := aa / bb;
       -- nk1 := k0 - h*k1; Clear(k0); k0 := k1; k1 := nk1;
       -- nl1 := l0 - h*l1; Clear(l0); l0 := l1; l1 := nl1;
        Copy(k1,tmp); k1 := k0 - h*k1; Copy(tmp,k0);
        Copy(l1,tmp); l1 := l0 - h*l1; Copy(tmp,l0);
        Clear(h); Clear(tmp);
        Copy(bb,aa); Copy(r,bb);
        Clear(r);
      end loop;
      k := k1; l := l1; d := bb;
    end if;
 -- exception
 --   when others => put_line("exception in pos_gcd"); raise;
  end pos_gcd;

  procedure gcd ( a,b : in number; k,l,d : out number ) is

    mina,minb,kk,ll : number;

  begin
    if Equal(a,zero) then
      if b < zero
       then d := -b; k := create(0); l := create(1); Min(l); 
       else Copy(b,d); k := create(0); l := create(1);
      end if;
      return;
    elsif Equal(b,zero) then
      if a < zero
       then d := -a; k := create(1); Min(k); l := create(0); 
       else Copy(a,d); k := create(1); l := create(0);
      end if;
      return;
    end if;
    if a < b then
      if b < zero then
        minb := -b; mina := -a;
        pos_gcd(mina,minb,kk,ll,d);
        Clear(minb); Clear(mina);
        Min(kk); k := kk;
        Min(ll); l := ll;
      elsif a < zero then
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
      if a < zero then
        mina := -a; minb := -b;
        pos_gcd(minb,mina,ll,kk,d);
        Clear(mina); Clear(minb);
        Min(kk); k := kk;
        Min(ll); l := ll;
      elsif b < zero then
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

end Greatest_Common_Divisors;
