with unchecked_deallocation;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;

package body Standard_Laurent_Homotopy is

-- INTERNAL DATA -- to store evaluable forms of the polynomials

  type homtype is (nat,art);  -- natural or artificial parameter homotopy

  -- In the homdata, the integer parameters are as follows:
  -- nq : the number of polynomials, nv : the number of variables,
  -- and nv1 = nv + 1.

  type homdata ( ht : homtype; nq,nv,nv1 : integer32 ) is record

    p : Laur_Sys(1..nq);                      -- target system
    pe : Eval_Laur_Sys(1..nq);                -- evaluable form of target
    dh : Jaco_Mat(1..nq,1..nv1);              -- Jacobian matrix of homotopy
    dhe : Eval_Jaco_Mat(1..nq,1..nv1);        -- evaluable form of dh

    case ht is
      when nat =>
        i : integer32;                        -- which variable is parameter
      when art =>
        q,h : Laur_Sys(1..nq);                -- start system and homotopy
        qe,he : Eval_Laur_Sys(1..nq);         -- evaluable form of q and h
        dpe,dqe : Eval_Jaco_Mat(1..nq,1..nv); -- evaluable Jacobians
        k : natural32;                        -- relaxation power
        gamma,beta : Vector(1..nq);           -- accessibility constants
        linear : boolean;                     -- linear in t if true
      end case;

  end record;
 
  type homtp is access homdata;
  hom : homtp;

-- GENERAL AUXILIARIES :

  function Mul ( v : Vector; p : Laur_Sys ) return Laur_Sys is

  -- DESCRIPTION :
  --   Multiplies the ith polynomial with the ith entry in v.

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := v(i)*p(i);
    end loop;
    return res;
  end Mul;

  function Mul ( m : Matrix; v : Vector ) return Matrix is

  -- DESCRIPTION :
  --   Multiplies the ith row of m with the ith entry in v.

    res : Matrix(m'range(1),m'range(2));

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        res(i,j) := v(i)*m(i,j);
      end loop;
    end loop;
    return res;
  end Mul;

-- AUXILIARIES TO THE CONSTRUCTORS :

  function Linear_Start_Factor
             ( n : in integer32; k : in natural32;
               a : in Complex_Number ) return Poly is

  -- DESCRIPTION :
  --   Returns a*(1-t)^k, with t as the (n+1)-th variable in the polynomial.

    res,tmp : Poly;
    t : Term;

  begin
    t.cf := a;
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+1 => 0);
    res := Create(t);                                           -- res = a
    t.cf := Create(1.0);
    tmp := Create(t);                                           -- tmp = 1
    t.dg(n+1) := 1;
    Sub(tmp,t);                                                 -- tmp = 1-t
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(t.dg));
    for i in 1..k loop
      Mul(res,tmp);
    end loop;
    Clear(tmp);
    return res;
  end Linear_Start_Factor;

  function Nonlinear_Start_Factor
             ( n : in integer32; k : in natural32;
               a : in Complex_Number ) return Poly is

  -- DESCRIPTION :
  --   Returns (1 - (t - t*(1-t)*a))^k, with t as the (n+1)-th variable
  --   in the polynomial.

    res,tmp : Poly;
    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+1 => 0);
    res := Create(t);                                         -- res = 1
    tmp := Create(t);                                         -- tmp = 1
    t.dg(n+1) := 1;
    Sub(tmp,t);                                           -- tmp = 1-t
    Mul(tmp,t);                                           -- tmp = t*(1-t)
    Mul(tmp,-a);                                          -- tmp = -a*t*(1-t)
    Add(tmp,t);                                       -- tmp = t - a*t*(1-t)
    Sub(res,tmp);                                 -- res = 1 - t + a*t*(1-t)
    if k > 1 then
      Copy(res,tmp);
      for i in 1..(k-1) loop
        Mul(res,tmp);
      end loop;
    end if;
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(t.dg));
    Clear(tmp);
    return res;
  end Nonlinear_Start_Factor;

  function Linear_Target_Factor
             ( n : in integer32; k : in natural32 ) return Poly is

  -- DESCRIPTION :
  --   Returns t^k, with t as the (n+1)-th variable in the polynomial.

    res : Poly;
    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+1 => 0);
    t.dg(n+1) := integer32(k);
    res := Create(t);
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(t.dg));
    return res;
  end Linear_Target_Factor;

  function Linear_Target_Factor
             ( n : in integer32; k : in natural32;
               a : in Complex_Number ) return Poly is

  -- DESCRIPTION :
  --   Returns a*t^k, with t as the (n+1)-th variable in the polynomial.

    res : Poly;
    t : Term;

  begin
    t.cf := a;
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+1 => 0);
    t.dg(n+1) := integer32(k);
    res := Create(t);
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(t.dg));
    return res;
  end Linear_Target_Factor;

  function Nonlinear_Target_Factor
             ( n : in integer32; k : in natural32;
               a : in Complex_Number ) return Poly is

  -- DESCRIPTION :
  --   Returns (t - t*(1-t)*a)^k, with t as the (n+1)-th variable in
  --   the polynomial.

    res,tmp : Poly;
    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+1 => 0);
    tmp := Create(t);                                         -- tmp = 1
    t.dg(n+1) := 1;
    res := Create(t);                                         -- res = t
    Sub(tmp,t);                                           -- tmp = 1-t
    Mul(tmp,t);                                           -- tmp = t*(1-t)
    Mul(tmp,a);                                           -- tmp = a*t*(1-t)
    Sub(res,tmp);                                     -- res = t - a*t*(1-t)
    if k > 1 then
      Copy(res,tmp);
      for i in 1..(k-1) loop
        Mul(res,tmp);
      end loop;
    end if;
    Standard_Integer_Vectors.Clear
      (Standard_Integer_Vectors.Link_to_Vector(t.dg));
    Clear(tmp);
    return res;
  end Nonlinear_Target_Factor;

  function Plus_one_Unknown ( p : in Poly ) return Poly is

  -- DESCRIPTION :
  --   The returning polynomial has place for one additional unknown.

    res : Poly;

    procedure Plus_Unknown_In_Term ( t : in out Term; c : out boolean ) is

      n : constant integer32 := t.dg'length;
      temp : Standard_Integer_Vectors.Vector(1..(n+1));

    begin
      temp(1..n) := t.dg.all;
      temp(n+1) := 0;
      Standard_Integer_Vectors.Clear
        (Standard_Integer_Vectors.Link_to_Vector(t.dg));
      t.dg := new Standard_Integer_Vectors.Vector'(temp);
      c := true;
    end Plus_Unknown_In_Term;
    procedure Plus_Unknown_In_Terms is 
      new Changing_Iterator (process => Plus_Unknown_In_Term);

  begin
    Copy(p,res);
    Plus_Unknown_in_Terms(res);
    return res;
  end Plus_one_Unknown;

  function Linear_Homotopy
             ( p,q : in Laur_Sys; k : in natural32; a : in Complex_Number )
             return Laur_Sys is

    h : Laur_Sys(p'range);
    tempp,tempq,q_fac,p_fac : Poly;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    q_fac := Linear_Start_Factor(nv,k,a);
    p_fac := Linear_Target_Factor(nv,k);
    for i in h'range loop
      tempq := Plus_one_Unknown(q(i));
      tempp := Plus_one_Unknown(p(i));
      Mul(tempq,q_fac);
      Mul(tempp,p_fac);
      h(i) := tempq + tempp;
      Clear(tempq); Clear(tempp);
    end loop;
    Clear(q_fac); Clear(p_fac);
    return h;
  end Linear_Homotopy;

  function Linear_Homotopy
             ( p,q : in Laur_Sys; k : in natural32; a : in Vector )
             return Laur_Sys is

    h : Laur_Sys(p'range);
    tempp,tempq,q_fac,p_fac : Poly;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    for i in h'range loop
      q_fac := Linear_Start_Factor(nv,k,a(i));
      p_fac := Linear_Target_Factor(nv,k);
      tempq := Plus_one_Unknown(q(i));
      tempp := Plus_one_Unknown(p(i));
      Mul(tempq,q_fac);
      Mul(tempp,p_fac);
      h(i) := tempq + tempp;
      Clear(tempq); Clear(tempp);
      Clear(q_fac); Clear(p_fac);
    end loop;
    return h;
  end Linear_Homotopy;

  function Linear_Homotopy
             ( p,q : in Laur_Sys; k : in natural32; a,b : in Vector )
             return Laur_Sys is

    h : Laur_Sys(p'range);
    tempp,tempq,q_fac,p_fac : Poly;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    for i in h'range loop
      q_fac := Linear_Start_Factor(nv,k,a(i));
      p_fac := Linear_Target_Factor(nv,k,b(i));
      tempq := Plus_one_Unknown(q(i));
      tempp := Plus_one_Unknown(p(i));
      Mul(tempq,q_fac);
      Mul(tempp,p_fac);
      h(i) := tempq + tempp;
      Clear(tempq); Clear(tempp);
      Clear(q_fac); Clear(p_fac);
    end loop;
    return h;
  end Linear_Homotopy;

  function Nonlinear_Homotopy
             ( p,q : in Laur_Sys; k : in natural32; a,b : in Vector )
             return Laur_Sys is

    h : Laur_Sys(p'range);
    tempp,tempq,q_fac,p_fac : Poly;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    for i in h'range loop
      q_fac := Nonlinear_Start_Factor(nv,k,a(i));
      p_fac := Nonlinear_Target_Factor(nv,k,b(i));
      tempq := Plus_one_Unknown(q(i));
      tempp := Plus_one_Unknown(p(i));
      Mul(tempq,q_fac);
      Mul(tempp,p_fac);
      h(i) := tempq + tempp;
      Clear(tempq); Clear(tempp);
      Clear(q_fac); Clear(p_fac);
    end loop;
    return h;
  end Nonlinear_Homotopy;

  procedure Create ( p,q : in Laur_Sys; k : in natural32;
                     a : in Complex_Number ) is

    nq : constant integer32 := p'length;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    dp,dq : Jaco_Mat(1..nq,1..nv);
    ho : homdata(art,nq,nv,nv+1);

  begin
    Copy(p,ho.p); Copy(q,ho.q);
    ho.h := Linear_Homotopy(p,q,k,a);
    ho.pe := Create(ho.p);
    ho.qe := Create(ho.q);
    ho.he := Create(ho.h);
    dp := Create(ho.p);
    dq := Create(ho.q);
    ho.dh := Create(ho.h);
    ho.dpe := Create(dp);
    ho.dqe := Create(dq);
    ho.dhe := Create(ho.dh);
    Clear(dp); Clear(dq);
    ho.k := k;
    for i in 1..nq loop
      ho.gamma(i) := a;
    end loop;
    for i in 1..nq loop
      ho.beta(i) := Create(1.0);
    end loop;
    ho.linear := true;
    hom := new homdata'(ho);
  end Create; 

  procedure Create ( p,q : in Laur_Sys; k : in natural32; a : in Vector ) is

    nq : constant integer32 := p'length;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    dp,dq : Jaco_Mat(1..nq,1..nv);
    ho : homdata(art,nq,nv,nv+1);

  begin
    Copy(p,ho.p); Copy(q,ho.q);
    ho.h := Linear_Homotopy(p,q,k,a);
    ho.pe := Create(ho.p);
    ho.qe := Create(ho.q);
    ho.he := Create(ho.h);
    dp := Create(ho.p);
    dq := Create(ho.q);
    ho.dh := Create(ho.h);
    ho.dpe := Create(dp);
    ho.dqe := Create(dq);
    ho.dhe := Create(ho.dh);
    Clear(dp); Clear(dq);
    ho.k := k;
    ho.gamma := a;
    for i in 1..nq loop
      ho.beta(i) := Create(1.0);
    end loop;
    ho.linear := true;
    hom := new homdata'(ho);
  end Create;

  procedure Create ( p,q : in Laur_Sys; k : in natural32; a,b : in Vector;
                     linear : in boolean ) is

    nq : constant integer32 := p'length;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    dp,dq : Jaco_Mat(1..nq,1..nv);
    ho : homdata(art,nq,nv,nv+1);

  begin
    Copy(p,ho.p); Copy(q,ho.q);
    ho.linear := linear;
    if linear
     then ho.h := Linear_Homotopy(p,q,k,a,b);
     else ho.h := Nonlinear_Homotopy(p,q,k,a,b);      
    end if;
    ho.pe := Create(ho.p);
    ho.qe := Create(ho.q);
    ho.he := Create(ho.h);
    dp := Create(ho.p);
    dq := Create(ho.q);
    ho.dh := Create(ho.h);
    ho.dpe := Create(dp);
    ho.dqe := Create(dq);
    ho.dhe := Create(ho.dh);
    Clear(dp); Clear(dq);
    ho.k := k;
    ho.gamma := a;
    ho.beta := b;
    hom := new homdata'(ho);
  end Create;

  procedure Create ( p : in Laur_Sys; k : in integer32 ) is

    nq : constant integer32 := p'last-p'first+1; 
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    ho : homdata(nat,nq,nv,nv+1);

  begin
    Copy(p,ho.p);
    ho.pe := Create(ho.p);
    ho.dh := Create(ho.p);
    ho.dhe := Create(ho.dh);
    ho.i := k;
    hom := new homdata'(ho);
  end Create;

-- SELECTOR :

  function Homotopy_System return Laur_Sys is

    ho : homdata renames hom.all;

  begin
    case ho.ht is
      when nat => return ho.p;
      when art => return ho.h;
    end case;
  end Homotopy_System;

-- SYMBOLIC ROUTINES :

  function Eval ( t : Complex_Number ) return Laur_Sys is

    ho : homdata renames hom.all;

  begin
    case ho.ht is
      when nat =>  -- t = x(ho.i)
        return Eval(ho.p,t,ho.i);
      when art =>  -- compute : a * ((1 - t)^k) * q + (t^k) * p
        declare
          p_factor,q_factor : Vector(1..ho.nq);
          res,tmp : Laur_Sys(1..ho.nq);
        begin
          if ho.linear then
            if AbsVal(t) = 0.0 then
              return Mul(ho.gamma,ho.q);
            elsif abs(REAL_PART(t) - 1.0 ) + 1.0 = 1.0 
                and then abs(IMAG_PART(t)) + 1.0 = 1.0 then
              return Mul(ho.beta,ho.p);
            else
              for i in 1..ho.nq loop
                q_factor(i) := ho.gamma(i);
                p_factor(i) := ho.beta(i);
                for i in 1..integer32(ho.k) loop
                  q_factor(i) := (Create(1.0)-t) * q_factor(i);
                  p_factor(i) := t * p_factor(i);
                end loop;
              end loop;
              res := Mul(p_factor,ho.p); 
              tmp := Mul(q_factor,ho.q);
              Add(res,tmp);
              Clear(tmp);
              return res;
            end if;
          else
            return res;  -- still to do !
          end if;
        end;
    end case;
  end Eval;

  function Diff ( t : Complex_Number ) return Laur_Sys is

    ho : homdata renames hom.all;

  begin
    case ho.ht is
      when nat =>  -- t = x(ho.i)
        return Diff(ho.p,ho.i);
      when art =>  -- compute  - a*k*(1 - t)^(k-1)*q + b*k*t^(k-1)*p
        declare
          q_factor,p_factor : Vector(1..ho.nq);
          tmp : Laur_Sys(1..ho.nq);
          res : Laur_Sys(1..ho.nq);
        begin
          if ho.linear then
            for i in 1..ho.nq loop
              q_factor(i) := (-ho.gamma(i)) * Create(double_float(ho.k));
              p_factor(i) := ho.beta(i) * Create(double_float(ho.k));
            end loop;
            if AbsVal(t) = 0.0 then
              if ho.k = 1 then
                res := Mul(p_factor,ho.p);
                tmp := Mul(q_factor,ho.q);
                Add(res,tmp);
                Clear(tmp);
                return res;
              else
                return Mul(q_factor,ho.q);
              end if;
            elsif abs( REAL_PART(t) - 1.0 ) + 1.0 = 1.0
                 and then abs(IMAG_PART(t)) + 1.0 = 1.0 then
              return Create(double_float(ho.k)) * ho.p;
            else
              for i in 1..(ho.k-1) loop
                q_factor := (Create(1.0)-t) * q_factor;
                p_factor := t * p_factor;
              end loop;
              res := Mul(p_factor,ho.p);
              tmp := Mul(q_factor,ho.q);
              Add(res,tmp);
              Clear(tmp);
              return res;
            end if;
          else
            return res;   -- still left to do !!!
          end if;
        end;
    end case;
  end Diff;

-- NUMERIC ROUTINES :

  function Eval ( x : Vector; t : Complex_Number ) return Vector is

    ho : homdata renames hom.all;

  begin
    case ho.ht is
      when nat =>
        declare
          y : Vector(x'first..x'last+1);
        begin
          y(1..ho.i-1) := x(1..ho.i-1);
          y(ho.i) := t;
          y(ho.i+1..y'last) := x(ho.i..x'last);
          return Eval(ho.pe,y);
        end;
      when art =>
        if AbsVal(t) + 1.0 = 1.0 then
          if ho.linear
           then return ho.gamma * Eval(ho.qe,x);
           else return Eval(ho.qe,x);
          end if;
        elsif abs( REAL_PART(t) - 1.0 ) + 1.0 = 1.0
             and then abs(IMAG_PART(t)) + 1.0 = 1.0 then
          if ho.linear
           then return ho.beta * Eval(ho.pe,x);
           else return Eval(ho.pe,x);
          end if;
        else
          declare
            y : Vector(x'first..x'last+1);
          begin
            y(x'range) := x;
            y(y'last) := t;
            return Eval(ho.he,y);
          end;
        end if;
    end case;
  end Eval;

  function Diff ( x : Vector; t : Complex_Number ) return Vector is

    n : constant integer32 := x'length;

  begin
    case hom.ht is
      when nat => return Diff(x,t,hom.i);
      when art => return Diff(x,t,n+1);
    end case;
  end Diff;
 
  function Diff ( x : Vector; t : Complex_Number ) return Matrix is

    ho : homdata renames hom.all;
    nq : integer32 renames ho.nq;
    nv : integer32 renames ho.nv;

  begin
    case ho.ht is
      when nat =>
        declare
          m : Matrix(1..nq,1..nv);
          y : Vector(1..nv+1);
        begin
          y(1..ho.i-1) := x(1..ho.i-1);
          y(ho.i) := t;
          y(ho.i+1..nv+1) := x(ho.i..nv);
          for i in 1..nq loop
            for j in 1..nv loop
              m(i,j) := Eval(ho.dhe(i,j),y);
            end loop;
          end loop;
          return m;
        end;
      when art =>
        if AbsVal(t) + 1.0 = 1.0 then
          declare
            m : constant Matrix(1..nq,1..nv) := Eval(ho.dqe,x);
          begin
            if ho.linear
             then return Mul(m,ho.gamma);
             else return m;
            end if;
          end;
        elsif abs( REAL_PART(t) - 1.0 ) + 1.0 = 1.0
             and then abs(IMAG_PART(t)) + 1.0 = 1.0 then
          declare
            m : constant Matrix(1..nq,1..nv) := Eval(ho.dpe,x);
          begin
            if ho.linear
             then return Mul(m,ho.beta);
             else return m;
            end if;
          end;
        else
          declare
            m : Matrix(1..nq,1..nv);
            y : Vector(1..nv+1);
          begin
            y(1..nv) := x;
            y(nv+1) := t;
            for i in 1..nq loop
              for j in 1..nv loop
                m(i,j) := Eval(ho.dhe(i,j),y);
              end loop;
            end loop;
            return m;
          end;
        end if;
    end case;
  end Diff;

  function Diff ( x : Vector; t : Complex_Number;
                  k : integer32 ) return Vector is

    ho : homdata renames hom.all;
    nq : integer32 renames ho.nq;
    nv : integer32 renames ho.nv;
    y : Vector(1..nv+1);
    res : Vector(1..nq);

  begin
    case ho.ht is
      when nat => y(1..ho.i-1) := x(1..ho.i-1);
                  y(ho.i) := t;
                  y(ho.i+1..nv+1) := x(ho.i..nv);
      when art => y(1..nv) := x;
                  y(nv+1) := t;
    end case;
    for i in 1..nq loop
      res(i) := Eval(ho.dhe(i,k),y);
    end loop;
    return res;
  end Diff;

  function Diff ( x : Vector; t : Complex_Number;
                  k : integer32 ) return Matrix is

    ho : homdata renames hom.all;
    nq : integer32 renames ho.nq;
    nv : integer32 renames ho.nv;
    y : Vector(1..nv+1);
    res : Matrix(1..nq,1..nv);

  begin
    case ho.ht is
      when nat => y(1..ho.i-1) := x(1..ho.i-1);
                  y(ho.i) := t;
                  y(ho.i+1..nv+1) := x(ho.i..nv); 
      when art => y(1..nv) := x;
                  y(nv+1) := t;
    end case;
    for j in 1..(k-1) loop
      for i in 1..nq loop
        res(i,j) := Eval(ho.dhe(i,j),y);
      end loop;
    end loop;
    for j in (k+1)..(nv+1) loop
      for i in 1..nq loop
        res(i,j-1) := Eval(ho.dhe(i,j),y);
      end loop;
    end loop;
    return res;
  end Diff;

-- DESTRUCTOR :

  procedure free is new unchecked_deallocation (homdata,homtp);

  procedure Clear is
  begin
    if hom /= null then
      declare
        ho : homdata renames hom.all;
      begin
        Clear(ho.p);  Clear(ho.pe);
        Clear(ho.dh); Clear(ho.dhe);
        case ho.ht is
          when nat => null;
          when art =>
            Clear(ho.q);   Clear(ho.qe);
            Clear(ho.h);   Clear(ho.he);
            Clear(ho.dpe); Clear(ho.dqe);
        end case;
      end;
      free(hom);
    end if;
  end Clear;

end Standard_Laurent_Homotopy;
