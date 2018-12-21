with unchecked_deallocation;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;

package body Generic_Polynomial_Vectors is

-- CONSTRUCTORS :

  procedure Power_Update ( s : in out System ) is
  begin
    s.deg1 := Compute_Deg1(s.pols);
    Largest_Exponents(s.pols,s.maxexp);
    if not s.deg1
     then s.powtab := Allocate_Power_Table(s.maxexp);
    end if;
  end Power_Update;

  procedure Power_Update ( s : in out Link_to_System ) is
  begin
    if s /= null then
      s.deg1 := Compute_Deg1(s.pols);
      Largest_Exponents(s.pols,s.maxexp);
      if not s.deg1
       then s.powtab := Allocate_Power_Table(s.maxexp);
      end if;
    end if;
  end Power_Update;

  function Compute_Deg1 ( p : Polynomial_Vector ) return boolean is
  begin
    for i in p'range loop
      if not p(i).deg1
       then return false;
      end if;
    end loop;
    return true;
  end Compute_Deg1;

  function Compute_Deg1 ( p : Link_to_Polynomial_Vector ) return boolean is
  begin
    if p = null
     then return true;
     else return Compute_Deg1(p.all);
    end if;
  end Compute_Deg1;

  function Compute_Deg1 ( s : System ) return boolean is
  begin
    return Compute_Deg1(s.pols);
  end Compute_Deg1;

  function Compute_Deg1 ( s : Link_to_System ) return boolean is
  begin
    if s = null
     then return true;
     else return Compute_Deg1(s.pols);
    end if;
  end Compute_Deg1;

  procedure Largest_Exponents
              ( p : in Polynomial_Vector;
                e : out Standard_Natural_Vectors.Vector ) is

    w : Standard_Natural_Vectors.Vector(e'range);

  begin
    e := (e'range => 0);
    for i in p'range loop
      w := p(i).maxexp;
      for j in w'range loop
        if w(j) > e(j)
         then e(j) := w(j);
        end if;
      end loop;
    end loop;
  end Largest_Exponents;

  procedure Largest_Exponents
              ( p : in Link_to_Polynomial_Vector;
                e : out Standard_Natural_Vectors.Vector ) is
  begin
    if p /= null
     then Largest_Exponents(p.all,e);
    end if;
  end Largest_Exponents;

  procedure Largest_Exponents
              ( s : in System;
                e : out Standard_Natural_Vectors.Vector ) is
  begin
    Largest_Exponents(s.pols,e);
  end Largest_Exponents;

  procedure Largest_Exponents
              ( s : in Link_to_System;
                e : out Standard_Natural_Vectors.Vector ) is
  begin
    if s /= null
     then Largest_Exponents(s.pols,e);
    end if;
  end Largest_Exponents;

  function Allocate_Power_Table
              ( maxexp : Standard_Natural_Vectors.Vector )
              return VecVecs.Link_to_VecVec is
  begin
    return Polynomials.Allocate_Power_Table(maxexp);
  end Allocate_Power_Table;

  procedure Allocate_Power_Table
              ( powtab : in out VecVecs.VecVec;
                maxexp : in Standard_Natural_Vectors.Vector ) is
  begin
    Polynomials.Allocate_Power_Table(powtab,maxexp);
  end Allocate_Power_Table;

  procedure Update_Power_Table
              ( powtab : in out VecVecs.VecVec; x : in Vectors.Vector ) is
  begin
    Polynomials.Update_Power_Table(powtab,x);
  end Update_Power_Table;

  procedure Update_Power_Table
              ( powtab : in VecVecs.Link_to_VecVec; x : in Vectors.Vector ) is
  begin
    Polynomials.Update_Power_Table(powtab,x);
  end Update_Power_Table;

-- EVALUATION and DIFFERENTIATION :

  function Eval ( p : Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector is

    res : Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomials.Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Link_to_Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector is

    res : Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomials.Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( s : System; x : Vectors.Vector ) return Vectors.Vector is
  begin
    return Eval(s.pols,x);
  end Eval;

  function Eval ( s : Link_to_System;
                  x : Vectors.Vector ) return Vectors.Vector is
  begin
    return Eval(s.pols,x);
  end Eval;

  procedure Diff ( p : in Polynomial_Vector; x : in Vectors.Vector;
                   m : out Matrices.Matrix ) is

    yd,wrk : Vectors.Vector(x'range);

  begin
    for i in p'range loop
      Polynomials.Diff(p(i).mons,x,yd,wrk);
      for j in yd'range loop
        Ring.Copy(yd(j),m(i,j));
      end loop;
    end loop;
  end Diff;

  procedure Diff ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                   m : out Matrices.Matrix ) is
  begin
    if p /= null
     then Diff(p.all,x,m);
    end if;
  end Diff;

  procedure Diff ( s : in System; x : in Vectors.Vector;
                   m : out Matrices.Matrix ) is
  begin
    Diff(s.pols,x,m);
  end Diff;

  procedure Diff ( s : in Link_to_System; x : in Vectors.Vector;
                   m : out Matrices.Matrix ) is
  begin
    if s /= null
     then Diff(s.pols,x,m);
    end if;
  end Diff;

  procedure Speel ( s : in System; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    if s.deg1 then
      Speel_on_Product(s,x,y,m,yd,wrk);
    else
      Update_Power_Table(s.powtab,x);
      Speel(s.pols,x,s.powtab,y,m,yd,wrk);
    end if;
  end Speel;

  procedure Speel ( s : in System; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix ) is
  begin
    if s.deg1 then
      Speel_on_Product(s,x,y,m);
    else
      Update_Power_Table(s.powtab,x);
      Speel(s.pols,x,s.powtab,y,m);
    end if;
  end Speel;

  procedure Speel ( s : in Link_to_System; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    if s /= null then
      if s.deg1 then
        Speel_on_Product(s,x,y,m,yd,wrk);
      else
        Update_Power_Table(s.powtab,x);
        Speel(s.pols,x,s.powtab,y,m,yd,wrk);
      end if;
    end if;
  end Speel;

  procedure Speel ( s : in Link_to_System; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix ) is
  begin
    if s /= null then
      if s.deg1 then
        Speel_on_Product(s,x,y,m);
      else
        Update_Power_Table(s.powtab,x);
        Speel(s.pols,x,s.powtab,y,m);
      end if;
    end if;
  end Speel;

  procedure Speel_on_Product
              ( s : in System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector ) is
  begin
    for i in s.pols'range loop
      Polynomials.Speel_on_Product(s.pols(i),x,y(i),yd,wrk);
      for j in yd'range loop
        Ring.Copy(yd(j),m(i,j));
      end loop;
    end loop;
  end Speel_on_Product;

  procedure Speel_on_Product
              ( s : in System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix ) is

    yd,wrk : Vectors.Vector(x'range);

  begin
    Speel_on_Product(s,x,y,m,yd,wrk);
  end Speel_on_Product;

  procedure Speel_on_Product
              ( s : in Link_to_System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector ) is
  begin
    for i in s.pols'range loop
      Polynomials.Speel_on_Product(s.pols(i),x,y(i),yd,wrk);
      for j in yd'range loop
        Ring.Copy(yd(j),m(i,j));
      end loop;
    end loop;
  end Speel_on_Product;

  procedure Speel_on_Product
              ( s : in Link_to_System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix ) is

    yd,wrk : Vectors.Vector(x'range);

  begin
    Speel_on_Product(s,x,y,m,yd,wrk);
  end Speel_on_Product;

  procedure Speel ( p : in Polynomial_Vector; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    for i in p'range loop
      Polynomials.Speel(p(i).mons,x,powtab,y(i),yd,wrk);
      Ring.add(y(i),p(i).cff0);
      for j in yd'range loop
        Ring.Copy(yd(j),m(i,j));
      end loop;
    end loop;
  end Speel;

  procedure Speel ( p : in Polynomial_Vector; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : out Vectors.Vector; m : out Matrices.Matrix ) is

    yd,wrk : Vectors.Vector(x'range);

  begin
    Speel(p,x,powtab,y,m,yd,wrk);
  end Speel;

  procedure Speel ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    if p /= null
     then Speel(p.all,x,powtab,y,m,yd,wrk);
    end if;
  end Speel;

  procedure Speel ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : out Vectors.Vector; m : out Matrices.Matrix ) is
  begin
    if p /= null
     then Speel(p.all,x,powtab,y,m);
    end if;
  end Speel;

  procedure Speel_without_System_Cache
              ( s : in System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector ) is
  begin
    Speel_without_System_Cache(s.pols,x,y,m,yd,wrk);
  end Speel_without_System_Cache;

  procedure Speel_without_System_Cache
              ( s : in System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix ) is
  begin
    Speel_without_System_Cache(s.pols,x,y,m);
  end Speel_without_System_Cache;

  procedure Speel_without_System_Cache
              ( s : in Link_to_System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector ) is
  begin
    Speel_without_System_Cache(s.pols,x,y,m,yd,wrk);
  end Speel_without_System_Cache;

  procedure Speel_without_System_Cache
              ( s : in Link_to_System; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix ) is
  begin
    Speel_without_System_Cache(s.pols,x,y,m);
  end Speel_without_System_Cache;

  procedure Speel_without_System_Cache
              ( p : in Polynomial_Vector; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector ) is
  begin
    for i in p'range loop
      Polynomials.Speel(p(i),x,y(i),yd,wrk);
      for j in yd'range loop
        Ring.Copy(yd(j),m(i,j));
      end loop;
    end loop;
  end Speel_without_System_Cache;

  procedure Speel_without_System_Cache
              ( p : in Polynomial_Vector; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix ) is

    yd,wrk : Vectors.Vector(x'range);

  begin
    Speel_without_System_Cache(p,x,y,m,yd,wrk);
  end Speel_without_System_Cache;

  procedure Speel_without_System_Cache
              ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix;
                yd,wrk : in out Vectors.Vector ) is
  begin
    if p /= null
     then Speel_without_System_Cache(p.all,x,y,m,yd,wrk);
    end if;
  end Speel_without_System_Cache;

  procedure Speel_without_System_Cache
              ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                y : out Vectors.Vector; m : out Matrices.Matrix ) is
  begin
    if p /= null
     then Speel_without_System_Cache(p.all,x,y,m);
    end if;
  end Speel_without_System_Cache;

-- DESTRUCTORS :

  procedure Clear ( p : in out Polynomial_Vector ) is
  begin
    for i in p'range loop
      Polynomials.Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Polynomial_Vector ) is

    procedure free is
      new unchecked_deallocation(Polynomial_Vector,Link_to_Polynomial_Vector);

  begin
    if p /= null then
      Clear(p.all);
      free(p);
    end if;
  end Clear;

  procedure Clear ( s : in out System ) is
  begin
    for i in 1..s.nbr loop
      Polynomials.Clear(s.pols(i));
    end loop;
  end Clear;

  procedure Clear ( s : in out Link_to_System ) is

    procedure free is new unchecked_deallocation(System,Link_to_System);

  begin
    if s /= null then
      Clear(s.all);
      free(s);
    end if;
  end Clear;

end Generic_Polynomial_Vectors;
