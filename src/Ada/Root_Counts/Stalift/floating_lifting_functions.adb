with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Floating_VecVecs;          use Standard_Floating_VecVecs;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;

package body Floating_Lifting_Functions is

-- AUXILIARIES :

  function Flt2Cmplx ( x : Standard_Floating_Vectors.Vector )
                     return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with complex entries.

    res : Standard_Complex_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := Create(x(i));
    end loop;
    return res;
  end Flt2Cmplx;

-- RANDOM FLOATING-POINT LIFTING :

  function Random_Lift ( lflow,lfupp : double_float ) return double_float is

    res : double_float := random;                          -- in [-1,1]

  begin
    res := ((1.0+res)/2.0)*lflow + ((1.0-res)/2.0)*lfupp;  -- in [lflow,lfupp]
    return res;
  end Random_Lift;

  function Random_Lift ( v : Vector; lflow,lfupp : double_float )
                       return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v; 
    res(res'last) := Random_Lift(lflow,lfupp);
    return res;
  end Random_Lift;

  function Random_Lift ( L : List; lflow,lfupp : double_float ) return List is

    res,res_last : List;
    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Random_Lift(Head_Of(tmp).all,lflow,lfupp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Random_Lift;

  function Random_Lift ( L : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                         lflow,lfupp : Vector )
                       return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Random_Lift(L(i),lflow(i),lfupp(i));
    end loop;
    return res;
  end Random_Lift;

-- LINEAR LIFTING FUNCTIONS :

  function Linear_Lift ( x,v : Vector ) return Vector is

    res : Vector(x'first..x'last+1);

  begin
    res(x'range) := x;
    res(res'last) := x*v;
    return res;
  end Linear_Lift;

  function Linear_Lift ( f : Face; v : Vector ) return Face is

    res : constant Face := new VecVec(f'range);

  begin
    for i in res'range loop
      res(i) := new Standard_Floating_Vectors.Vector'(Linear_Lift(f(i).all,v));
    end loop;
    return res;
  end Linear_Lift;

  function Linear_Lift ( L : List; v : Vector ) return List is

  -- DESCRIPTION :
  --   Returns a linearly lifted list of points.

    res,res_last : List;
    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Linear_Lift(Head_Of(tmp).all,v));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Linear_Lift;

  function Linear_Lift ( f : Faces; v : Vector ) return Faces is

    res,res_last : Faces;
    tmp : Faces := f;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Linear_Lift(Head_Of(tmp),v));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Linear_Lift;

-- RANDOM FLOATING-POINT LINEAR LIFTING :

  function Random ( n : integer32; lflow,lfupp : double_float )
                  return Vector is

    res : Vector(1..n);

  begin
    for i in res'range loop
      res(i) := Random_Lift(lflow,lfupp);
    end loop;
    return res;
  end Random;

-- POLYNOMIAL LIFTING FUNCTIONS :

  function Polynomial_Lift
             ( lf : Standard_Complex_Polynomials.Poly;
               x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

    use Standard_Complex_Poly_Functions;

    res : Vector(x'first..x'last+1);

  begin
    res(x'range) := x;
    res(res'last) := REAL_PART(Eval(lf,Flt2Cmplx(x)));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Laurentials.Poly;
               x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

    use Standard_Complex_Laur_Functions;

    res : Vector(x'first..x'last+1);

  begin
    res(x'range) := x;
    res(res'last) := REAL_PART(Eval(lf,Flt2Cmplx(x)));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Poly_Functions.Eval_Poly;
               x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

    use Standard_Complex_Poly_Functions;

    res : Vector(x'first..x'last+1);

  begin
    res(x'range) := x;
    res(res'last) := REAL_PART(Eval(lf,Flt2Cmplx(x)));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Laur_Functions.Eval_Poly;
               x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

    use Standard_Complex_Laur_Functions;

    res : Vector(x'first..x'last+1);

  begin
    res(x'range) := x;
    res(res'last) := REAL_PART(Eval(lf,Flt2Cmplx(x)));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Polynomials.Poly;
               L : List ) return List is

    res,res_last,tmp : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Laurentials.Poly;
               L : List ) return List is

    res,res_last,tmp : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Poly_Functions.Eval_Poly;
               L : List ) return List is

    res,res_last,tmp : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Laur_Functions.Eval_Poly;
               L : List ) return List is

    res,res_last,tmp : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Standard_Complex_Poly_Systems.Poly_Sys;
                             L : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Polynomial_Lift(lf(i),L(i));
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Standard_Complex_Laur_Systems.Laur_Sys;
                             L : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Polynomial_Lift(lf(i),L(i));
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Eval_Poly_Sys; L : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Polynomial_Lift(lf(i),L(i));
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Eval_Laur_Sys; L : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Polynomial_Lift(lf(i),L(i));
    end loop;
    return res;
  end Polynomial_Lift;

-- FOR STABLE MIXED VOLUMES :

  function Max_Degree ( p : in Standard_Complex_Poly_Systems.Poly_Sys )
                      return integer32 is

    use Standard_Complex_Polynomials;

    max : integer32 := Degree(p(p'first));
    deg : integer32;

  begin
    for i in p'first+1..p'last loop
      deg := Degree(p(i));
      if deg > max
       then max := deg;
      end if;
    end loop;
    return max;
  end Max_Degree;

  function Max_Degree ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys )
                      return integer32 is

    use DoblDobl_Complex_Polynomials;

    max : integer32 := Degree(p(p'first));
    deg : integer32;

  begin
    for i in p'first+1..p'last loop
      deg := Degree(p(i));
      if deg > max
       then max := deg;
      end if;
    end loop;
    return max;
  end Max_Degree;

  function Max_Degree ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys )
                      return integer32 is

    use QuadDobl_Complex_Polynomials;

    max : integer32 := Degree(p(p'first));
    deg : integer32;

  begin
    for i in p'first+1..p'last loop
      deg := Degree(p(i));
      if deg > max
       then max := deg;
      end if;
    end loop;
    return max;
  end Max_Degree;

  function Max_Degree ( p : in Standard_Complex_Laur_Systems.Laur_Sys )
                      return integer32 is

    use Standard_Complex_Laurentials;

    max : integer32 := Degree(p(p'first));
    deg : integer32;

  begin
    for i in p'first+1..p'last loop
      deg := Degree(p(i));
      if deg > max
       then max := deg;
      end if;
    end loop;
    return max;
  end Max_Degree;

  function Max_Degree ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys )
                      return integer32 is

    use DoblDobl_Complex_Laurentials;

    max : integer32 := Degree(p(p'first));
    deg : integer32;

  begin
    for i in p'first+1..p'last loop
      deg := Degree(p(i));
      if deg > max
       then max := deg;
      end if;
    end loop;
    return max;
  end Max_Degree;

  function Max_Degree ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys )
                      return integer32 is

    use QuadDobl_Complex_Laurentials;

    max : integer32 := Degree(p(p'first));
    deg : integer32;

  begin
    for i in p'first+1..p'last loop
      deg := Degree(p(i));
      if deg > max
       then max := deg;
      end if;
    end loop;
    return max;
  end Max_Degree;

  function Lifting_Bound ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float is

    use Standard_Complex_Polynomials;

    res : double_float;
    n : constant integer32 := p'last;
    d : constant integer32 := Max_Degree(p);

  begin
    res := double_float(n*(n+1));
    for i in 1..n loop
      res := res*double_float(d);
      if res > max
       then return max;
      end if;
    end loop;
    return res;
  end Lifting_Bound;

  function Lifting_Bound ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float is

    use DoblDobl_Complex_Polynomials;

    res : double_float;
    n : constant integer32 := p'last;
    d : constant integer32 := Max_Degree(p);

  begin
    res := double_float(n*(n+1));
    for i in 1..n loop
      res := res*double_float(d);
      if res > max
       then return max;
      end if;
    end loop;
    return res;
  end Lifting_Bound;

  function Lifting_Bound ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float is

    use QuadDobl_Complex_Polynomials;

    res : double_float;
    n : constant integer32 := p'last;
    d : constant integer32 := Max_Degree(p);

  begin
    res := double_float(n*(n+1));
    for i in 1..n loop
      res := res*double_float(d);
      if res > max
       then return max;
      end if;
    end loop;
    return res;
  end Lifting_Bound;

  function Lifting_Bound ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float is

    use Standard_Complex_Laurentials;

    res : double_float;
    n : constant integer32 := p'last;
    d : constant integer32 := Max_Degree(p);

  begin
    res := double_float(n*(n+1));
    for i in 1..n loop
      res := res*double_float(d);
      if res > max
       then return max;
      end if;
    end loop;
    return res;
  end Lifting_Bound;

  function Lifting_Bound ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float is

    use DoblDobl_Complex_Laurentials;

    res : double_float;
    n : constant integer32 := p'last;
    d : constant integer32 := Max_Degree(p);

  begin
    res := double_float(n*(n+1));
    for i in 1..n loop
      res := res*double_float(d);
      if res > max
       then return max;
      end if;
    end loop;
    return res;
  end Lifting_Bound;

  function Lifting_Bound ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float is

    use QuadDobl_Complex_Laurentials;

    res : double_float;
    n : constant integer32 := p'last;
    d : constant integer32 := Max_Degree(p);

  begin
    res := double_float(n*(n+1));
    for i in 1..n loop
      res := res*double_float(d);
      if res > max
       then return max;
      end if;
    end loop;
    return res;
  end Lifting_Bound;

  function Stable_Lift ( L : List; b : double_float ) return List is

    res : List := Random_Lift(L,0.0,1.0);
    n : constant integer32 := Head_Of(L)'last;
    origin : constant Vector(1..n) := (1..n => 0.0);
    lifted : Vector(1..n+1);
    lpt : Link_to_Vector;

  begin
    if not Is_In(L,origin) then
      lifted(1..n) := (1..n => 0.0);
      lifted(n+1) := b;
      lpt := new Vector'(lifted);
      Construct(lpt,res);
    end if;
    return res;
  end Stable_Lift;

  function Stable_Lift ( L : Array_of_Lists; b : double_float )
                       return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Stable_Lift(L(i),b);
    end loop;
    return res;
  end Stable_Lift;

  function Is_Origin ( lpt : Link_to_Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all except the last value of the lifted point
  --   are equal to zero.

  begin
    for i in lpt'first..lpt'last-1 loop
      if lpt(i) /= 0.0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Origin;

  function Lifting_Bound
              ( lifted : Arrays_of_Floating_Vector_Lists.Array_of_Lists ) 
              return double_float is

  -- DESCRIPTION :
  --   Returns the maximal lifting value of the origins.

    res : double_float := 0.0;
    tmp : List;
    lpt : Link_to_Vector;

  begin
    for i in lifted'range loop
      tmp := lifted(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        if Is_Origin(lpt) then
          if lpt(lpt'last) > res
           then res := lpt(lpt'last);
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Lifting_Bound;

end Floating_Lifting_Functions;
