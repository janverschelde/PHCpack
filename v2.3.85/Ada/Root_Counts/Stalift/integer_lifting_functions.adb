with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Random_Numbers;

package body Integer_Lifting_Functions is

-- AUXILIARIES :

  function Convert ( v : Standard_Integer_Vectors.Vector )
                   return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the given vector into a vector with complex entries.

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Create(double_float(v(i)));
    end loop;
    return res;
  end Convert;

  function Random_Vector ( m,low,upp : integer32 ) return Vector is

  -- DESCRIPTION :
  --   Returns a random vector of range 1..m of randomly generated integer
  --   values between low and upp.

    res : Vector(1..m);

  begin
    for k in res'range loop
      res(k) := Standard_Random_Numbers.Random(low,upp);
    end loop;
    return res;
  end Random_Vector;

  function Random ( vec : Standard_Integer_Vectors.Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns a random number from the given vector.

    index : constant integer32
          := Standard_Random_Numbers.Random(vec'first,vec'last);

  begin
    return vec(index);
  end Random;

-- LINEAR LIFTING :

  function Linear_Lift ( lf,v : Vector ) return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := lf*v;
    return res;
  end Linear_Lift;

  function Linear_Lift ( lf : Vector; L : List ) return List is

    res,res_last : List;
    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Linear_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Linear_Lift;

  function Linear_Lift ( lf : VecVec;
                         L : Array_of_Lists ) return Array_of_Lists is

    res : Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Linear_Lift(lf(i).all,L(i));
    end loop;
    return res;
  end Linear_Lift;   

-- POLYNOMIAL LIFTING :

  function Polynomial_Lift
             ( lf : Standard_Complex_Polynomials.Poly;
               v : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    use Standard_Complex_Poly_Functions;

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := integer32(REAL_PART(Eval(lf,Convert(v))));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Laurentials.Poly;
               v : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    use Standard_Complex_Laur_Functions;

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := integer32(REAL_PART(Eval(lf,Convert(v))));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Poly_Functions.Eval_Poly;
               v : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    use Standard_Complex_Poly_Functions;

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := integer32(REAL_PART(Eval(lf,Convert(v))));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Laur_Functions.Eval_Poly;
               v : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    use Standard_Complex_Laur_Functions;

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := integer32(REAL_PART(Eval(lf,Convert(v))));
    return res;
  end Polynomial_Lift;
 
  function Polynomial_Lift
             ( lf : Standard_Complex_Polynomials.Poly;
               L : List ) return List is
 
    res : List;
    elf : Standard_Complex_Poly_Functions.Eval_Poly
        := Standard_Complex_Poly_Functions.Create(lf);

  begin
    res := Polynomial_Lift(elf,L);
    Standard_Complex_Poly_Functions.Clear(elf);
    return res;
  end Polynomial_Lift;
 
  function Polynomial_Lift
             ( lf : Standard_Complex_Laurentials.Poly;
               L : List ) return List is
 
    res : List;
    elf : Standard_Complex_Laur_Functions.Eval_Poly
        := Standard_Complex_Laur_Functions.Create(lf);

  begin
    res := Polynomial_Lift(elf,L);
    Standard_Complex_Laur_Functions.Clear(elf);
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Standard_Complex_Poly_Functions.Eval_Poly;
               L : List ) return List is

    res,res_last : List;
    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift 
             ( lf : Standard_Complex_Laur_Functions.Eval_Poly;
               L : List ) return List is

    res,res_last : List;
    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Poly_Sys; L : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(L'range);
    elf : Eval_Poly_Sys(lf'range) := Create(lf);
  
  begin
    res := Polynomial_Lift(elf,L);
    Standard_Complex_Poly_SysFun.Clear(elf);
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Laur_Sys; L : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(L'range);
    elf : Eval_Laur_Sys(lf'range) := Create(lf);
  
  begin
    res := Polynomial_Lift(elf,L);
    Standard_Complex_Laur_SysFun.Clear(elf);
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Eval_Poly_Sys; 
               L : Array_of_Lists ) return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Polynomial_Lift(lf(i),L(i));
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift
             ( lf : Eval_Laur_Sys; 
               L : Array_of_Lists ) return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Polynomial_Lift(lf(i),L(i));
    end loop;
    return res;
  end Polynomial_Lift;

-- RANDOM LIFTING :

  function Random_Lift ( lflow,lfupp : integer32; v : Vector ) return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := Standard_Random_Numbers.Random(lflow,lfupp);
    return res;
  end Random_Lift;

  function Random_Lift ( lflow,lfupp : integer32; L : List ) return List is

    res,res_last : List;
    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Random_Lift(lflow,lfupp,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Random_Lift;

  function Random_Lift ( lflow,lfupp : Vector; L : Array_of_Lists )
                       return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Random_Lift(lflow(i),lfupp(i),L(i));
    end loop;
    return res;
  end Random_Lift;

-- RANDOM LINEAR LIFTING :

  function Random_Linear_Lift ( lflow,lfupp : integer32; v : Vector )
                              return Vector is

    lf : constant Vector(v'range) := Random_Vector(v'last,lflow,lfupp);
 
  begin
    return Linear_Lift(lf,v);
  end Random_Linear_Lift;

  function Random_Linear_Lift ( lflow,lfupp : integer32; L : List )
                              return List is
  begin
    if Is_Null(L) then
      return L;
    else
      declare
        n : constant integer32 := Head_Of(L)'length;
        lf : constant Vector(Head_Of(l)'range) := Random_Vector(n,lflow,lfupp);
      begin
        return Linear_Lift(lf,L);
      end;
    end if;
  end Random_Linear_Lift;

  function Random_Linear_Lift ( lflow,lfupp : Vector; L : Array_of_Lists )
                              return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Random_Linear_Lift(lflow(i),lfupp(i),L(i));
    end loop;
    return res;
  end Random_Linear_Lift;

-- POINT-WISE LIFTING :

  function Point_Lift ( lf : integer32; v : Vector ) return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := lf;
    return res;
  end Point_Lift;

  function Point_Lift ( lf : Vector; L : List ) return List is

    res,res_last : List;
    tmp : List := L;
    ind : integer32 := lf'first;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Point_Lift(lf(ind),Head_Of(tmp).all));
      ind := ind + 1;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Point_Lift;

  function Point_Lift ( lf : VecVec;
                        L : Array_of_Lists ) return Array_of_Lists is

    res : Array_of_Lists(L'range);

  begin
    for i in res'range loop
      res(i) := Point_Lift(lf(i).all,L(i));
    end loop;
    return res;
  end Point_Lift;

end Integer_Lifting_Functions;
