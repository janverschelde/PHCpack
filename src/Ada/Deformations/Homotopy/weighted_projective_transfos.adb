with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Weighted_Projective_Transfos is

  function Projective_Transformation ( p : Poly; q : Vector ) return Poly is
  
    deg : Degrees := new Vector(q'first..q'last+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural := 0;

    begin
      ht.cf := t.cf;
      deg(t.dg'range) := t.dg.all;
      deg(deg'last) := t.dg.all*q;
      ht.dg := deg;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(deg);
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation ( p : in out Poly; q : in Vector ) is
  
    res : Poly := Projective_Transformation(p,q);

  begin
    Copy(res,p); Clear(res);
  end Projective_Transformation;

  function Projective_Transformation
             ( p : Laur_Sys; q : Vector ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k),q);
    end loop;
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation ( p : in out Laur_Sys; q : in Vector ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k),q);
    end loop;
  end Projective_Transformation;

  function Affine_Transformation ( s : Solution; q : Vector ) return Solution is

    n : natural := s.n;
    r : Solution(n-1);
    absvn : double_float := AbsVal(s.v(n));

  begin
    for i in 1..(n-1) loop
      r.v(i) := s.v(i);
      for j in 1..q(i) loop
        Mul(r.v(i),s.v(n));
      end loop;
     end loop;
     r.t := s.t;
     r.m := s.m;
     r.err := s.err;
     r.rco := s.rco;
     r.res := s.res;
     return r;
  exception
    when numeric_error => r.v(1..(n-1)) := (1..(n-1) => Create(10.0**10));
                          return r;
  end Affine_Transformation;

  function Affine_Transformation
             ( sols : Solution_List; q : Vector ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Affine_Transformation(ls.all,q));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Affine_Transformation;

  procedure Affine_Transformation
              ( sols : in out Solution_List; q : in Vector ) is
  begin
    if Is_Null(sols)
     then null;
     else declare
            n : natural := Head_Of(sols).n;
            s1 : Solution(n);
            s2 : Solution(n-1);
            temp : Solution_List := sols;
            l : Link_To_Solution;
          begin
            while not Is_Null(temp) loop
              l := Head_Of(temp);
              s1 := l.all;
              s2 := Affine_Transformation(s1,q);
              Clear(l);
              l := new Solution'(s2);
              Set_Head(temp,l);
              temp := Tail_Of(temp);
            end loop;
          end;
    end if;
  end Affine_Transformation;

end Weighted_Projective_Transfos;
