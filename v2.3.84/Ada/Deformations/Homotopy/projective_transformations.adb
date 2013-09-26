with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;

package body Projective_Transformations is

  function Projective_Transformation ( p : Poly ) return Poly is
  
    deg : constant integer32 := Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := natural32(deg)-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation ( p : in out Poly ) is
  
    res : Poly := Projective_Transformation(p);

  begin
    Copy(res,p); Clear(res);
  end Projective_Transformation;

  function Projective_Transformation ( p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation ( p : in out Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  function Projective_Transformation ( s : Solution ) return Solution is

    n : constant integer32 := s.n;
    r : Solution(n+1);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(1.0);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : Solution_List ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation ( sols : in out Solution_List ) is
  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  function Affine_Transformation ( s : Solution ) return Solution is

    n : constant integer32 := s.n;
    r : Solution(n-1);
    absvn : constant double_float := AbsVal(s.v(n));

  begin
    for i in 1..(n-1) loop
      if absvn + Create(1.0) = Create(1.0)
       then r.v(i) := Create(10.0**10);
       else r.v(i) := s.v(i) / s.v(n);
      end if;
     end loop;
     r.t := s.t;
     r.m := s.m;
     r.err := s.err;
     r.rco := s.rco;
     r.res := s.res;
     return r;
  exception
    when constraint_error =>
       r.v(1..(n-1)) := (1..(n-1) => Create(10.0**10));
       return r;
  end Affine_Transformation;

  procedure Affine_Transformation ( sols : in out Solution_List ) is
  begin
    if Is_Null(sols) then
      null;
    else
      declare
        n : constant integer32 := Head_Of(sols).n;
        s1 : Solution(n);
        s2 : Solution(n-1);
        temp : Solution_List := sols;
        ls : Link_To_Solution;
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s1 := ls.all;
          s2 := Affine_Transformation(s1);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Affine_Transformation;

end Projective_Transformations;
