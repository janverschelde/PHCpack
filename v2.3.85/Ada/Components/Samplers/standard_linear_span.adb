with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Linear_Solvers;   use Standard_Complex_Linear_Solvers;
with Planes_and_Polynomials;            use Planes_and_Polynomials;

package body Standard_Linear_Span is

-- AUXILIARY :

  function Is_In ( v : Standard_Integer_Vectors.Vector; k : integer32 )
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if there is an index i in v such that v(i) = k.

  begin
    for i in v'range loop
      if v(i) = k
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

-- TARGET ROUTINES :

  function Create ( sols : Solution_List ) return Matrix is

    n : constant integer32 := Head_Of(sols).n;
    len : constant integer32 := integer32(Length_Of(sols))-1;
    res : Matrix(1..len,1..n);
    first : constant Link_to_Solution := Head_Of(sols);
    tmp : Solution_List := Tail_Of(sols);
    ptr : Link_to_Solution;

  begin
    for i in 1..len loop
      ptr := Head_Of(tmp);
      for j in 1..n loop
        res(i,j) := ptr.v(j) - first.v(j);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  procedure Rank ( v : in out Matrix;
                   tol : in double_float; r : out natural32 ) is

    col : integer32 := 1;

  begin
    Triangulate(v,tol,v'last(1),v'last(2));
    r := 0;
    for i in v'range(1) loop
       while (col <= v'last(2)) and then (AbsVal(v(i,col)) < tol) loop
         col := col + 1;
       end loop;
       exit when (col > v'last(2));
       r := r + 1;
    end loop;
  end Rank;

  function Pivots ( v : Matrix; tol : double_float; r : natural32 )
                  return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(r));
    col : integer32 := 1;

  begin
    for i in v'range(1) loop
      while (col <= v'last(2)) and then (AbsVal(v(i,col)) < tol) loop
        col := col + 1;
      end loop;
      exit when (col > v'last(2));
      res(i) := col;
    end loop;
    return res;
  end Pivots;

  function Kernel ( v : Matrix; tol : double_float; r : natural32;
                    pivots : Standard_Integer_Vectors.Vector;
                    point : Standard_Complex_Vectors.Vector )
                  return Matrix is

    n : constant integer32 := v'last(2);
    res : Matrix(1..n-integer32(r),0..n);
    mat : Matrix(1..integer32(r),1..integer32(r));
    rhs : Standard_Complex_Vectors.Vector(1..integer32(r));
    ipvt : Standard_Integer_Vectors.Vector(1..integer32(r));
    info : integer32;
    ind : integer32 := 0;
    eva : Complex_Number;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        mat(i,j) := v(i,pivots(j));
      end loop;
    end loop;
    lufac(mat,integer32(r),ipvt,info);
    for j in 1..n loop
      if not Is_In(pivots,j) then
        ind := ind + 1;
        res(ind,j) := Create(1.0);
        for i in 1..integer32(r) loop
          rhs(i) := -v(i,j);
        end loop;
        lusolve(mat,integer32(r),ipvt,rhs);
        for i in 1..integer32(r) loop
          res(ind,pivots(i)) := rhs(i);
        end loop;
        eva := Create(0.0);
        for k in point'range loop
          eva := eva + point(k)*res(ind,k);
        end loop;
        res(ind,0) := -eva;
      end if;
    end loop;
    for i in res'range(1) loop
      for j in res'range(2) loop
        if AbsVal(res(i,j)) < tol
         then res(i,j) := Create(0.0);
        end if;
      end loop;
    end loop;
    return res;
  end Kernel;

  function Equations ( m : Matrix ) return Poly_Sys is

    res : Poly_Sys(m'range(1));
    h : Standard_Complex_Vectors.Vector(0..m'last(2));

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        h(j) := m(i,j);
      end loop;
      res(i) := Hyperplane(h);
    end loop;
    return res;
  end Equations;

  function Eliminators ( kernel_eqs : Matrix;
                         pivots : Standard_Integer_Vectors.Vector )
                       return Poly_Sys is

    res : Poly_Sys(kernel_eqs'range(1));
    hyp : Standard_Complex_Vectors.Vector(0..pivots'last);

  begin
    for i in kernel_eqs'range(1) loop
      hyp(0) := -kernel_eqs(i,0);
      for j in pivots'range loop
        hyp(j) := -kernel_eqs(i,pivots(j));
      end loop;
      res(i) := Hyperplane(hyp);
    end loop;
    return res;
  end Eliminators;

  function Eliminate_non_Pivots
             ( t : Term; pivots : Standard_Integer_Vectors.Vector;
               elm : Poly_Sys ) return Poly is

    res : Poly;
    rt : Term;
    piv_ind : integer32 := 0;
    elm_ind : integer32 := 0;
    first : boolean := true;
    fac : Poly := Null_Poly;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Natural_Vectors.Vector(pivots'range);
    for i in t.dg'range loop
      if Is_In(pivots,i) then
        piv_ind := piv_ind + 1;
        rt.dg(piv_ind) := t.dg(i);
      else
        elm_ind := elm_ind + 1;
        if t.dg(i) > 0 then
          if first then
            first := false;
            Copy(elm(elm_ind),fac);
            for j in 2..t.dg(i) loop
              Mul(fac,elm(elm_ind));
            end loop;
          else
            for j in 1..t.dg(i) loop
              Mul(fac,elm(elm_ind));
            end loop;
          end if;
        end if;
      end if;
    end loop;
    if first
     then res := Create(rt);
     else res := rt*fac;
    end if;
    Clear(rt.dg); Clear(fac);
    return res;
  end Eliminate_non_Pivots;

  function Eliminate_non_Pivots
             ( p : Poly; pivots : Standard_Integer_Vectors.Vector;
               elm : Poly_Sys ) return Poly is

    res : Poly := Null_Poly;

    procedure Eliminate_in_Term ( t : in Term; continue : out boolean ) is

      et : Poly := Eliminate_non_Pivots(t,pivots,elm);
      
    begin
      Add(res,et);
      Standard_Complex_Polynomials.Clear(et);
      continue := true;
    end Eliminate_in_Term;
    procedure Eliminate_in_Terms is
      new Visiting_Iterator(Eliminate_in_Term);

  begin
    Eliminate_in_Terms(p);
    return res;
  end Eliminate_non_Pivots;

  function Eliminate_non_Pivots
             ( p : Poly_Sys; pivots : Standard_Integer_Vectors.Vector;
               elm : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Eliminate_non_Pivots(p(i),pivots,elm);
    end loop;
    return res;
  end Eliminate_non_Pivots;

  function Filter ( p : Poly; tol : double_float ) return Poly is

    res : Poly := Null_Poly;

    procedure Filter_Term ( t : in Term; continue : out boolean ) is
    begin
      if AbsVal(t.cf) > tol
       then Add(res,t);
      end if;
      continue := true;
    end Filter_Term;
    procedure Filter_Terms is new Visiting_Iterator(Filter_Term);

  begin
    Filter_Terms(p);
    return res;
  end Filter;

  function Filter ( p : Poly_Sys; tol : double_float ) return Poly_Sys is

    res : Poly_Sys(p'range);
    cnt : integer32 := 0;

  begin
    for i in p'range loop
      declare
        fp : constant Poly := Filter(p(i),tol);
      begin
        if fp /= Null_Poly then
          cnt := cnt + 1;
          res(cnt) := fp;
        end if;
      end;
    end loop;
    return res(res'first..cnt);
  end Filter;

  function Eliminate_non_Pivots
             ( v : Standard_Complex_Vectors.Vector;
               pivots : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(pivots'range);

  begin
    for i in pivots'range loop
      res(i) := v(pivots(i));
    end loop;
    return res;
  end Eliminate_non_Pivots;

  function Eliminate_non_Pivots
             ( s : Solution; pivots : Standard_Integer_Vectors.Vector )
             return Solution is

    res : Solution(pivots'last);

  begin
    res.t := s.t;
    res.m := s.m;
    res.v := Eliminate_non_Pivots(s.v,pivots);
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Eliminate_non_Pivots;

  function Eliminate_non_Pivots
             ( sols : Solution_List;
               pivots : Standard_Integer_Vectors.Vector )
             return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Append(res,res_last,Eliminate_non_Pivots(ls.all,pivots));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Eliminate_non_Pivots;

end Standard_Linear_Span;
