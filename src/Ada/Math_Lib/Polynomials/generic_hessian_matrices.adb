with unchecked_deallocation;

package body Generic_Hessian_Matrices is 

-- CREATORS :

  function Create ( p : Poly ) return Hessian is

    dim : constant integer32 := integer32(Number_of_Unknowns(p));
    res : Hessian(1..dim,1..dim);
    dpi : Poly;

  begin
    for i in 1..dim loop
      dpi := Diff(p,i);
      for j in 1..(i-1) loop
        Copy(res(j,i),res(i,j));
      end loop;
      for j in i..dim loop
        res(i,j) := Diff(dpi,j);
      end loop;
      Clear(dpi);
    end loop;
    return res;
  end Create;

  function Create ( p : Poly ) return Link_to_Hessian is

    res : constant Link_to_Hessian := new Hessian'(Create(p));

  begin
    return res;   
  end Create;

  function Create ( p : Poly_Sys ) return Array_of_Hessians is

    res : Array_of_Hessians(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Create;

  function Create ( p : Poly; k : integer32 ) return Hessian is

    dim : constant integer32 := integer32(Number_of_Unknowns(p));
    dm1 : constant integer32 := dim-1;
    res : Hessian(1..dm1,1..dm1);
    dpi : Poly;

  begin
    for i in 1..dim loop
      if i < k then
        dpi := Diff(p,i);
        for j in 1..(i-1) loop
          if j < k then
            Copy(res(j,i),res(i,j));
          elsif j > k then
            Copy(res(j-1,i),res(i,j-1));
          end if;
        end loop;
        for j in i..dim loop
          if j < k then
            res(i,j) := Diff(dpi,j);
          elsif j > k then
            res(i,j-1) := Diff(dpi,j);
          end if;
        end loop;
      elsif i > k then
        dpi := Diff(p,i);
        for j in 1..(i-1) loop
          if j < k then
            Copy(res(j,i-1),res(i-1,j));
          elsif j > k then
            Copy(res(j-1,i-1),res(i-1,j-1));
          end if;
        end loop;
        for j in i..dim loop
          if j < k then
            res(i-1,j) := Diff(dpi,j);
          elsif j > k then
            res(i-1,j-1) := Diff(dpi,j);
          end if;
        end loop;
      end if;
      Clear(dpi);
    end loop;
    return res;
  end Create;

  function Create ( p : Poly; k : integer32 ) return Link_to_Hessian is

    res : constant Link_to_Hessian := new Hessian'(Create(p,k));

  begin
    return res;
  end Create;

  function Create ( p : Poly_Sys; k : integer32 ) return Array_of_Hessians is

    res : Array_of_Hessians(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i),k);
    end loop;
    return res;
  end Create;

-- EVALUATORS :

  function Eval ( h : Hessian; x : Vector ) return Matrix is

    res : Matrix(h'range(1),h'range(2));

  begin
    for i in h'range(1) loop
      for j in h'range(2) loop
        res(i,j) := zero;
      end loop;
    end loop;
    for i in h'range(1) loop
      for j in h'first(2)..(i-1) loop
        res(i,j) := res(j,i);
      end loop;
      for j in i..h'last(1) loop
        res(i,j) := Eval(h(i,j),x);
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( h : Link_to_Hessian; x : Vector ) return Matrix is

    res : Matrix(h'range(1),h'range(2));

  begin
    for i in h'range(1) loop
      for j in h'range(2) loop
        res(i,j) := zero;
      end loop;
    end loop;
    for i in h'range(1) loop
      for j in h'first(2)..(i-1) loop
        res(i,j) := res(j,i);
      end loop;
      for j in i..h'last(1) loop
        res(i,j) := Eval(h(i,j),x);
      end loop;
    end loop;
    return res;
  end Eval;
  
-- DESTRUCTORS :

  procedure Clear ( h : in out Hessian ) is
  begin
    for i in h'range(1) loop
      for j in h'range(2) loop
        Clear(h(i,j));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( h : in out Link_to_Hessian ) is

    procedure free is
      new unchecked_deallocation(Hessian,Link_to_Hessian);

  begin
    if h /= null then
      Clear(h.all);
      free(h);
    end if;
  end Clear;

  procedure Clear ( h : in out Array_of_Hessians ) is
  begin
    for i in h'range loop
      Clear(h(i));
    end loop;
  end Clear;

  procedure Clear ( h : in out Link_to_Array_of_Hessians ) is

    procedure free is
      new unchecked_deallocation
            (Array_of_Hessians,Link_to_Array_of_Hessians);

  begin
    if h /= null then
      Clear(h.all);
      free(h);
    end if;
  end Clear;

end Generic_Hessian_Matrices;
