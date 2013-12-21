with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Binomial_Varieties;       use Standard_Binomial_Varieties;

package body Standard_Monomial_Map_Solvers is

  function Free_Monomial_Map
              ( n,d : integer32;
                free : Standard_Integer_Vectors.Vector )
              return Monomial_Map is

    res : Monomial_Map(n);
    ind_free : integer32 := 0;

  begin
    res.d := d;
    for i in 1..n loop
      declare
        v : Standard_Integer_Vectors.Vector(1..d) := (1..d => 0);
      begin
        if free(i) = 0 then
          res.c(i) := Create(0.0);
        else
          res.c(i) := Create(1.0);
          ind_free := ind_free + 1;
          v(ind_free) := 1;
        end if;
        res.v(i) := new Standard_Integer_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Free_Monomial_Map;

  function Monomial_Map_Solution
              ( n,d : integer32;
                M : Standard_Integer_Matrices.Matrix; c : Solution )
              return Monomial_Map is

    res : Monomial_Map(n);

  begin
    res.d := d;
    res.c := Transform_Coefficients(d,M,c.v);
    for j in 1..n loop
      declare
        v : Standard_Integer_Vectors.Vector(1..d);
      begin
        for i in 1..d loop
          v(i) := M(i,j);
        end loop;
        res.v(j) := new Standard_Integer_Vectors.Vector'(v);
      end; 
    end loop;
    return res;
  end Monomial_Map_Solution;

  function Monomial_Map_Solution
              ( n,d : integer32;
                s : Standard_Integer_Vectors.Vector;
                M : Standard_Integer_Matrices.Matrix; c : Solution )
              return Monomial_Map is

    res : Monomial_Map(n);
    ind : integer32 := 0;
    tc : constant Standard_Complex_Vectors.Vector
       := Transform_Coefficients(d,M,c.v);

  begin
    res.d := d;
    for j in 1..n loop
      declare
        v : Standard_Integer_Vectors.Vector(1..d);
      begin
        if s(j) = 1 then
          v := (1..d => 0);      
          res.c(j) := Create(0.0);
        else
          ind := ind + 1;
          res.c(j) := tc(ind); 
          for i in 1..d loop
            v(i) := M(i,ind);
          end loop;
        end if;
        res.v(j) := new Standard_Integer_Vectors.Vector'(v);
      end; 
    end loop;
    return res;
  end Monomial_Map_Solution;

  function Monomial_Map_Solution
              ( n,d,e : integer32;
                s,f : Standard_Integer_Vectors.Vector;
                M : Standard_Integer_Matrices.Matrix; c : Solution )
              return Monomial_Map is

    res : Monomial_Map(n);
    ind : integer32 := 0;
    cnt_free : integer32 := 0;
    tc : constant Standard_Complex_Vectors.Vector
       := Transform_Coefficients(d,M,c.v);

  begin
    res.d := d+e;
    for j in 1..n loop
      declare
        v : Standard_Integer_Vectors.Vector(1..d+e);
      begin
        if s(j) = 1 then
          v := (1..d+e => 0);      
          res.c(j) := Create(0.0);
        elsif f(j) = 1 then
          cnt_free := cnt_free + 1;
          v := (1..d+e => 0);
          v(d+cnt_free) := 1;        -- put free variables last
          res.c(j) := Create(1.0);
        else
          ind := ind + 1;
          res.c(j) := tc(ind); 
          for i in 1..d loop
            v(i) := M(i,ind);
          end loop;
          for i in d+1..d+e loop
            v(i) := 0;
          end loop;
        end if;
        res.v(j) := new Standard_Integer_Vectors.Vector'(v);
      end; 
    end loop;
    return res;
  end Monomial_Map_Solution;

  function Monomial_Map_Solutions 
              ( n,d : integer32;
                M : Standard_Integer_Matrices.Matrix; c : Solution_List )
              return Monomial_Map_Array is

    res : Monomial_Map_Array(1..integer32(Length_Of(c)));
    tmp : Solution_List := c;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      declare
        map : constant Monomial_Map(n)
            := Monomial_Map_Solution(n,d,M,ls.all);
      begin
        res(i) := new Monomial_Map'(map);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Monomial_Map_Solutions;

  function Monomial_Map_Solutions 
              ( n,d : integer32;
                s : Standard_Integer_Vectors.Vector;
                M : Standard_Integer_Matrices.Matrix; c : Solution_List )
              return Monomial_Map_Array is

    res : Monomial_Map_Array(1..integer32(Length_Of(c)));
    tmp : Solution_List := c;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      declare
        map : constant Monomial_Map(n)
            := Monomial_Map_Solution(n,d,s,M,ls.all);
      begin
        res(i) := new Monomial_Map'(map);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Monomial_Map_Solutions;

  function Monomial_Map_Solutions 
              ( n,d,e : integer32;
                s,f : Standard_Integer_Vectors.Vector;
                M : Standard_Integer_Matrices.Matrix; c : Solution_List )
              return Monomial_Map_Array is

    res : Monomial_Map_Array(1..integer32(Length_Of(c)));
    tmp : Solution_List := c;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      declare
        map : constant Monomial_Map(n)
            := Monomial_Map_Solution(n,d,e,s,f,M,ls.all);
      begin
        res(i) := new Monomial_Map'(map);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Monomial_Map_Solutions;

  function Toric_Solve ( p : Laur_Sys ) return Link_to_Monomial_Map_Array is

    res : Link_to_Monomial_Map_Array;
    d : integer32;
    M : Standard_Integer_Matrices.Link_to_Matrix;
    c : Solution_List;
    fail : boolean;

  begin
    Black_Box_Solver(p,fail,d,M,c);
    if Length_Of(c) > 0 then
      declare
        maps : Monomial_Map_Array(1..integer32(Length_Of(c)));
      begin
        maps := Monomial_Map_Solutions(M'last(1),d,M.all,c);
        res := new Monomial_Map_Array'(maps);
      end;
    end if;
    return res;
  end Toric_Solve;

  function Affine_Solve
             ( p : Laur_Sys; s : Standard_Integer_Vectors.Vector ) 
             return Link_to_Monomial_Map_Array is

    res : Link_to_Monomial_Map_Array;
    d : integer32;
    M : Standard_Integer_Matrices.Link_to_Matrix;
    c : Solution_List;
    fail : boolean;

  begin
    Black_Box_Solver(p,fail,d,M,c);
    if Length_Of(c) > 0 then
      declare
        maps : Monomial_Map_Array(1..integer32(Length_Of(c)));
      begin
        maps := Monomial_Map_Solutions(s'last,d,s,M.all,c);
        res := new Monomial_Map_Array'(maps);
      end;
    end if;
    return res;
  end Affine_Solve;

  function Affine_Solve
             ( p : Laur_Sys; nbf : integer32;
               s,f : Standard_Integer_Vectors.Vector ) 
             return Link_to_Monomial_Map_Array is

    res : Link_to_Monomial_Map_Array;
    d : integer32;
    M : Standard_Integer_Matrices.Link_to_Matrix;
    c : Solution_List;
    fail : boolean;

  begin
    Black_Box_Solver(p,fail,d,M,c);
    if Length_Of(c) > 0 then
      declare
        maps : Monomial_Map_Array(1..integer32(Length_Of(c)));
      begin
        maps := Monomial_Map_Solutions(s'last,d,nbf,s,f,M.all,c);
        res := new Monomial_Map_Array'(maps);
      end;
    end if;
    return res;
  end Affine_Solve;

end Standard_Monomial_Map_Solvers;
