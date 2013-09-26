-- with text_io,integer_io;  use text_io,integer_io;
-- with Standard_Integer_Vectors_io; use Standard_Integer_Vectors_io;
-- with Standard_Integer_Matrices_io; use Standard_Integer_Matrices_io;

with Standard_Common_Divisors;

package body Standard_Monomial_Map_Ideals is

  function One_Variable_Equation ( n,i : integer32 ) return Poly is

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    t.dg(i) := 1;
    t.cf := Create(1.0);
    res := Create(t);
    Clear(t);
    return res;
  end One_Variable_Equation;

  function Is_Constant ( map : Monomial_Map; i : integer32 ) return boolean is

    lv : constant Standard_Integer_Vectors.Link_to_Vector := map.v(i);

  begin
    for j in lv'range loop
      if lv(j) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Constant;

  function One_Variable_Constant_Equation
             ( n,i : integer32; c : Complex_Number ) return Poly is

    res : Poly := One_Variable_Equation(n,i);
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    t.cf := c;
    Sub(res,t);
    Clear(t);
    return res;
  end One_Variable_Constant_Equation;

  function Variable_Difference_Equation
              ( n,i,j : integer32; ci,cj : Complex_Number ) return Poly is

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    t.dg(i) := 1;
    t.cf := cj;
    res := Create(t);
    t.dg(i) := 0;
    t.dg(j) := 1;
    t.cf := ci;
    Sub(res,t);
    Clear(t);
    return res;
  end Variable_Difference_Equation;

  function Same_Exponents
              ( map : Monomial_Map; i,j : integer32 ) return boolean is

    ilv : constant Standard_Integer_Vectors.Link_to_Vector := map.v(i);
    jlv : constant Standard_Integer_Vectors.Link_to_Vector := map.v(j);

  begin
    for k in ilv'range loop
      if ilv(k) /= jlv(k)
       then return false;
      end if;
    end loop;
    return true;
  end Same_Exponents;

  function Variable_Differences ( map : Monomial_Map ) return Poly_List is

    res,res_last : Poly_List;
    n : constant integer32 := map.n;
    done : array(1..n,1..n) of boolean;

  begin
    for i in done'range(1) loop
      for j in i..done'last(2) loop
        done(i,j) := false;
      end loop;
    end loop;
    for i in 1..n loop
      if not Is_Zero(map.c(i)) and not Is_Constant(map,i) then
        for j in (i+1)..n loop
          if not done(i,j) then
            if not Is_Zero(map.c(j)) and not Is_Constant(map,j) then
              if Same_Exponents(map,i,j) then
                declare
                  p : constant Poly
                    := Variable_Difference_Equation(n,i,j,map.c(i),map.c(j));
                begin
                  Append(res,res_last,p);
                  done(i,j) := true;
                  for k in (i+1)..(j-1) loop  -- apply transitivity
                    if done(i,k)
                     then done(k,j) := true;
                    end if;
                  end loop;
                end;
              end if;
            end if;
          end if;
        end loop;
      end if;
    end loop;
    return res;
  end Variable_Differences;

  function Share_One_Same_Parameter
              ( map : Monomial_Map; i,j : integer32 ) return integer32 is

    res : integer32 := 0;
    ilv : constant Standard_Integer_Vectors.Link_to_Vector := map.v(i);
    jlv : constant Standard_Integer_Vectors.Link_to_Vector := map.v(j);

  begin
   -- put("ilv : "); put(ilv); new_line;
   -- put("jlv : "); put(jlv); new_line;
    for k in 1..map.d loop
      if ilv(k) /= jlv(k) then
        if res = 0
         then res := k;      -- first common parameter
         else return 0;      -- second common parameter
        end if;
      end if;
    end loop;
    return res;
  end Share_One_Same_Parameter;

  function Variable_Relation_Same_Sign
              ( n,i,j,a,b,e : integer32;
                c,d : Complex_Number ) return Poly is

    res : Poly;
    ka : constant integer32 := e/a;
    kb : constant integer32 := e/b;
    t1,t2 : Term;

  begin
    t1.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    t1.dg(i) := kb;
    t1.cf := d;
    for k in 2..kb loop
      t1.cf := t1.cf*d;
    end loop;
    t2.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    t2.dg(j) := ka;
    t2.cf := c;
    for k in 2..ka loop
      t2.cf := t2.cf*c;
    end loop;
    t2.cf := t2.cf/t1.cf;  -- normalize the equation
    t1.cf := Create(1.0);
    if t1.dg(i) < 0 then
      t1.dg(i) := -t1.dg(i);
      t2.dg(j) := -t2.dg(j);
      t2.cf := Create(1.0)/t2.cf;
    end if;
    res := Create(t1); Sub(res,t2);
    Clear(t1); Clear(t2);
    return res;
  end Variable_Relation_Same_Sign;

  function Variable_Relation_Opposite_Sign
              ( n,i,j,a,b,e : integer32;
                c,d : Complex_Number ) return Poly is

    res : Poly;
    ka : constant integer32 := e/a;
    kb : constant integer32 := e/b;
    t1,t2 : Term;

  begin
    t1.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    t1.dg(i) := abs(ka);
    t1.dg(j) := abs(kb);
    t1.cf := Create(1.0);
    t2.dg := new Standard_Integer_Vectors.Vector'(1..n => 0);
    t2.cf := c;
    for k in 2..t1.dg(i) loop
      t2.cf := t2.cf*c;
    end loop;
    for k in 1..t1.dg(j) loop
      t2.cf := t2.cf*d;
    end loop;
    res := Create(t1); Sub(res,t2);
    Clear(t1); Clear(t2);
    return res;
  end Variable_Relation_Opposite_Sign;

  function Variable_Relation
              ( n,i,j,a,b : integer32;
                c,d : Complex_Number ) return Poly is

    e : constant integer32 := Standard_Common_Divisors.LCM(abs(a),abs(b));

  begin
    if a*b > 0
     then return Variable_Relation_Same_Sign(n,i,j,a,b,e,c,d);
     else return Variable_Relation_Opposite_Sign(n,i,j,a,b,e,c,d);
    end if;
  end Variable_Relation;

  function Paired_Variables ( map : Monomial_Map ) return Poly_List is

    res,res_last : Poly_List;
    n : constant integer32 := map.n;
    done : array(1..n,1..n) of boolean;
    k : integer32;

  begin
    for i in done'range(1) loop
      for j in i..done'last(2) loop
        done(i,j) := false;
      end loop;
    end loop;
    for i in 1..n loop
      if not Is_Zero(map.c(i)) and not Is_Constant(map,i) then
        for j in (i+1)..n loop
          if not done(i,j) then
            if not Is_Zero(map.c(j)) and not Is_Constant(map,j) then
              k := Share_One_Same_Parameter(map,i,j);
              if k /= 0 then
                declare
                  a : constant integer32 := map.v(j)(k);
                  b : constant integer32 := map.v(i)(k);
                  p : constant Poly
                    := Variable_Relation(map.n,i,j,a,b,map.c(i),map.c(j));
                begin
                  Append(res,res_last,p);
                  done(i,j) := true;
                  for k in (i+1)..(j-1) loop  -- apply transitivity
                    if done(i,k)
                     then done(k,j) := true;
                    end if;
                  end loop;
                end;
              end if;
            end if;
          end if;
        end loop;
      end if;
    end loop;
    return res;
  end Paired_Variables;

  function Number_of_Nonzeroes
             ( v : Standard_integer_Vectors.Vector ) return natural32 is

    res : natural32 := 0;

  begin
    for i in v'range loop
      if v(i) /= 0
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Nonzeroes;

  function Basis_Parameter
             ( map : Monomial_Map; i : integer32 ) return integer32 is

    res : integer32 := 0;
    lv : Standard_Integer_Vectors.Link_to_Vector;

  begin
    if Is_One(map.c(i)) then
      lv := map.v(i);
      for k in lv'range loop
        if lv(k) /= 0 then
          if res = 0
           then res := k;   -- record position of first parameter
           else return 0;   -- we have at least 2 parameters
          end if;
        end if;
      end loop;
    end if;
    return res;
  end Basis_Parameter;

  function Is_In ( v : Standard_Integer_Vectors.Vector;
                   e : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there is a k so that v(k) = e. 

  begin
    for k in v'range loop
      if v(k) = e
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  procedure Exponents_of_Basis
               ( map : in Monomial_Map; i,m : in integer32;
                 exp : out Standard_Integer_Matrices.Matrix;
                 pars,vars : out Standard_Integer_Vectors.Vector ) is
 
  -- DESCRIPTION :
  --   Returns the exponents of the parameters for the i-th variable
  --   in the top column of the m-by-m matrix on return,
  --   with in the other rows the other basis elements.

  -- REQUIRED : m = Number_of_Nonzeroes(map.v(i)) > 1,
  --   exp is an (m+1)-by-m matrix, pars'range = 1..m = vars'range.

  -- ON ENTRY :
  --   map       a monomial map;
  --   i         current variable under consideration;
  --   m         the number of nonzero elements in map.v(i),
  --             which equals the number of parameters used for x(i).

  -- ON RETURN :
  --   exp       first row contains the nonzero elements of map.v(i),
  --             the next rows contains the monomial definitions for
  --             the variables that are basis elements and with their
  --             parameter occurring in map.v(i);
  --   pars      if pars(k) = 1, then the k-th parameter was selected;
  --   vars      if vars(k) = 1, then the k-th variable was selected.

    lvi : constant Standard_Integer_Vectors.Link_to_Vector := map.v(i);
    row,col,k : integer32;

  begin
    pars := (1..m => 0);
    vars := (1..m => 0);
   -- put("lvi = "); put(lvi); new_line;
    row := 1; col := 0;
    for k in 1..map.d loop
      if lvi(k) /= 0 then
        col := col + 1;
        exp(row,col) := lvi(k);
      end if;
    end loop;
    col := 0;
    for j in 1..map.n loop
      if j /= i then
        k := Basis_Parameter(map,j);   -- is x(j) a basis variable ?
       -- put("  j = "); put(j,1);
       -- put("  k = "); put(k,1);
       -- put("  col = "); put(col,1); new_line;
        if k /= 0 then
          if not Is_In(pars(1..col),k) then -- use parameter k only once
            if lvi(k) /= 0 then        -- occurs among parameters of x(i) ?
              col := col + 1;
              pars(col) := k;
              vars(col) := j;
              row := row + 1;
              for kk in 1..m loop
                exp(row,kk) := 0;
              end loop;
              exp(row,col) := map.v(j)(k);
            end if;
          end if;
        end if;
      end if;
      exit when (col = m);
    end loop;
  end Exponents_of_Basis;

  procedure Multiply ( e : in out Standard_Integer_Matrices.Matrix;
                       f : out integer32 ) is

    m,p : integer32;

  begin
    f := 1;
    for k in e'range(2) loop
      if e(1,k) mod e(k+1,k) /= 0 then
        m := Standard_Common_Divisors.LCM(abs(e(1,k)),abs(e(k+1,k)));
        p := m/e(1,k);
        for j in e'range(2) loop
          e(1,j) := p*e(1,j);
        end loop;
        f := f*p;
      end if;
    end loop;
  end Multiply;

  function Basis_Equation ( map : Monomial_Map; i : integer32 ) return Poly is

    res : Poly := Null_Poly;
    lvi : constant Standard_Integer_Vectors.Link_to_Vector := map.v(i);
    dim : constant integer32 := integer32(Number_of_Nonzeroes(lvi.all));

  begin
    if dim > 1 then
     -- put("dimension : "); put(dim,1); new_line;
      declare
        bas : Standard_Integer_Matrices.Matrix(1..dim+1,1..dim);       
        par : Standard_Integer_Vectors.Vector(1..dim);
        var : Standard_Integer_Vectors.Vector(1..dim);
        fac : integer32;
        t1,t2 : Term;
      begin
        Exponents_of_Basis(map,i,dim,bas,par,var);
       -- put_line("the exponent matrix of the basis : "); put(bas);
       -- put("parameters : "); put(par); new_line;
       -- put(" variables : "); put(var); new_line;
        Multiply(bas,fac);
       -- put("after multiplying with "); put(fac,1);
       -- put_line(" :"); put(bas);
        t1.dg := new Standard_Integer_Vectors.Vector'(1..map.n => 0);
        t2.dg := new Standard_Integer_Vectors.Vector'(1..map.n => 0);
        t1.dg(i) := fac;
        t1.cf := Create(1.0);
        for k in 1..dim loop
          if bas(k+1,k) > 0
           then t2.dg(var(k)) := bas(1,k)/bas(k+1,k);
           else t1.dg(var(k)) := bas(1,k)/(-bas(k+1,k));
          end if;
        end loop;
        t2.cf := Create(1.0);
        res := Create(t1);
        Sub(res,t2);
        Clear(t1); Clear(t2);
      end;
    end if;
    return res;
  end Basis_Equation;

  function Equations ( map : Monomial_Map ) return Poly_List is

    res,res_last,pl : Poly_List;

  begin
    for i in 1..map.n loop
      if Is_Zero(map.c(i)) then
        Append(res,res_last,One_Variable_Equation(map.n,i));
      elsif Is_Constant(map,i) then
        Append(res,res_last,One_Variable_Constant_Equation(map.n,i,map.c(i)));
      end if;
    end loop;
    pl := Variable_Differences(map);
    Concatenate(res,res_last,pl); Shallow_Clear(pl);
    pl := Paired_Variables(map);
    Concatenate(res,res_last,pl); Shallow_Clear(pl);
    for i in 1..map.n loop
      declare
        p : constant Poly := Basis_Equation(map,i);
      begin
        if p /= Null_Poly
         then Append(res,res_last,p);
        end if;
      end;
    end loop;
    return res;
  end Equations;

  function Equations ( map : Monomial_Map ) return Laur_Sys is

    eqs : Poly_List := Equations(map);
    res : constant Laur_Sys := Create(eqs);

  begin
    Shallow_Clear(eqs);
    return res;
  end Equations;

end Standard_Monomial_Map_Ideals;
