with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;       use Standard_Integer_VecVecs_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Vectors;           use QuadDobl_Random_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;     use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;      use QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;

procedure ts_evalform is

-- DESCRIPTION :
--   This procedure reads a polynomial system and computes the format
--   for evaluation of the system and all its partial derivatives.

--  function Append_Derivatives ( p : Poly_Sys ) return Poly_Sys is
--
--  -- DESCRIPTION :
--  --   Appends to the system all partial derivatives of the Jacobian
--  --   matrix of p, row after row.
--
--  -- REQUIRED : the first polynomial of p is not null
--  --   and has as many variables as all other polynomials.
--
--    nv : constant natural := Number_of_Unknowns(p(p'first));
--    nq : constant natural := p'last;
--    n : constant natural := nq + nq*nv;
--    res : Poly_Sys(1..n);
--    ind : natural := 0;
--
--  begin
--    res(p'range) := p;
--    ind := nq;
--    for i in 1..nq loop
--      for j in 1..nv loop
--        ind := ind+1;
--        res(ind) := Diff(p(i),j);
--      end loop;
--    end loop;
--    return res;
--  end Append_Derivatives;

  procedure Append_Derivatives 
              ( nv,nq : in integer32;
                p : in Poly_Sys; q : out Poly_Sys ) is

  -- DESCRIPTION :
  --   Appends to the system all partial derivatives of the Jacobian
  --   matrix of p, row after row.

  -- ON ENTRY :
  --   nv       the number of variables in p;
  --   nq       the number of equations in p;

  -- ON RETURN :
  --   q        contains p and all its partial derivatives.

  -- REQUIRED : the range of q can hold p and its Jacobian matrix.

    ind : integer32 := 0;

  begin
    q(p'range) := p;
    ind := nq;
    for i in 1..nq loop
      for j in 1..nv loop
        ind := ind+1;
        q(ind) := Diff(p(i),j);
      end loop;
    end loop;
  end Append_Derivatives;

  function Support ( p : Poly ) return List is

    res,res_last : List;

    procedure Visit ( t : in Term; continue : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      Append(res,res_last,v);
      continue := true;
    end Visit;
    procedure Visit_Terms is new Visiting_Iterator(Visit);

  begin
    Visit_Terms(p);
    return res;
  end Support;

  function Distinct_Supports ( s : Array_of_Lists ) return List is

  -- DESCRIPTION :
  --   Returns the list of distinct elements in s.

    res,res_last,tmp : List;

  begin
    for i in s'range loop
      tmp := s(i);
      while not Is_Null(tmp) loop
        declare
          v : constant Standard_Integer_Vectors.Link_to_Vector
            := Head_Of(tmp);  
        begin
          if not Is_In(res,v.all)
           then Append(res,res_last,v.all);
          end if;
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Distinct_Supports;

  function Lex_Greater ( a,b : Standard_Integer_Vectors.Link_to_Vector )
                       return boolean is

  -- DESCRIPTION :
  --   Returns true if a(i) > b(i) and a(j) = b(j) for all j < i.

  -- REQUIRED : a'range = b'range.

  begin
    for i in a'range loop
       if a(i) > b(i) then
         return true;
       elsif a(i) < b(i) then
         return false;
       end if;
    end loop;
    return false; -- a and b are equal
  end Lex_Greater;

  function Sort ( n : integer32; s : List )
                return Standard_Integer_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the sorted list in descending lexicographic order.

  -- REQUIRED : Length_Of(s) = n > 0.

    res : Standard_Integer_VecVecs.VecVec(1..n);
    tmp : List := s;
    v : Standard_Integer_Vectors.Link_to_Vector;
    k : integer32;

  begin
    v := Head_Of(s);
    res(1) := v;
    tmp := Tail_Of(s);
    for i in 2..n loop
      v := Head_Of(tmp);
      k := i;
      for j in 1..(i-1) loop
        if Lex_Greater(v,res(j))
         then k := j;
        end if;
        exit when (k /= i);
      end loop;
      if k < i then
        for j in reverse (k+1)..i loop 
          res(j) := res(j-1);
        end loop;
      end if;
      res(k) := v;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Sort;

  procedure Write_Coefficients ( file : in file_type; p : in Poly ) is

    procedure Visit ( t : in Term; continue : out boolean ) is
    begin
      put(file,t.cf); new_line(file);
      continue := true;
    end Visit;
    procedure Visit_Terms is new Visiting_Iterator(Visit);

  begin
    Visit_Terms(p);
  end Write_Coefficients;

  function Index ( v : Standard_Integer_VecVecs.VecVec; d : Degrees )
                 return integer32 is

  -- DESCRIPTION :
  --   Returns the index k for which v(k) = x, otherwise 0 is returned.

    equal : boolean;

  begin
    for k in v'range loop
      equal := true;
      for i in d'range loop
        equal := (v(k)(i) = integer32(d(i)));
        exit when not equal;
      end loop;
      if equal
       then return k;
      end if;
    end loop;
    return 0;
  end Index;

  function Monomial_Indices
              ( p : Poly; v : Standard_Integer_VecVecs.VecVec )
              return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --    Returns a vector of indices for the monomials in p
  --    the i-th index on return equals the position of the i-th
  --    monomial of p in the vector v.

    nt : constant integer32 := integer32(Number_of_Terms(p));
    res : Standard_Natural_Vectors.Vector(1..nt);
    ind : integer32 := 0;

    procedure Visit ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := natural32(Index(v,t.dg));
      continue := true;
    end Visit;
    procedure Visit_Terms is new Visiting_Iterator(Visit);

  begin
    res := (res'range => 0);
    Visit_Terms(p);
    return res;
  end Monomial_Indices;

  function Monomial_Indices
              ( p : Poly_Sys; v : Standard_Integer_VecVecs.VecVec )
              return Standard_Natural_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns a vector of vectors with the indices of the monomials.

    res : Standard_Natural_VecVecs.VecVec(p'range);

  begin
    for i in p'range loop
      declare
        ind : constant Standard_Natural_Vectors.Vector
            := Monomial_Indices(p(i),v);
      begin
        res(i) := new Standard_Natural_Vectors.Vector'(ind);
      end;
    end loop;
    return res;
  end Monomial_Indices;

  function Eval ( v : in Standard_Integer_VecVecs.VecVec;
                  x : QuadDobl_Complex_Vectors.Vector ) 
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in res'range loop
      res(i) := Create(integer(1));
      for j in v(i)'range loop
        if v(i)(j) > 0 then
          for k in 1..v(i)(j) loop
            res(i) := res(i)*x(j);
          end loop;
        end if;
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( v : in Standard_Integer_VecVecs.VecVec;
                  m : in Standard_Natural_VecVecs.VecVec;
                  c : in QuadDobl_Complex_VecVecs.VecVec;
                  x : QuadDobl_Complex_Vectors.Vector ) 
                return QuadDobl_Complex_Vectors.Vector is

    xv : constant QuadDobl_Complex_Vectors.Vector(v'range) := Eval(v,x);
    res : QuadDobl_Complex_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := Create(integer(0));
      for j in c(i)'range loop
        res(i) := res(i) + c(i)(j)*xv(integer32(m(i)(j)));
      end loop;
    end loop;
    return res;
  end Eval;

  function Coeffs ( p : in Poly ) return QuadDobl_Complex_Vectors.Vector is

    nt : constant integer32 := integer32(Number_of_Terms(p));
    res : QuadDobl_Complex_Vectors.Vector(1..nt);
    ind : integer32 := 0;

    procedure Visit ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := t.cf;
      continue := true;
    end Visit;
    procedure Visit_Terms is new Visiting_Iterator(Visit);

  begin
    Visit_Terms(p);
    return res;
  end Coeffs;

  procedure Evaluate_at_Random_Point
               ( p : in Poly_Sys; nv : in integer32;
                 v : in Standard_Integer_VecVecs.VecVec;
                 m : in Standard_Natural_VecVecs.VecVec ) is

    x : constant QuadDobl_Complex_Vectors.Vector(1..nv) := Random_Vector(1,nv);
    y,z : QuadDobl_Complex_Vectors.Vector(p'range);
    c : QuadDobl_Complex_VecVecs.VecVec(p'range);

  begin
    put_line("evaluating the system once at a random point ...");
    y := Eval(p,x);
    put_line("evaluating the system once more ...");
    for i in p'range loop
      declare
        cff : constant QuadDobl_Complex_Vectors.Vector := Coeffs(p(i));
      begin
        c(i) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    z := Eval(v,m,c,x);
    put_line("comparing values : "); 
    for i in y'range loop
      put(y(i),3); put(" - "); put(z(i),3);
      put(" = "); put(y(i)-z(i),3); new_line;
    end loop;
  end Evaluate_at_Random_Point;

  procedure Write_to_File
               ( file : in file_type; nq : in integer32;
                 s : in Array_of_Lists; p : in Poly_Sys;
                 v : in Standard_Integer_VecVecs.VecVec;
                 m : in Standard_Natural_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Writes the results to file.

  begin
    put(file,nq,1); new_line(file);
    put(file,v'last,1); new_line(file);
    for i in p'range loop
      put(file," "); put(file,Length_Of(s(i)),1);
    end loop;
    new_line(file);
    new_line(file);
    for i in p'range loop
      Write_Coefficients(file,p(i));
      new_line(file);
    end loop;
    put(file,v);
    new_line(file);
    for i in m'range loop
      put(file,m(i).all); new_line(file);
    end loop;
  end Write_to_File;

  procedure Supports ( file : in file_type;
                       nv,nq : in integer32; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Computes the supports of the polynomials in p.

    s : Array_of_Lists(p'range);
    L : List;
    n : integer32;

  begin
    for i in p'range loop
      s(i) := Support(p(i));
      put(" "); put(Length_Of(s(i)),1);
    end loop;
    new_line;
    L := Distinct_Supports(s);
    n := integer32(Length_Of(L));
    put("The number of distinct monomials : "); put(n,1); put_line(".");
    new_line;
    put_line("sorting the monomials in lexicographic descending order...");
    declare
      v : constant Standard_Integer_VecVecs.VecVec(1..n) := Sort(n,L);
      m : Standard_Natural_VecVecs.VecVec(p'range);
    begin
      put_line("the lexicographically sorted list of exponents : "); put(v);
      put_line("computing the positions of the monomials...");
      m := Monomial_Indices(p,v);
      put_line("writing the results to file...");
      Write_to_File(file,nq,s,p,v,m);
      put_line("check: evaluate system at random point...");
      Evaluate_at_Random_Point(p,nv,v,m);
    end;
  end Supports;

  procedure Main is

    p : Link_to_Poly_Sys;
    ans : character;
    file : file_type;

  begin
    new_line;
    put_line("Reading a system with quad double complex coefficients...");
    get(p);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    declare
      nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
      nq : constant integer32 := p'last;
      n : constant integer32 := nq + nq*nv;
      q : Poly_Sys(1..n);
     -- q : constant Poly_Sys := Append_Derivatives(p.all);
    begin
      put("Read a system of "); put(nq,1); put(" equations in ");
      put(nv,1); put_line(" unknowns.");
      put("Do you want to see the polynomials ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then put_line("your system : "); put_line(p.all);
      end if; 
      new_line;
      put_line("computing all partial derivatives ...");
      Append_Derivatives(nv,nq,p.all,q);
      put("Do you want to see the system and all its derivatives ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put_line("your system with all partial derivatives : ");
        put_line(q);
      end if;
      Supports(file,nv,nq,q);
    end;
  end Main;

begin
  Main;
end ts_evalform;
