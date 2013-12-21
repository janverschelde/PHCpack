with text_io;                           use text_io;
with Standard_Random_Numbers;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Affine_Binomials;
with Affine_Binomial_Iterator;

procedure ts_affpol is

-- DESCRIPTION :
--   Looking for affine solution sets of sparse polynomial systems.

  function Reduce_Exponents ( t : Term ) return Term is

  -- DESCRIPTION :
  --   Exponents higher than 1 are reduced to 1.
  --   The coefficients are made at random.
  
    res : Term;

  begin
    res.cf := Standard_Random_Numbers.Random1;
    res.dg := new Standard_Integer_Vectors.Vector(t.dg'range);
    for i in t.dg'range loop
      if t.dg(i) = 0 then
        res.dg(i) := 0;
      elsif t.dg(i) > 0 then
        res.dg(i) := 1;
      else
        res.dg(i) := -1;
      end if;
    end loop;
    return res;
  end Reduce_Exponents;

  function Reduce_Exponents ( p : Poly ) return Poly is

  -- DESCRIPTION :
  --   Returns a polynomials with reduced exponents
  --   and random coefficients.

    res : Poly := Null_Poly;

    procedure Reduce_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Reduce_Exponents(t);

    begin
      Add(res,rt);
      continue := true;
      Clear(rt);
    end Reduce_Term;
    procedure Reduce_Terms is new Visiting_Iterator(Reduce_Term);

  begin
    Reduce_Terms(p);
    return res;
  end Reduce_Exponents;

  function Reduce_Exponents ( p : Laur_Sys ) return Laur_Sys is

  -- DESCRIPTION :
  --   Returns a polynomial system with reduced exponents.

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Reduce_Exponents(p(i));
    end loop;
    return res;
  end Reduce_Exponents;

  function Lengths ( p : Laur_Sys ) return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the number of monomials in every polynomial of p.
  --   The vector on return has range p'range.

    res : Standard_Natural_Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Number_of_Terms(p(i));
    end loop;
    return res;
  end Lengths;

  function Adjacency_Matrix
               ( p : Laur_Sys; m,n : integer32 )
               return Standard_Integer_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a matrix of m rows and n columns with the adjacency
  --   matrix M[x^a,x_k] = 1 if x_k occurs in x^a, or equivalently,
  --   the rows of the matrix on return contain the reduced exponents
  --   of the system.

  -- ON ENTRY :
  --   p        a Laurent system with reduced exponents;
  --   m        total number of monomials in p;
  --   n        number of variables in each polynomial of p.

    res : Standard_Integer_Matrices.Matrix(1..m,1..n);
    ind : integer32 := 0;

    procedure Add ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      for i in t.dg'range loop
        res(ind,i) := t.dg(i);
      end loop;
      continue := true;
    end Add;
    procedure Add_Exponents is new Visiting_Iterator(Add);

  begin
    for i in p'range loop
      Add_Exponents(p(i));
    end loop;
    return res;
  end Adjacency_Matrix;

  generic

    with procedure Report
            ( setzero : in Standard_Integer_Vectors.Vector;
              cntzero : in integer32;
              nonzero : in Standard_Integer_Vectors.Vector;
              cntrest : in integer32;
              continue : out boolean );

    -- DESCRIPTION :
    --   Each time a new choice is found, Report is called.
   
    -- ON ENTRY :
    --   setzero    indices of variables selected to be zero;
    --   cntzero    number of variables selected to be zero;
    --   nonzero    indices of variables that should not be zero,
    --              because of skipped binomial equations.
    --   cntrest    number of remaining equations.
  
    -- ON RETURN :
    --   continue   true if the enumeration should go on,
    --              false if the enumerate must terminate.

  procedure Enumerate ( A : in Standard_Integer_Matrices.Matrix;
                        len : in Standard_Natural_Vectors.Vector;
                        s0_max : in integer32 );

  -- DESCRIPTION :
  --   Enumerates choices of variables to be set to zero,
  --   going over the rows of A.

  -- ON ENTRY :
  --   A        rows of A are indexed by the monomials,
  --            the columns by the variables 
  --   A[i,j] = 1 : the j-th variable has positive power in monomial i,
  --          = 0 : zeroing the j-th variable does nothing to monomial i;
  --   len      indicates the length of each support, 
  --            and the sum of the elements in len equals A'last(1);
  --   s0_max   maximum number of variables which may be set to zero,
  --            to limit the depth of the enumeration tree.

  function Start_Indices
              ( len : Standard_Natural_Vectors.Vector ) 
              return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 1..len'last+1
  --   which contains the start position of support k
  --   provided len[k] equals the length of support k,
  --   for all k in len'range. 

    res : Standard_Natural_Vectors.Vector(len'first..len'last+1);

  begin
    res(res'first) := 1;
    for i in len'range loop
      res(i+1) := res(i) + len(i);
    end loop;
    return res;
  end Start_Indices;

  function Update_Equation_Index
                ( p : Standard_Natural_Vectors.Vector;
                  k,e : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   If p(e+1) = k, then the updated equation index is returned,
  --   otherwise e is returned.

  begin
    if k = integer32(p(e+1))
     then return e+1;
     else return e;
    end if;
  end Update_Equation_Index;
 
  procedure Enumerate ( A : in Standard_Integer_Matrices.Matrix;
                        len : in Standard_Natural_Vectors.Vector;
                        s0_max : in integer32 ) is

    s0 : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    -- s0 contains the selected variables set to be zero
    s0_cnt : integer32 := 0; -- counts number of ones in s0
    s1 : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    -- s1 contains the variables which should not be set to zero
    eq_cnt : integer32 := 0; -- counts number of ones in eq
    pos : constant Standard_Natural_Vectors.Vector := Start_Indices(len);
    continue : boolean := true;

    procedure Enum ( k,e : in integer32 ) is

    -- DESCRIPTION :
    --   Examines the k-th monomial of polynomial e.

    begin
     -- put("k : "); put(k,1); put("  e : "); put(e,1); new_line;
      if k > A'last(1) then
        Report(s0,s0_cnt,s1,eq_cnt,continue);
      elsif Affine_Binomial_Iterator.Set_to_Zero(A,k,s0) then
        Enum(k+1,Update_Equation_Index(pos,k+1,e));
      else
        if k = integer32(pos(e)) then -- can we skip polynomial ?
          put("check if polynomial "); put(e,1);
          put_line(" may be skipped...");
        end if;
        if continue then
          for j in A'range(2) loop
            if A(k,j) > 0 and s1(j) = 0 then
              if s0_cnt < s0_max then
                s0(j) := 1; s0_cnt := s0_cnt + 1;
                Enum(k+1,Update_Equation_Index(pos,k+1,e));
                s0(j) := 0; s0_cnt := s0_cnt - 1;
              end if;
            end if;
            exit when not continue;
          end loop;
        end if;
      end if;
    end Enum;

  begin
    Enum(A'first(1),1);
  end Enumerate;

  procedure Start_Enumeration 
              ( n : in integer32;
                A : in Standard_Integer_Matrices.Matrix;
                len : in Standard_Natural_Vectors.Vector ) is

    cnt : integer32 := 0;
    s0_max : integer32 := n;

    procedure Write ( s0 : in Standard_Integer_Vectors.Vector;
                      s0cnt : in integer32;
                      s1 : in Standard_Integer_Vectors.Vector;
                      eqcnt : in integer32;
                      continue : out boolean ) is
    begin
      cnt := cnt + 1; put(cnt,3); put(" : ");
      put(" s[1] : "); put(s1); 
      put(" #eq : "); put(eqcnt,1);
      put(" s[0] : "); put(s0); 
      put(" # : "); put(s0cnt,1); new_line;
      continue := true;
    end Write;
    procedure Enum is new Enumerate(Report=>Write);
  
  begin
    put("current max is "); put(s0_max,1); new_line;
    put("give new max : "); get(s0_max); 
    Enum(A,len,s0_max);
  end Start_Enumeration;

  procedure Search_for_Affine_Solutions ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Launches an enumeration in search of affine solution sets of p.
  --
    rp : constant Laur_Sys(p'range) := Reduce_Exponents(p);
    ln : constant Standard_Natural_Vectors.Vector := Lengths(rp);
    m : constant integer32 := integer32(Standard_Natural_Vectors.Sum(ln));
    n : constant integer32 := integer32(Number_of_Unknowns(rp(rp'first)));
    A : constant Standard_Integer_Matrices.Matrix(1..m,1..n)
      := Adjacency_Matrix(rp,m,n);

  begin
    put_line("The reduced polynomials : "); put(rp);
    put("The lengths : "); put(ln);
    put(" with sum = "); put(m,1); new_line;
    put_line("The adjacency matrix : "); put(A);
    Start_Enumeration(n,A,ln);
  end Search_for_Affine_Solutions;

  procedure Main is

    lp : Link_to_Laur_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    put_line("The polynomials : "); put(lp.all);
    Search_for_Affine_Solutions(lp.all);
  end Main;

begin
  Main;
end ts_affpol;
