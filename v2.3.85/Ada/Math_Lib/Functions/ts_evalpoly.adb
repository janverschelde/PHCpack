with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Symbol_Table;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors; 
with Standard_Natural_Matrices;         use Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;      use Standard_Natural_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Random_Polynomials;       use Standard_Random_Polynomials;

procedure ts_evalpoly is

-- DESCRIPTION :
--   Development code of a nested Horner scheme for multivariate polynomials
--   with complex coefficients.  The algorithm first lexicographically sorts
--   the tableau format of the polynomial and then proceeds recursively.

  procedure Tab ( p : in Poly;
                  c : out Standard_Complex_Vectors.Vector;
                  d : out Standard_Natural_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Returns the tableau representation of the polynomial p
  --   as a coefficient vector c and degree matrix d.

  -- REQUIRED : c'range = d'range(1) = 1..Number_of_Terms(p)
  --   and d'range(2) = 1..Number_of_Unknowns(p).

  -- ON INPUT :
  --   p        a polynomial in n unknowns with m terms.

  -- ON RETURN :
  --   c        the vector with the m coefficients of p;
  --   d        the degree matrix of the polynomials,
  --            with in its rows the exponents of the monomials.

   -- m : constant natural := d'last(1);
   -- n : constant natural := d'last(2);
    ind : integer32 := 0;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      c(ind) := t.cf;
      for k in t.dg'range loop
        d(ind,k) := t.dg(k);
      end loop;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
  end Tab;

  function Lexicographically_Smaller
              ( d : Standard_Natural_Matrices.Matrix; i,j : integer32 )
              return boolean is

  -- DESCRIPTION :
  --   Returns true if row i of d is lexicographically smaller than row j.
  --   Returns false if otherwise (which includes equal rows).

  -- REQUIRED : i and j are both in d'range(1).

  begin
    for k in d'range(2) loop
      if d(i,k) > d(j,k) then
        return false;
      elsif d(i,k) < d(j,k) then
        return true;
      end if;
    end loop;
    return false;
  end Lexicographically_Smaller;

  procedure Swap ( c : in out Standard_Complex_Vectors.Vector;
                   d : in out Standard_Natural_Matrices.Matrix;
                   i,j : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps entries i and j in c and rows i and j of d.

  -- REQUIRED : i and j are both in d'range(1).

    ct : constant Complex_Number := c(i);
    dt : natural32;

  begin
    c(i) := c(j); c(j) := ct; 
    for k in d'range(2) loop
      dt := d(i,k); d(i,k) := d(j,k); d(j,k) := dt;
    end loop;
  end Swap;

  procedure Lex ( c : in out Standard_Complex_Vectors.Vector;
                  d : in out Standard_Natural_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Sorts the coefficients and exponents in lexicographic order.

  -- REQUIRED : c'range = d'range(1).

    max : integer32;

  begin
    for i in d'range(1) loop
      max := i;                -- assume i is lexicographically largest
      for j in i+1..d'last(1) loop
        if Lexicographically_Smaller(d,max,j)
         then max := j;        -- found new maximum
        end if;
      end loop;
      if max /= i              -- swap if necessary 
       then Swap(c,d,i,max);
      end if;
    end loop;
  end Lex;

  function Evaluate
             ( c : Standard_Complex_Vectors.Vector;
               d : Standard_Natural_Matrices.Matrix;
               x : Standard_Complex_Vectors.Vector )
             return Complex_Number is

  -- DESCRIPTION :
  --   Evaluates the polynomial in a straightforward manner
  --   via the inner product of the coefficient vector and
  --   the evaluated monomials.

    res : Complex_Number := Create(0.0);
    trm : Complex_Number;

  begin
    for i in d'range(1) loop
      trm := c(i);
      for j in d'range(2) loop
        for k in 1..d(i,j) loop
          trm := trm*x(j);
        end loop;
      end loop;
      res := res + trm;
    end loop;
    return res;
  end Evaluate;

  function Number_of_Monomials ( n,d : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of monomials of a dense polynomial 
  --   in n variables of degree d.

  begin
    if d = 0 or n = 0 then
      return 1;
    elsif n = 1 then
      return d + 1;
    else
      return Number_of_Monomials(n,d-1) + Number_of_Monomials(n-1,d);
    end if;
  end Number_of_Monomials;

  function Evaluate_Dense_Polynomial
             ( c : Standard_Complex_Vectors.Vector;
               d : Standard_Natural_Matrices.Matrix;
               x : Standard_Complex_Vectors.Vector )
             return Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the polynomial represented in tableau format.
  --   Assumes the tableau comes from a general dense polynomial.

  -- REQUIRED : d is lexicographically sorted in descending order.
  --   x'range = d'range(2), c'range = d'range(1).

    res : Complex_Number;
    n : constant integer32 := d'last(2);
    m : constant integer32 := d'last(1);
    wrk : Standard_Complex_Vectors.Vector(c'range) := c;
    ind,ptr : integer32;

  begin
    if n = 1 then                       -- plain Horner
      res := c(1);
      for i in 2..m loop
        res := res*x(1) + c(i);
      end loop;
    else
      ind := n;                         -- use x(n) to evaluate
      while ind < m loop
        for j in ind+1..ind+integer32(d(ind,n)) loop
          wrk(ind) := wrk(ind)*x(n) + c(j);
        end loop;
        ind := ind + integer32(d(ind,n)) + 1;
        while ind < m and d(ind,n) = 0 loop
          ind := ind + 1;
        end loop;
      end loop;
      for k in reverse 2..(n-1) loop    -- use x(k), 1 < k < n
        ind := k;
        while ind < m loop
          if d(ind-1,k) < d(ind,k) then
            ptr := ind;
            while ptr < m loop
              if d(ptr,k) > d(ptr+1,k)
               then wrk(ind) := wrk(ind)*x(k) + wrk(ptr+1); 
              end if;
              ptr := ptr + 1;
              exit when d(ptr,k) = 0;
            end loop;
            ind := ptr;
          end if;
          ind := ind + 1;
          while ind < m and d(ind,k) = 0 loop
            ind := ind + 1;
          end loop;
        end loop;
      end loop;
      res := wrk(1);                    -- use x(1) to evaluate
      ind := 1;
      for i in 1..d(1,1) loop
        res := res*x(1) + wrk(ind+1);
        ind := ind + integer32(Number_of_Monomials(natural32(n)-1,i));
      end loop;
    end if;
    return res;
  end Evaluate_Dense_Polynomial;

  function Evaluate_Any_Polynomial
             ( c : Standard_Complex_Vectors.Vector;
               d : Standard_Natural_Matrices.Matrix;
               x : Standard_Complex_Vectors.Vector )
             return Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the polynomial represented in tableau format.
  --   Makes no assumptions on density of the polynomial.

  -- REQUIRED : d is lexicographically sorted in descending order.
  --   x'range = d'range(2), c'range = d'range(1).

    res : Complex_Number;
    n : constant integer32 := d'last(2);
    m : constant integer32 := d'last(1);
    wrk : Standard_Complex_Vectors.Vector(c'range) := c;
    gap,ind,ptr : integer32;
    stop,cut : boolean;

  begin
    if n = 1 then
      res := c(1);
      for i in 2..m loop
        gap := integer32(d(i-1,1)) - integer32(d(i,1));
        for j in 1..gap loop
          res := res*x(1);
        end loop;
        res := res + c(i);
      end loop;
      for j in 1..d(m,1) loop
        res := res*x(1);
      end loop;
    else
      ind := 1;                      -- use x(n) to evaluate
      if d(1,n) > 0 then             -- use x(n) at special start
        while ind < m loop
          stop := (d(ind,n) <= d(ind+1,n));
          exit when stop;
          for j in reverse 1..(n-1) loop
            stop := (d(ind,j) /= d(ind+1,j));
            exit when stop;
          end loop;
          exit when stop;
          gap := integer32(d(ind,n)) - integer32(d(ind+1,n));
          for j in 1..gap loop
            wrk(1) := wrk(1)*x(n);
          end loop;
          ind := ind + 1;
          wrk(1) := wrk(1) + c(ind);
        end loop;
        for j in 1..d(ind,n) loop
          wrk(1) := wrk(1)*x(n);
        end loop;
      end if;
      ptr := ind + 1;
      while ptr <= m loop               -- use x(n) in general
        while ptr <= m and then d(ptr,n) = 0 loop
          ptr := ptr + 1;
        end loop;
        exit when (ptr > m);
        cut := false;
        for j in reverse 1..n-1 loop
          cut := (d(ptr-1,j) /= d(ptr,j));
          exit when cut;
        end loop;
        if cut then
          ind := ptr;
          while ind < m loop
            stop := (d(ind,n) <= d(ind+1,n));
            exit when stop;
            for j in reverse 1..(n-1) loop
              stop := (d(ind,j) /= d(ind+1,j));
              exit when stop;
            end loop;
            exit when stop;
            gap := integer32(d(ind,n)) - integer32(d(ind+1,n));
            for j in 1..gap loop
              wrk(ptr) := wrk(ptr)*x(n);
            end loop;
            ind := ind + 1;
            wrk(ptr) := wrk(ptr) + c(ind);
          end loop;
          for j in 1..d(ind,n) loop
            wrk(ptr) := wrk(ptr)*x(n);
          end loop;
        end if;
        ptr := ptr + 1;
        if ptr < ind
         then ptr := ind;
        end if;
      end loop;
      for k in reverse 2..(n-1) loop  -- use x(k), 1 < k < n
        ind := 1;
        if d(1,k) > 0 then
          stop := (ind = m);
          while ind < m loop
            stop := (d(ind,k) < d(ind+1,k));
            exit when stop;
            for j in reverse 1..(k-1) loop
              stop := (d(ind,j) /= d(ind+1,j));
              exit when stop;
            end loop;
            exit when stop;
            gap := integer32(d(ind,k)) - integer32(d(ind+1,k));
            if gap > 0 then
              for j in 1..gap loop
                wrk(1) := wrk(1)*x(k);
              end loop;
              wrk(1) := wrk(1) + wrk(ind+1);
            end if;
            ind := ind + 1;
            exit when d(ind,k) = 0;
          end loop;
          if stop or ind = m then
            for j in 1..d(ind,k) loop
              wrk(1) := wrk(1)*x(k);
            end loop;
          end if;
        end if;
        ptr := ind + 1;
        while ptr <= m loop
          while ptr <= m and then d(ptr,k) = 0 loop
            ptr := ptr + 1;
          end loop;
          exit when (ptr > m);
          cut := false;
          for j in reverse 1..(k-1) loop
            cut := (d(ptr-1,j) /= d(ptr,j));
            exit when cut;
          end loop;
          if cut then
            ind := ptr;
            stop := (ind = m);
            while ind < m loop
              stop := (d(ind,k) < d(ind+1,k));
              exit when stop;
              for j in reverse 1..(k-1) loop
                stop := (d(ind,j) /= d(ind+1,j));
                exit when stop;
              end loop;
              exit when stop;
              gap := integer32(d(ind,k)) - integer32(d(ind+1,k));
              if gap > 0 then
                for j in 1..gap loop
                  wrk(ptr) := wrk(ptr)*x(k);
                end loop;
                wrk(ptr) := wrk(ptr) + wrk(ind+1);
              end if;
              ind := ind + 1;
              exit when d(ind,k) = 0;
            end loop;
            if stop or ind = m then
              for j in 1..d(ind,k) loop
                wrk(ptr) := wrk(ptr)*x(k);
              end loop;
            end if;
          end if;
          ptr := ptr + 1;
          if ptr < ind
           then ptr := ind;
          end if;
        end loop;
      end loop;
      res := wrk(1);                  -- use x(1) to evaluate
      ind := 1;
      while ind < m loop
        gap := integer32(d(ind,1)) - integer32(d(ind+1,1));
        if gap > 0 then
          for j in 1..gap loop
            res := res*x(1);
          end loop;
          res := res + wrk(ind+1);
        end if;
        ind := ind + 1;
      end loop;
      for j in 1..d(m,1) loop
        res := res*x(1);
      end loop;
    end if;
    return res;
  end Evaluate_Any_Polynomial;

 -- function Evaluate_Forward
 --            ( c : Standard_Complex_Vectors.Vector;
 --              d : Standard_Natural_Matrices.Matrix;
 --              x : Standard_Complex_Vectors.Vector )
 --            return Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the polynomial stored in the lexicographically
  --   sorted tableau format.

 --   res : constant Complex_Number := Create(0.0);
 --   m : constant integer32 := c'last;
 --   first,next : Standard_Natural_Vectors.Vector(x'range);
 --   ind,current : integer32;

 -- begin
 --   first := (x'range => 1);
 --   current := x'last;
 --   for i in next'range loop
 --     ind := integer32(first(i)) + 1 ;
 --     while ind < m and d(ind,i) = d(integer32(first(i)),i) loop
 --       ind := ind + 1;
 --     end loop;
 --     next(i) := natural32(ind);
 --     if next(i) = first(i) + 1
 --      then current := i;
 --     end if;
 --     exit when (current < x'last);
 --   end loop;
 --   return res;
 -- end Evaluate_Forward;

  function Maximal_Degrees
             ( d : Standard_Natural_Matrices.Matrix )
             return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   The vector on return has range equal to d'range(2)
  --   and its k-th entry contains the maximal power of x(k) in d.

    res : Standard_Natural_Vectors.Vector(d'range(2)) := (d'range(2) => 0);

  begin
    for i in d'range(1) loop
      for k in d'range(2) loop
        if d(i,k) > res(k)
         then res(k) := d(i,k);
        end if;
      end loop;
    end loop;
    return res;
  end Maximal_Degrees;

  function Maximum ( d : Standard_Natural_Vectors.Vector ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximum of the elements in d.

    res : natural32 := d(d'first);

  begin
    for i in d'first+1..d'last loop
      if d(i) > res
       then res := d(i);
      end if;
    end loop;
    return res;
  end Maximum;

  procedure Build_Power_Matrix
              ( m : out Standard_Complex_Matrices.Matrix;
                x : in Standard_Complex_Vectors.Vector;
                p : in Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns in m(i,j) the powers of x(i) for j in 1..p(i).

  -- REQUIRED : m'range(1) = x'range = p'range
  --   and m'range(2) = 1..Maximum(p).

  begin
    for i in x'range loop
      m(i,1) := x(i);
      for j in 2..integer32(p(i)) loop
        m(i,j) := m(i,j-1)*x(i);
      end loop;
    end loop;
  end Build_Power_Matrix;

  function Evaluate_Sparse_Polynomial
             ( c : Standard_Complex_Vectors.Vector;
               d : Standard_Natural_Matrices.Matrix;
               x : Standard_Complex_Vectors.Vector )
             return Complex_Number is

  -- DESCRIPTION :
  --   Evaluates a sparse polynomial stored in tableau format,
  --   using cache for the powers.

    md : constant Standard_Natural_Vectors.Vector(x'range)
       := Maximal_Degrees(d);
    mp : constant natural32 := Maximum(md);
    xp : Standard_Complex_Matrices.Matrix(x'range,1..integer32(mp));
    tx : Complex_Number;
    res : Complex_Number := Create(0.0);

  begin
    Build_Power_Matrix(xp,x,md);
    for i in c'range loop
      tx := c(i);
      for k in x'range loop
        if d(i,k) > 0 
         then tx := tx*xp(k,integer32(d(i,k)));
        end if;
      end loop;
      res := res + tx;
    end loop;
    return res;
  end Evaluate_Sparse_Polynomial;

  procedure Horner ( p : in Poly; dense : in boolean ) is

  -- DESCRIPTION :
  --   Sorts the tableau representation of p for the Horner scheme.
  --   The boolean variable dense is true if p is a dense polynomial.

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    m : constant integer32 := integer32(Number_of_Terms(p));
    c : Standard_Complex_Vectors.Vector(1..m);
    d : Standard_Natural_Matrices.Matrix(1..m,1..n);
    x : Standard_Complex_Vectors.Vector(1..n);
    y,z,e,s : Complex_Number;
    ans : character;

  begin
    Tab(p,c,d);
    put_line("The coefficient vector : "); put_line(c);
    put_line("The degree matrix : "); put(d);
    Lex(c,d);
    put_line("The lexicographically sorted degree matrix : "); put(d);
    put_line("The corresponding coefficient vector : "); put_line(c);
    loop
      x := Random_Vector(1,n);
      y := Eval(p,x);
      put(" Evaluated at random : "); put(y); new_line;
      e := Evaluate(c,d,x);
      put("    straight tableau : "); put(e); new_line;
      if dense then
        z := Evaluate_Dense_Polynomial(c,d,x);
        put("  value with tableau : "); put(z); new_line;
      else
        z := Evaluate_Any_Polynomial(c,d,x);
        put("  value with tableau : "); put(z); new_line;
        s := Evaluate_Sparse_Polynomial(c,d,x);
        put("sparse tableau value : "); put(s); new_line;
      end if;
      put("More random points ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Horner;

  procedure Generate_Polynomial ( p : out Poly; dense : out boolean ) is

  -- DESCRIPTION :
  --   Interactive generation of a polynomial p.  If dense is true
  --   on return, then a general dense polynomial was generated.

    n,d,m : natural32 := 0;

  begin
    new_line;
    put("Give number of variables : "); get(n);
    put("Give the highest degree  : "); get(d);
    put("Give number of monomials (0 for dense) : "); get(m);
    dense := (m = 0);
    Symbol_Table.Init(n);
    if dense then
      p := Random_Dense_Poly(n,d,0);
      put("A random dense polynomial of degree ");
    else
      p := Random_Sparse_Poly(n,d,m,0); 
      put("A random polynomial with at most "); put(m,1);
      put_line(" monomials,"); put(" of degree ");
    end if;
    put(d,1); put(" in "); put(n,1); put_line(" variables :");
    put_line(p);
  end Generate_Polynomial;

  procedure Horner_Test is

    p : Poly;
    dense : boolean;

  begin
    Generate_Polynomial(p,dense);
    Horner(p,dense);
  end Horner_Test;

  procedure Performance_Test is

    p : Poly;
    dense : boolean;
    nb : natural32 := 0;
    timer : timing_widget;

  begin
    Generate_Polynomial(p,dense);
    put("Give number of evaluations : "); get(nb);
    declare
      ep : constant Eval_Poly := Create(p);
      n : constant integer32 := integer32(Number_of_Unknowns(p));
      m : constant integer32 := integer32(Number_of_Terms(p));
      c : Standard_Complex_Vectors.Vector(1..m);
      d : Standard_Natural_Matrices.Matrix(1..m,1..n);
      x : Standard_Complex_Vectors.Vector(1..n);
      y,z : Complex_Number;
    begin
      Tab(p,c,d);
      Lex(c,d);
      x := Random_Vector(1,n);
      tstart(timer);
      for i in 1..nb loop
        y := Eval(p,x);
      end loop;
      tstop(timer);
      put("user cpu time for list evaluation : ");
      print_hms(standard_output,Elapsed_User_Time(timer)); new_line;
      tstart(timer);
      for i in 1..nb loop
        y := Eval(ep,x);
      end loop;
      tstop(timer);
      put("cpu time for eval_poly evaluation : ");
      print_hms(standard_output,Elapsed_User_Time(timer)); new_line;
      tstart(timer);
      for i in 1..nb loop
        z := Evaluate(c,d,x);
      end loop;
      tstop(timer);
      put("  for straight tableau evaluation : ");
      print_hms(standard_output,Elapsed_User_Time(timer)); new_line;
      if dense then
        tstart(timer);
        for i in 1..nb loop
          z := Evaluate_Dense_Polynomial(c,d,x);
        end loop;
        tstop(timer);
      else
        tstart(timer);
        for i in 1..nb loop
          z := Evaluate_Any_Polynomial(c,d,x);
        end loop;
        tstop(timer);
      end if;
      put("  cpu time for tableau evaluation : ");
      print_hms(standard_output,Elapsed_User_Time(timer)); new_line;
      tstart(timer);
      for i in 1..nb loop
        z := Evaluate_Sparse_Polynomial(c,d,x);
      end loop;
      tstop(timer);
      put("    for sparse tableau evaluation : ");
      print_hms(standard_output,Elapsed_User_Time(timer)); new_line;
    end;
  end Performance_Test;

  procedure Main is

    ans : character;
 
  begin
    new_line;
    put_line("MENU to test evaluation of multivariate polynomials :");
    put_line("  1. test new and compare with original evaluation;");
    put_line("  2. performance test between old and new methods.");
    put("Type 1 or 2 to select : "); Ask_Alternative(ans,"12");
    if ans = '1'
     then Horner_Test;
     else Performance_Test;
    end if;
  end Main;

begin
  Main;
end ts_evalpoly;
