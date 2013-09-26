with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Multprec_Random_Numbers;           use Multprec_Random_Numbers;
with Multprec_Random_Vectors;           use Multprec_Random_Vectors;
with Standard_Complex_NesVecs;
with Standard_Complex_NesVecs_io;       use Standard_Complex_NesVecs_io;
with Multprec_Complex_NesVecs;
with Multprec_Complex_NesVecs_io;       use Multprec_Complex_NesVecs_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Functions;
with Standard_Univariate_Interpolators;
with Standard_Nvariate_Interpolators;
with Multprec_Univariate_Interpolators;
with Multprec_Nvariate_Interpolators;

procedure ts_newint is

-- DESCRIPTION :
--   Test facility for univariate and n-variate Newton interpolation.

-- UNIVARIATE NEWTON INTERPOLATION :

  procedure Standard_Univariate_Interpolation ( d : in integer32 ) is

  -- DESCRIPTION :
  --   Conducts a random test on the univariate interpolation method:
  --     1) generate a random polynomial and sample it;
  --     2) interpolate at the samples and reconstruct coefficient vector.

    use Standard_Complex_Numbers,Standard_Complex_Vectors;
    use Standard_Univariate_Interpolators;

    c : constant Vector(0..d) := Random_Vector(0,d); -- coefficient vector
    x : constant Vector(0..d) := Random_Vector(0,d); -- interpolation points
    y : Vector(0..d);                        -- function values
    f : Vector(0..d);                        -- divided differences
    e : Vector(0..d);                        -- expanded interpolator
    point : Complex_Number;
    ans : character;

  begin
    put_line("Sampling a random polynomial...");
    for i in x'range loop
      y(i) := Evalc(c,x(i));
    end loop;
    put_line("Computing divided differences...");
    f := Create(x,y);
    put_line("Evaluating at interpolation points :");
    for i in x'range loop
      put("At interpolation point "); put(i,1); put_line(" :");
      put(y(i)); new_line;
      put(Evalf(f,x,x(i))); new_line;
    end loop;
    e := Expand(f,x);
    put_line("Coefficients of generated and computed polynomial : ");
    for i in c'range loop
      put(c(i)); new_line;
      put(e(i)); new_line;
    end loop; 
    put_line("Evaluating at random points :");
    loop
      point := Random1;
      put(Evalc(c,point)); new_line;
      put(Evalf(f,x,point)); new_line;
      put("Do you want to generate another random point ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Univariate_Interpolation;

  procedure Multprec_Univariate_Interpolation
              ( d : in integer32; deci : in natural32 ) is

  -- DESCRIPTION :
  --   Conducts a random test on the univariate interpolation method:
  --     1) generate a random polynomial and sample it;
  --     2) interpolate at the samples and reconstruct coefficient vector.

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;
    use Multprec_Univariate_Interpolators;

    size : constant natural32 := Decimal_to_Size(deci);
    c : Vector(0..d) := Random_Vector(0,d,size);  -- exact coefficient vector
    x : Vector(0..d) := Random_Vector(0,d,size);  -- interpolation points
    y : Vector(0..d);                             -- function values
    f : Vector(0..d);                             -- divided differences
    e : Vector(0..d);                             -- expanded interpolator
    point : Complex_Number;
    ans : character;

  begin
    put_line("Sampling a random polynomial...");
    for i in x'range loop
      y(i) := Evalc(c,x(i));
    end loop;
    put_line("Computing divided differences...");
    f := Create(x,y);
    put_line("Evaluating at interpolation points :");
    for i in x'range loop
      put("At interpolation point "); put(i,1); put_line(" :");
      put(y(i)); new_line;
      put(Evalf(f,x,x(i))); new_line;
    end loop;
    e := Expand(f,x);
    put_line("Coefficients of generated and computed polynomial : ");
    for i in c'range loop
      put(c(i)); new_line;
      put(e(i)); new_line;
    end loop; 
    put_line("Evaluating at random points :");
    loop
      point := Random(size);
      put(Evalc(c,point)); new_line;
      put(Evalf(f,x,point)); new_line;
      Clear(point);
      put("Do you want to generate another random point ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Clear(c); Clear(x); Clear(y); Clear(f); Clear(e);
  end Multprec_Univariate_Interpolation;

-- MULTIVARIATE NEWTON INTERPOLATION -- UTILITIES :

  procedure Initialize ( v : in out Standard_Complex_NesVecs.NesVec ) is

  -- DESCRIPTION :
  --   The multi-dimensional matrix is initialized with zeros.
  --   All submatrices have the same dimension.

    use Standard_Complex_Numbers,Standard_Complex_NesVecs;

  begin
    if v.n = 1 then
      for i in v.v'range loop
        v.v(i) := Create(0.0);
      end loop;
    else
      for i in v.w'range loop
        v.w(i) := new NesVec(v.n-1,v.a,v.b);
        Initialize(v.w(i).all);
      end loop;
    end if;
  end Initialize;

  procedure Initialize ( v : in out Multprec_Complex_NesVecs.NesVec ) is

  -- DESCRIPTION :
  --   The multi-dimensional matrix is initialized with zeros.
  --   All submatrices have the same dimension.

    use Multprec_Complex_Numbers,Multprec_Complex_NesVecs;

  begin
    if v.n = 1 then
      for i in v.v'range loop
        v.v(i) := Create(integer(0));
      end loop;
    else
      for i in v.w'range loop
        v.w(i) := new NesVec(v.n-1,v.a,v.b);
        Initialize(v.w(i).all);
      end loop;
    end if;
  end Initialize;

  function Sample ( p : Standard_Complex_Poly_Functions.Eval_Poly;
                    x : Standard_Complex_VecVecs.VecVec )
                  return Standard_Complex_NesVecs.NesVec is

  -- DESCRIPTION :
  --   Returns the values of the polynomial p at the grid x,
  --   defined by the products of all vectors in x.

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors,Standard_Complex_NesVecs;
    use Standard_Complex_Poly_Functions;

    fst : constant Link_to_Vector := x(x'first);
    res : NesVec(natural32(x'last),fst'first,fst'last);
    acc : Vector(x'range);

    procedure Recursive_Sample
                ( k : in integer32; accu : in out Vector;
                  y : in out NesVec ) is

    -- DESCRIPTION :
    --   Accumulates grid points from level k to the next level,
    --   and evaluates when the accumulator is full.

    begin
      for j in x(k)'range loop
        accu(k) := x(k)(j);
        if k < x'last
         then Recursive_Sample(k+1,accu,y.w(j).all);
         else y.v(j) := Eval(p,accu);
        end if;
      end loop;
    end Recursive_Sample;

  begin
    Initialize(res);
   -- put_line("The initialized multi-dimensional matrix : "); put(res);
    Recursive_Sample(x'first,acc,res);
   -- put_line("The multi-dimensional matrix of samples : "); put(res);
    return res;
  end Sample;

  function Sample ( p : Multprec_Complex_Poly_Functions.Eval_Poly;
                    x : Multprec_Complex_VecVecs.VecVec )
                  return Multprec_Complex_NesVecs.NesVec is

  -- DESCRIPTION :
  --   Returns the values of the polynomial p at the grid x,
  --   defined by the products of all vectors in x.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors,Multprec_Complex_NesVecs;
    use Multprec_Complex_Poly_Functions;

    fst : constant Link_to_Vector := x(x'first);
    res : NesVec(natural32(x'last),fst'first,fst'last);
    acc : Vector(x'range);

    procedure Recursive_Sample
                ( k : in integer32; accu : in out Vector;
                  y : in out NesVec ) is

    -- DESCRIPTION :
    --   Accumulates grid points from level k to the next level,
    --   and evaluates when the accumulator is full.

    begin
      for j in x(k)'range loop
        accu(k) := x(k)(j);
        if k < x'last
         then Recursive_Sample(k+1,accu,y.w(j).all);
         else y.v(j) := Eval(p,accu);
        end if;
      end loop;
    end Recursive_Sample;

  begin
    Initialize(res);
    put_line("The initialized multi-dimensional matrix : "); put(res);
    Recursive_Sample(x'first,acc,res);
    put_line("The multi-dimensional matrix of samples : "); put(res);
    return res;
  end Sample;

  function Retrieve ( v : in Standard_Complex_NesVecs.NesVec;
                      i : in Standard_Integer_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the number in v corresponding to the given index vector.

  -- REQUIRED : dimensions of v and i must match.

  begin
    if v.n = 1
     then return v.v(i(i'first));
     else return Retrieve(v.w(i(i'first)).all,i(i'first+1..i'last));
    end if;
  end Retrieve;

  procedure Query_Samples ( p : in Standard_Complex_Poly_Functions.Eval_Poly;
                            x : in Standard_Complex_VecVecs.VecVec;
                            y : in Standard_Complex_NesVecs.NesVec ) is

  -- DESCRIPTION :
  --   This is an interactive testing facility to compare the samples
  --   stored in the matrix y with the actual function values of p at x.

    use Standard_Complex_Poly_Functions;

    index : Standard_Integer_Vectors.Vector(x'range);
    point : Standard_Complex_Vectors.Vector(x'range);
    ans : character;

  begin
    loop
      put("Reading "); put(x'last,1); put_line(" indices : ");
      for i in index'range loop
        loop
          put("  index "); put(i,1); put(" : ");
          index(i) := 0; get(index(i));
          if (index(i) >= x(i)'first) and (index(i) <= x(i)'last) then
            exit;
          else
            put("Index must be between ");
            put(x(i)'first,1); put(" and ");
            put(x(i)'last,1); put_line(".  Please try again");
          end if;
        end loop;
      end loop;
      put("The indices read : "); put(index); new_line;
      for i in index'range loop
        point(i) := x(i)(index(i));
      end loop;
      put("Value retrieved : "); put(Retrieve(y,index)); new_line;
      put("Value computed  : "); put(Eval(p,point)); new_line;
      put("Do you want more retrievals ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Query_Samples;

  procedure Write_Triangular_Grid
              ( v : in Standard_Complex_NesVecs.NesVec ) is

  -- DESCRIPTION :
  --   Writes the triangular portion of the grid, used for interpolation.

    use Standard_Complex_NesVecs;

    procedure Rec_Write ( nv : in NesVec; nb : in integer32 ) is
    begin
      put(nv.n,1);
      put(" "); put(nv.a,1);
      put(" "); put(nb,1); new_line;
      if nv.n = 1 then
        for i in nv.v'first..nb loop
          put(nv.v(i)); new_line;
        end loop;
      else
        for i in nv.w'first..nb loop
          Rec_Write(nv.w(i).all,nb-i);
        end loop;
      end if;
    end Rec_Write;

  begin
    Rec_Write(v,v.b);
  end Write_Triangular_Grid;

  procedure Write_Triangular_Grid
              ( v : in Multprec_Complex_NesVecs.NesVec ) is

  -- DESCRIPTION :
  --   Writes the triangular portion of the grid, used for interpolation.

    use Multprec_Complex_NesVecs;

    procedure Rec_Write ( nv : in NesVec; nb : in integer32 ) is
    begin
      put(nv.n,1);
      put(" "); put(nv.a,1);
      put(" "); put(nb,1); new_line;
      if nv.n = 1 then
        for i in nv.v'first..nb loop
          put(nv.v(i)); new_line;
        end loop;
      else
        for i in nv.w'first..nb loop
          Rec_Write(nv.w(i).all,nb-i);
        end loop;
      end if;
    end Rec_Write;

  begin
    Rec_Write(v,v.b);
  end Write_Triangular_Grid;

  procedure Standard_Test_Evaluation 
               ( n,d : in integer32;
                 f : in Standard_Complex_NesVecs.NesVec;
                 x : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   This is an interactive testing facility to evaluate
  --   the Newton form of the interpolating polynomial.

    use Standard_Complex_Numbers,Standard_Complex_Vectors;
    use Standard_Nvariate_Interpolators;

    pt : Vector(1..n) := (1..n => Create(1.0));
    eva0,eva1,diff : Complex_Number;
    absdiff : double_float;
    ans : character;

  begin
    put_line("The formula for the plain evaluation : ");
    eva0 := Eval0(Standard_Output,f,x,pt);
    put("Give a character to continue : "); get(ans);
    put_line("The formula for evaluation with Horner's scheme : ");
    eva1 := Eval(Standard_Output,f,x,pt);
    put("Give a character to continue : "); get(ans);
    put_line("Computing the sum of the coefficients : ");
    put("  plain evaluation   : "); put(eva0); new_line;
    put("  with Horner scheme : "); put(eva1); new_line;
    diff := eva0 - eva1;
    absdiff := AbsVal(diff);
    put("Absolute Value of the Difference : "); put(absdiff); new_line;
    put_line("Evaluation at random points :");
    loop
      pt := Random_Vector(1,n);
      eva0 := Eval0(f,x,pt);
      eva1 := Eval(f,x,pt);
      put("  plain evaluation   : "); put(eva0); new_line;
      put("  with Horner scheme : "); put(eva1); new_line;    
      diff := eva0 - eva1;
      absdiff := AbsVal(diff);
      put("Absolute Value of the Difference : "); put(absdiff); new_line;
      put("Do you want another random point ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Test_Evaluation;

-- MULTIVARIATE NEWTON INTERPOLATION -- MAIN TESTING PROCEDURES :

  procedure Standard_Multivariate_Interpolation ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   This procedure tests the interpolation of a polynomial function
  --   in n variables of degree d, over the standard complex numbers.
  --   The user given polynomial is reconstructed by interpolation.

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors,Standard_Complex_VecVecs;
    use Standard_Complex_NesVecs,Standard_Nvariate_Interpolators;
    use Standard_Complex_Polynomials,Standard_Complex_Poly_Functions;

    x : VecVec(1..n);
    y,f : NesVec(natural32(n),0,d);
    p,ip : Poly;
    ep : Eval_Poly;
    pt : Vector(1..n) := (1..n => Create(1.0));
    maxerr : double_float;
    ans : character;

  begin
    Symbol_Table.Init(natural32(n));
    put("Give a polynomial in "); put(n,1); 
    put(" variables of degree "); put(d,1); put_line(" :");
    get(p);
    put("Your polynomial : "); put(p); new_line;
    for i in 1..n loop
      x(i) := new Vector'(Random_Vector(0,d));
    end loop;
    ep := Create(p);
    y := Sample(ep,x);
    put_line("The triangular grid : ");
    Write_Triangular_Grid(y);
    Query_Samples(ep,x,y);
    put("Compute extra higher order divided differences ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      f := Create_on_Square(natural32(n),natural32(d),x,y);
      put_line("The generalized divided differences : "); put(f);
    else
      f := Create_on_Square(natural32(n),natural32(d),x,y);
      put_line("The generalized divided differences : ");
      Write_Triangular_Grid(f);
    end if;
    Standard_Test_Evaluation(n,d,f,x);
    Eval_Grid(Standard_Output,x,y,f,maxerr);
    put("Maximal error on grid recomputed : ");
    put(Maximal_Error(x,y,f)); new_line;
    ip := Expand(f,x);
    put_line("The Newton form expanded as multivariate polynomial : ");
    put_line(ip);
  end Standard_Multivariate_Interpolation;

  procedure Multprec_Multivariate_Interpolation
              ( n,d : in integer32; deci : in natural32 ) is

  -- DESCRIPTION :
  --   This procedure tests the interpolation of a polynomial function
  --   in n variables of degree d, over multi-precision complex numbers.
  --   The user given polynomial which is reconstructed by interpolation.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors,Multprec_Complex_VecVecs;
    use Multprec_Complex_NesVecs,Multprec_Nvariate_Interpolators;
    use Multprec_Complex_Polynomials,Multprec_Complex_Poly_Functions;

    size : constant natural32 := Decimal_to_Size(deci);
    x : VecVec(1..n);
    y,f : NesVec(natural32(n),0,d);
    p,ip : Poly;
    ep : Eval_Poly;
    pt : Vector(1..n) := (1..n => Create(integer(1)));
    eva : Complex_Number;
    maxerr : Floating_Number;
    ans : character;

  begin
    Symbol_Table.Init(natural32(n));
    put("Give a polynomial in "); put(n,1);
    put(" variables of degree "); put(d,1); put_line(" :");
    get(p);
    put("Your polynomial : "); put(p); new_line;
    for i in 1..n loop
      x(i) := new Vector'(Random_Vector(0,d,size));
    end loop;
    ep := Create(p);
    y := Sample(ep,x);
    put("Compute extra higher order divided differences ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      f := Create_on_Square(natural32(n),natural32(d),x,y);
      put_line("The generalized divided differences : "); put(f);
    else
      f := Create(natural32(n),natural32(d),x,y);
      put_line("The generalized divided differences : "); 
      Write_Triangular_Grid(f);
    end if;
    eva := Eval0(Standard_Output,f,x,pt);
    put("Give a character to continue : "); get(ans);
  --  put_line("The formula for evaluation with Horner's scheme : ");
  --  eva := Eval(Standard_Output,f,x,pt);
  --  put("Give a character to continue : "); get(ans);
  --  put_line("Computing the sum of the coefficients : ");
  --  put("  with output    : "); put(eva); new_line;
  --  eva := Eval(f,x,pt);
  --  put("  without output : "); put(eva); new_line;
    Eval_Grid(Standard_Output,x,y,f,maxerr);
    put("Maximal error on grid recomputed : ");
    put(Maximal_Error(x,y,f)); new_line;
    ip := Expand(f,x);
    put_line("The Newton form expanded as multivariate polynomial : ");
    put_line(ip);
  end Multprec_Multivariate_Interpolation;

  procedure Main is

    n,d : integer32 := 0;
    deci : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line
      ("Univariate and n-Variate Interpolation with Divided Differences.");
    new_line;
    loop
      put("Give the number of variables : "); get(n);
      put("Give the degree : "); get(d);
      put("Give number of decimal places (<= 16 : standard) : "); get(deci);
      if deci <= 16 then
        if n = 1
         then Standard_Univariate_Interpolation(d);
         else Standard_Multivariate_Interpolation(n,d);
        end if;
      else
        if n = 1 
         then Multprec_Univariate_Interpolation(d,deci);
         else Multprec_Multivariate_Interpolation(n,d,deci);
        end if;
      end if;
      put("Do you want to test for another polynomial ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Main;

begin
  Main;
end ts_newint;
