with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
--with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
--with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
--with Standard_Natural64_Vectors;
with Standard_Natural64_Vectors_io;     use Standard_Natural64_Vectors_io;
with Standard_Natural_VecVecs;
-- with Standard_Natural_VecVecs_io;       use Standard_Natural_VecVecs_io;
with Standard_Natural64_VecVecs;
with Standard_Natural64_VecVecs_io;     use Standard_Natural64_VecVecs_io;
--with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;          use Standard_Complex_VecMats;
-- with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Query_Matrices;           use Standard_Query_Matrices;
with Standard_Complex_Polynomials; 
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with Monomial_Hashing;                  use Monomial_Hashing;
with Standard_Deflate_Singularities;
with Standard_Jacobian_Trees;
with Standard_Deflation_Matrices;
with Standard_Evaluate_Deflation;
with Standard_Evaluate_Deflation_io;
with DoblDobl_Deflation_Matrices;
with DoblDobl_Evaluate_Deflation;
with DoblDobl_Evaluate_Deflation_io;
with QuadDobl_Deflation_Matrices;
with QuadDobl_Evaluate_Deflation;
with QuadDobl_Evaluate_Deflation_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Multprec_Random_Numbers;           use Multprec_Random_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;
with Multprec_Complex_VecMats;
with Multprec_Query_Matrices;           use Multprec_Query_Matrices;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Jaco_Matrices;
with Multprec_Jacobian_Trees;
with Multprec_Deflate_Singularities;
--with Multprec_Deflation_Matrices;
with Multprec_Evaluate_Deflation;
with Multprec_Evaluate_Deflation_io;

procedure ts_defmat is

-- DESCRIPTION :
--   This is an attempt to deal efficiently with the "deflation matrices",
--   i.e.: the matrices which arise in the treatment of isolated singularities
--   in Newton's method with deflation.

-- MULTI-LOOP for deriving multi-block multi-variable matrices

  procedure Test_Multi_Loop is

    m : integer32 := 0;

    use Standard_Deflation_Matrices;

    procedure Loop_Body ( index : in Standard_Natural_Vectors.Vector ) is
    begin
      put("Loop body with index = "); put(index); new_line;
    end Loop_Body;
    procedure Run_Loop is new Multi_Loop(Loop_Body);

  begin
    put("Give number of blocks of variables : "); get(m);
    declare
      d,n : Standard_Natural_Vectors.Vector(0..m-1) := (0..m-1 => 0);
      sum : natural32 := 0;
    begin
      for i in 0..m-1 loop
        put("Give number of variables in block "); put(i,1);
        put(" : "); get(n(i));
      end loop;
      for i in 0..m-1 loop
        put("Give order of derivative on block "); put(i,1);
        put(" : "); get(d(i));
      end loop;
      sum := Standard_Natural_Vectors.Sum(d);
      put("Executing loop for s = "); put(sum,1);
      put(", d ="); put(d); put(" and n ="); put(n); put_line(" ...");
      Run_Loop(sum,d,n);
    end;
  end Test_Multi_Loop;

-- DERIVATIVES OF ORIGINAL JACOBIAN : 

  procedure Enumerate_Distinct_Monomials
              ( k,n : in natural32; nb : out natural32 ) is

  -- DESCRIPTION :
  --   Enumerates all distinct monomials of degree k in n variables.
  --   The "distinct" means that a*b is the same as b*a.
  --   On return is the total number of monomials enumerated.

    use Standard_Natural_Vectors;

    cnt : natural32 := 0;

    procedure Write_and_Count ( m : in Standard_Natural_Vectors.Vector ) is

    -- DESCRIPTION :
    --   Writes the monomial m to standard output and increases
    --   the counter by one.

      c : constant natural64 := Monomial_Code(k+1,m);

    begin
      put(m); put("  label : ");
      Standard_Natural_Numbers_io.put(c,1);
      new_line;
      cnt := cnt + 1;
    end Write_and_Count;
    procedure Enumerate is new Enumerate_Monomials(Write_and_Count);

  begin
    Enumerate(k,n);
    nb := cnt;
  end Enumerate_Distinct_Monomials;

  procedure Enumerate_All_Monomials
              ( k,n : in natural32; nb : out natural32 ) is

  -- DESCRIPTION :
  --   Enumerates all monomials of degree k in n variables.
  --   The "all" means that both a*b and b*a are listed.
  --   On return is the total number of monomials enumerated.

  --  use Standard_Natural_Vectors;

    cnt : natural32 := 0;

    procedure Write_and_Count ( m : in Standard_Natural_Vectors.Vector ) is

    -- DESCRIPTION :
    --   Writes the monomial m to standard output and increases
    --   the counter by one.

    begin
      put(m); put("  label : ");
      Standard_Natural_Numbers_io.put(Monomial_Code(k+1,m),1);
      new_line;
      cnt := cnt + 1;
    end Write_and_Count;
    procedure Enumerate is
      new Enumerate_Leaves_of_Monomial_Tree(Write_and_Count);

  begin
    Enumerate(k,n);
    nb := cnt;
  end Enumerate_All_Monomials;

  procedure Test_Hash_Search
              ( monkeys : in Standard_Natural64_VecVecs.VecVec;
                k,n : in natural32; nb : out natural32 ) is

  -- DESCRIPTION :
  --   Enumerates all distinct monomials of degree k in n variables,
  --   and lists their index in the hash table.
  --   On return is the total number of monomials enumerated.

    use Standard_Natural_Vectors;

    cnt : natural32 := 0;

    procedure Write_and_Count ( m : in Vector ) is

    -- DESCRIPTION :
    --   Writes the monomial m to standard output and increases
    --   the counter by one.

      code : constant natural64 := Monomial_Code(k+1,m);
      index : constant natural32 := Search(monkeys,m);

    begin
      put(m);
      put("  label : "); Standard_Natural_Numbers_io.put(code,1);
      put("  index : "); put(index,1); new_line;
      cnt := cnt + 1;
    end Write_and_Count;
    procedure Enumerate is new Enumerate_Monomials(Write_and_Count);

  begin
    Enumerate(k,n);
    nb := cnt;
  end Test_Hash_Search;

  procedure Enumerate_Jacobian_Matrices
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Natural_Vectors;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Jacobian_Trees;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : constant Jaco_Mat(p'range,1..n) := Create(p);
    nd : Node(n) := Initialize(jm);
    k : natural32 := 0;

    procedure Compute_Derivative ( m : Vector ) is

      da : Link_to_Jaco_Mat;

    begin
      Derivative(nd,m,da);
      put("The derivative w.r.t."); put(m); put(" is ");
      if da = null then
        put_line("0");
      else
        new_line;
        for i in da'range(1) loop
          for j in da'range(2) loop
            put("  "); put(da(i,j));
          end loop;
          new_line;
        end loop;
      end if;
    end Compute_Derivative;
    procedure Enumerate is new Enumerate_Monomials(Compute_Derivative);

  begin
    put("Give k : "); get(k);
    Enumerate(k,natural32(n));
  end Enumerate_Jacobian_Matrices;

-- SYMBOLIC CREATION OF DIRECTED ACYCLIC GRAPH OF DERIVATIVE OPERATORS :

  procedure Show_Unwinding_Tree ( k : in natural32 ) is

  -- DESCRIPTION :
  --   Shows the complete tree when unwinding the multipliers in
  --   the deflation matrix of the k-th stage.

    use Standard_Natural_Vectors;
    use Standard_Evaluate_Deflation;
    use Standard_Evaluate_Deflation_io;

    procedure Show ( d : in Vector; i,l : in natural32;
                     zero : in boolean; continue : out boolean ) is
    begin
      if zero
       then Write_Zero(l);
       else Write_Derivative_Operator(d,i,l); new_line;
      end if;
      continue := true;
    end Show;
    procedure Enumerate is new Enumerate_Multiplier_Derivatives(Show);

  begin
    Enumerate(k);
  end Show_Unwinding_Tree;

  procedure Show_Tree_with_Stack ( k : in natural32; size : out natural32 ) is

  -- DESCRIPTION :
  --   Maintains a stack of derivative operators so that not the whole
  --   tree when unwinding the multipliers of A(k) is displayed.
  --   The size of the stack is returned.

    use Standard_Natural_Vectors;
    use Standard_Natural_VecVecs;
    use Standard_Evaluate_Deflation;
    use Standard_Evaluate_Deflation_io;

    max : constant natural32 := 1000; --2*k*k + 1;
    sd : Standard_Natural_VecVecs.VecVec(1..integer32(max));
    sk : Standard_Natural_Vectors.Vector(1..integer32(max));
    sz : natural32 := 0;

    procedure Show ( d : in Vector; i,l : in natural32;
                     zero : in boolean; continue : out boolean ) is

      found : boolean := false;

    begin
      if zero then
        Write_Zero(l); continue := true;
      else
        Write_Derivative_Operator(d,i,l);
        Update_Stack(sd,sk,sz,d,i,found);
        Write_Spaces(i);
        if found
         then put("  already on stack."); continue := false;
         else put("  added to stack."); continue := true;
        end if;
        put("  Stack size = "); put(sz,1); put_line(".");
      end if;
    end Show;
    procedure Enumerate is new Enumerate_Multiplier_Derivatives(Show);

  begin
    Enumerate(k);
    size := sz;
  end Show_Tree_with_Stack;

  procedure Show_Tree_and_Dimensions
              ( k : in natural32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                stack_size,defmat_size,jacmat_size : out natural32 ) is

  -- DESCRIPTION :
  --   Displays the tree of derivative operators and the
  --   dimensions of all the different matrices handled by them.

  -- ON ENTRY :
  --   k        number of deflation stages;
  --   nv       of range 0..k, number of variables in the systems;
  --   nq       of range 0..k, number of equations in the systems;
  --   R1       of range 1..k, number of multiplier variables. 

  -- ON RETURN :
  --   stack_size is the size of the stack;
  --   defmat_size measures the space occupied by the deflation matrices;
  --   jacmat_size measures the space occupied by the Jacobian matrices.

    use Standard_Natural_Vectors;
    use Standard_Natural_VecVecs;
    use Standard_Deflation_Matrices;
    use Standard_Evaluate_Deflation;
    use Standard_Evaluate_Deflation_io;

    max : constant natural32 := 1000; --2*k*k + 1;
    sd : VecVec(1..integer32(max));
    sk : Vector(1..integer32(max));
    sz : natural32 := 0;
    cnt_defmat,cnt_jacmat : natural32 := 0;

    procedure Show ( d : in Vector; i,l : in natural32;
                     zero : in boolean; continue : out boolean ) is

      found : boolean := false;
      inc,moncnt,nc : natural32;

    begin
      if zero then
        Write_Zero(l); continue := true;
      else
        Write_Derivative_Operator(d,i,l);
        Update_Stack(sd,sk,sz,d,i,found);
        Write_Spaces(i);
        put("  stack# = "); put(sz,3);
        if found then
          put_line("  already on stack.");
        else
          put("  ");
          if ((d(d'first) = 0) or (i > 0)) then
            nc := Number_of_Columns(d,nv,R1,i);
            inc := nq(integer32(i))*nc;
            put(nq(integer32(i)),1); put("-by-"); put(nc,1); put(" matrix");
            put(" + "); put(inc,1); new_line;
            cnt_defmat := cnt_defmat + inc;
          else
           -- moncnt := Monomial_Count(d(d'first),nv(0));
            moncnt := nv(0)**integer(d(d'first));
            inc := moncnt*nq(integer32(i))*nv(integer32(i));
            put(moncnt,1); put(" ");
            put(nq(integer32(i)),1); put("-by-");
            put(nv(integer32(i)),1); put(" matrices");
            put(" + "); put(inc,1); new_line;
            cnt_jacmat := cnt_jacmat + inc;
          end if;
        end if;
        continue := not found;
      end if;
    end Show;
    procedure Enumerate is new Enumerate_Multiplier_Derivatives(Show);

  begin
    put("Column numbers : "); put(nv); new_line;
    put("Number of rows : "); put(nq); new_line;
    put("Lambda numbers : "); put(R1); new_line;
    Enumerate(k);
    stack_size := sz;
    defmat_size := cnt_defmat;
    jacmat_size := cnt_jacmat;
  end Show_Tree_and_Dimensions;

  procedure Calculate_Dimensions
              ( k : in natural32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                stack_size,defmat_size,jacmat_size : out natural32 ) is

  -- DESCRIPTION :
  --   Does not show the tree of derivative operators, but enumerates
  --   it to calculate dimensions of all the matrices.

  -- ON ENTRY :
  --   k        number of deflation stages;
  --   nv       of range 0..k, number of variables in the systems;
  --   nq       of range 0..k, number of equations in the systems;
  --   R1       of range 1..k, number of multiplier variables. 

  -- ON RETURN :
  --   stack_size is the size of the stack;
  --   defmat_size measures the space occupied by the deflation matrices;
  --   jacmat_size measures the space occupied by the Jacobian matrices.

    use Standard_Natural_Vectors;
    use Standard_Natural_VecVecs;
    use Standard_Evaluate_Deflation;

    max : constant natural32 := 1000; --2*k*k + 1;
    sd : VecVec(1..integer32(max));
    sk : Vector(1..integer32(max));
    sz : natural32 := 0;
    cnt_defmat,cnt_jacmat : natural32 := 0;

    procedure Enum ( d : in Vector; i,l : in natural32;
                     zero : in boolean; continue : out boolean ) is

      found : boolean := false;
      inc,moncnt : natural32;

    begin
      if zero then
        continue := true;
      else
        Update_Stack(sd,sk,sz,d,i,found);
        if not found then
          if ((d(d'first) = 0) or (i > 0)) then
            inc := nq(integer32(i))*nv(integer32(i));
            cnt_defmat := cnt_defmat + inc;
          else
            moncnt := Monomial_Count(d(d'first),nv(0));
            inc := moncnt*nq(integer32(i))*nv(integer32(i));
            cnt_jacmat := cnt_jacmat + inc;
          end if;
        end if;
        continue := not found;
      end if;
    end Enum;
    procedure Enumerate is new Enumerate_Multiplier_Derivatives(Enum);

  begin
    put("Column numbers : "); put(nv); new_line;
    put("Number of rows : "); put(nq); new_line;
    put("Lambda numbers : "); put(R1); new_line;
    Enumerate(k);
    stack_size := sz;
    defmat_size := cnt_defmat;
    jacmat_size := cnt_jacmat;
  end Calculate_Dimensions;

-- EVALUATION OF DEFLATION MATRICES :

  function Deflate ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                     k : integer32;
                     m : Standard_Natural_Vectors.Vector;
                     B : Standard_Complex_VecMats.VecMat;
                     h : Standard_Complex_VecVecs.VecVec )
                   return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Applies k deflations to the system p, using the random matrices B
  --   and the coefficients in h to scale the multipliers.
  --   The number of multipliers in the i-th stage is m(i).

    use Standard_Complex_Poly_Systems;
    use Standard_Deflate_Singularities;

  begin
    if k = 0 then
      return p;
    elsif k = 1 then
      return Deflate(p,m(1),B(1).all,h(1).all);
    else
      declare
        dp : constant Poly_Sys := Deflate(p,k-1,m,B,h);
      begin
        return Deflate(dp,m(k),B(k).all,h(k).all);
      end;
    end if;
  end Deflate;

  function Deflate ( p : Multprec_Complex_Poly_Systems.Poly_Sys;
                     k : integer32;
                     m : Standard_Natural_Vectors.Vector;
                     B : Multprec_Complex_VecMats.VecMat;
                     h : Multprec_Complex_VecVecs.VecVec )
                   return Multprec_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Applies k deflations to the system p, using the random matrices B
  --   and the coefficients in h to scale the multipliers.
  --   The number of multipliers in the i-th stage is m(i).

    use Multprec_Complex_Poly_Systems;
    use Multprec_Deflate_Singularities;

  begin
    if k = 0 then
      return p;
    elsif k = 1 then
      return Deflate(p,m(1),B(1).all,h(1).all);
    else
      declare
        dp : constant Poly_Sys := Deflate(p,k-1,m,B,h);
      begin
        return Deflate(dp,m(k),B(k).all,h(k).all);
      end;
    end if;
  end Deflate;

  function Standard_Evaluate_Symbolic_Deflation
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32;
                 nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                 B : in Standard_Complex_VecMats.VecMat;
                 h,x : in Standard_Complex_VecVecs.VecVec ) 
               return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   returns the value of the k-th deflation matrix at x.

  -- REQUIRED : nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   p         polynomial system to deflate k times;
  --   k         final deflation stage;
  --   nv        nv(i) is number of columns in i-th deflation matrix;
  --   nq        nq(i) is number of rows in i-th deflation matrix;
  --   R1        R1(i) is the number of multipliers in i-th stage;
  --   B         B(i) is the random matrix used in the i-th stage;
  --   h         h(i) is the random vector to scale the i-th multipliers;
  --   x         x(0) contains values for the original variables,
  --             x(i) contains values for the i-th multipliers.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Jaco_Matrices;

    res : Standard_Complex_Matrices.Matrix
            (1..integer32(nq(k)),1..integer32(nv(k)));
    dp : constant Poly_Sys(1..integer32(nq(k))) := Deflate(p,k,R1,B,h);
    jm : constant Jaco_Mat(1..integer32(nq(k)),1..integer32(nv(k)))
       := Create(dp);
    xx : Standard_Complex_Vectors.Vector(1..integer32(nv(k)));
    ind : integer32 := 0;

  begin
    for i in x'range loop
      for j in x(i)'range loop
        ind := ind+1;
        xx(ind) := x(i)(j);
      end loop;
    end loop;
    res := Eval(jm,xx);
    return res;
  end Standard_Evaluate_Symbolic_Deflation;

  function Multprec_Evaluate_Symbolic_Deflation
               ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32;
                 nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                 B : in Multprec_Complex_VecMats.VecMat;
                 h,x : in Multprec_Complex_VecVecs.VecVec ) 
               return Multprec_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   returns the value of the k-th deflation matrix at x.

  -- REQUIRED : nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   p         polynomial system to deflate k times;
  --   k         final deflation stage;
  --   nv        nv(i) is number of columns in i-th deflation matrix;
  --   nq        nq(i) is number of rows in i-th deflation matrix;
  --   R1        R1(i) is the number of multipliers in i-th stage;
  --   B         B(i) is the random matrix used in the i-th stage;
  --   h         h(i) is the random vector to scale the i-th multipliers;
  --   x         x(0) contains values for the original variables,
  --             x(i) contains values for the i-th multipliers.

    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Jaco_Matrices;

    res : Multprec_Complex_Matrices.Matrix
            (1..integer32(nq(k)),1..integer32(nv(k)));
    dp : constant Poly_Sys(1..integer32(nq(k))) := Deflate(p,k,R1,B,h);
    jm : constant Jaco_Mat(1..integer32(nq(k)),1..integer32(nv(k)))
       := Create(dp);
    xx : Multprec_Complex_Vectors.Vector(1..integer32(nv(k)));
    ind : integer32 := 0;

  begin
    for i in x'range loop
      for j in x(i)'range loop
        ind := ind+1;
        Copy(x(i)(j),xx(ind));
      end loop;
    end loop;
    res := Eval(jm,xx);
    return res;
  end Multprec_Evaluate_Symbolic_Deflation;

  procedure Standard_Interactive_Numeric_Evaluate_Deflation
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32;
                 nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                 B : in Standard_Complex_VecMats.VecMat;
                 h,x : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Calls the function to evaluate the k-th deflation matrix.
  --   This interactive procedure allows user input.

  -- REQUIRED : nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   p         polynomial system to deflate k times;
  --   k         final deflation stage;
  --   nv        nv(i) is number of columns in i-th deflation matrix;
  --   nq        nq(i) is number of rows in i-th deflation matrix;
  --   R1        R1(i) is the number of multipliers in i-th stage;
  --   B         B(i) is the random matrix used in the i-th stage;
  --   h         h(i) is the random vector to scale the i-th multipliers;
  --   x         x(0) contains values for the original variables,
  --             x(i) contains values for the i-th multipliers.

    use Standard_Complex_Jaco_Matrices;
    use Standard_Jacobian_Trees;
    use Standard_Evaluate_Deflation;
    use Standard_Evaluate_Deflation_io;

    eva : Standard_Complex_Matrices.Matrix
            (1..integer32(nq(k)),1..integer32(nv(k)));
    jm : constant Jaco_Mat(p'range,1..integer32(nv(0))) := Create(p);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..k);
    nd : Link_to_Eval_Node;
    evt : Link_to_Eval_Tree;
    ans : character;
    timer : Timing_Widget;
 
  begin
    put("creating a remember table for Jacobian matrices up to order ");
    put(k,1); put_line("...");
    Create_Remember_Derivatives(jm,k,nd);
    monkeys := Monomial_Keys(natural32(k),nv(0));
    put_line("The monomial keys : "); put(monkeys);
    put("creating a remember table for deflation matrices up to order ");
    put(k,1); put_line("...");
    evt := Create(natural32(k));
    put("Do you wish to see the tree of derivative operators ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Write(evt.all,nv,nq,R1);
    end if;
    put("Do you want intermediate output during evaluation ? (y/n) ");
    Ask_Yes_or_No(ans);
    put_line("calling the evaluation routine...");
    tstart(timer);
    if ans = 'y' then
      eva := Eval(standard_output,evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    else
      eva := Eval(evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    end if;
    tstop(timer);
    print_times(standard_output,timer,"numerical evaluation");
    Query_Matrix(eva);
    put("Do you want to compare with symbolic deflation ? (y/n) ");
    Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' then
      declare
        eva_s : Standard_Complex_Matrices.Matrix
                  (1..integer32(nq(k)),1..integer32(nv(k)));
        c_sum : double_float;
      begin
        tstart(timer);
        eva_s := Standard_Evaluate_Symbolic_Deflation(p,k,nv,nq,R1,B,h,x);
        tstop(timer);
        print_times(standard_output,timer,"symbolic evaluation");
        c_sum := Difference_of_Matrices(eva,eva_s);
        put("Componentwise Difference Sum : "); put(c_sum,3); new_line;
        put("Do you wish to see the difference element by element ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y'
         then Show_Differences_in_Matrices(eva,eva_s);
        end if;
        Query_Matrices(eva,eva_s);
      end;
    end if;
  end Standard_Interactive_Numeric_Evaluate_Deflation;

  procedure Multprec_Interactive_Numeric_Evaluate_Deflation
               ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32;
                 nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                 B : in Multprec_Complex_VecMats.VecMat;
                 h,x : in Multprec_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Calls the function to evaluate the k-th deflation matrix.
  --   This interactive procedure allows user input.

  -- REQUIRED : nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   p         polynomial system to deflate k times;
  --   k         final deflation stage;
  --   nv        nv(i) is number of columns in i-th deflation matrix;
  --   nq        nq(i) is number of rows in i-th deflation matrix;
  --   R1        R1(i) is the number of multipliers in i-th stage;
  --   B         B(i) is the random matrix used in the i-th stage;
  --   h         h(i) is the random vector to scale the i-th multipliers;
  --   x         x(0) contains values for the original variables,
  --             x(i) contains values for the i-th multipliers.

    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Jacobian_Trees;
    use Multprec_Evaluate_Deflation;
    use Multprec_Evaluate_Deflation_io;

    eva : Multprec_Complex_Matrices.Matrix
            (1..integer32(nq(k)),1..integer32(nv(k)));
    jm : constant Jaco_Mat(p'range,1..integer32(nv(0))) := Create(p);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..k);
    nd : Link_to_Eval_Node;
    evt : Link_to_Eval_Tree;
    ans : character;
 
  begin
    put("creating a remember table for Jacobian matrices up to order ");
    put(k,1); put_line("...");
    Create_Remember_Derivatives(jm,k,nd);
    monkeys := Monomial_Keys(natural32(k),nv(0));
    put_line("The monomial keys : "); put(monkeys);
    put("creating a remember table for deflation matrices up to order ");
    put(k,1); put_line("...");
    evt := Create(natural32(k));
    put("Do you wish to see the tree of derivative operators ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Write(evt.all,nv,nq,R1);
    end if;
    put("Do you want intermediate output during evaluation ? (y/n) ");
    Ask_Yes_or_No(ans);
    put_line("calling the evaluation routine...");
    if ans = 'y' then
      eva := Eval(standard_output,evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    else
      eva := Eval(evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    end if;
    Query_Matrix(eva);
    put("Do you want to compare with symbolic deflation ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        eva_s : constant Multprec_Complex_Matrices.Matrix
                  (1..integer32(nq(k)),1..integer32(nv(k)))
              := Multprec_Evaluate_Symbolic_Deflation(p,k,nv,nq,R1,B,h,x);
        c_sum : constant Floating_Number := Difference_of_Matrices(eva,eva_s);
      begin
        put("Componentwise Difference Sum : "); put(c_sum,3); new_line;
        put("Do you wish to see the difference element by element ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y'
         then Show_Differences_in_Matrices(eva,eva_s);
        end if;
        Query_Matrices(eva,eva_s);
      end;
    end if;
  end Multprec_Interactive_Numeric_Evaluate_Deflation;

  procedure Compare_Symbolic_Numeric_Evaluation
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32;
                 nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                 B : in Standard_Complex_VecMats.VecMat;
                 h,x : in Standard_Complex_VecVecs.VecVec;
                 check_sum : out double_float ) is

  -- DESCRIPTION :
  --   Compares the value of the k-th deflation matrix, obtained twice,
  --   once symbolically, and once numerically.
  --   No input of the user is required.

  -- REQUIRED : nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   p         polynomial system to deflate k times;
  --   k         final deflation stage;
  --   nv        nv(i) is number of columns in i-th deflation matrix;
  --   nq        nq(i) is number of rows in i-th deflation matrix;
  --   R1        R1(i) is the number of multipliers in i-th stage;
  --   B         B(i) is the random matrix used in the i-th stage;
  --   h         h(i) is the random vector to scale the i-th multipliers;
  --   x         x(0) contains values for the original variables,
  --             x(i) contains values for the i-th multipliers.
 
    use Standard_Complex_Jaco_Matrices;
    use Standard_Jacobian_Trees;
    use Standard_Evaluate_Deflation;

    num_val : Standard_Complex_Matrices.Matrix
                (1..integer32(nq(k)),1..integer32(nv(k)));
    sym_val : Standard_Complex_Matrices.Matrix
                (1..integer32(nq(k)),1..integer32(nv(k)));
    jm : constant Jaco_Mat(p'range,1..integer32(nv(0))) := Create(p);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..k);
    nd : Link_to_Eval_Node;
    evt : Link_to_Eval_Tree;
 
  begin
    Create_Remember_Derivatives(jm,k,nd);
    monkeys := Monomial_Keys(natural32(k),nv(0));
    evt := Create(natural32(k));
    num_val := Eval(evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    sym_val := Standard_Evaluate_Symbolic_Deflation(p,k,nv,nq,R1,B,h,x);
    check_sum := Difference_of_Matrices(num_val,sym_val);
  end Compare_Symbolic_Numeric_Evaluation;

  procedure Compare_Symbolic_Numeric_Evaluation
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32;
                 nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                 B : in Standard_Complex_VecMats.VecMat;
                 h,x : in Standard_Complex_VecVecs.VecVec;
                 check_sum : out double_float ) is

  -- DESCRIPTION :
  --   Compares the value of the k-th deflation matrix, obtained twice,
  --   once symbolically, and once numerically.
  --   No input of the user is required.

  -- REQUIRED : nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   file      for intermediate results and diagnostics;
  --   p         polynomial system to deflate k times;
  --   k         final deflation stage;
  --   nv        nv(i) is number of columns in i-th deflation matrix;
  --   nq        nq(i) is number of rows in i-th deflation matrix;
  --   R1        R1(i) is the number of multipliers in i-th stage;
  --   B         B(i) is the random matrix used in the i-th stage;
  --   h         h(i) is the random vector to scale the i-th multipliers;
  --   x         x(0) contains values for the original variables,
  --             x(i) contains values for the i-th multipliers.

    use Standard_Complex_Jaco_Matrices;
    use Standard_Jacobian_Trees;
    use Standard_Evaluate_Deflation;

    num_val : Standard_Complex_Matrices.Matrix
                (1..integer32(nq(k)),1..integer32(nv(k)));
    sym_val : Standard_Complex_Matrices.Matrix
                (1..integer32(nq(k)),1..integer32(nv(k)));
    jm : constant Jaco_Mat(p'range,1..integer32(nv(0))) := Create(p);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..k);
    nd : Link_to_Eval_Node;
    evt : Link_to_Eval_Tree;
 
  begin
    put(file,"creating a remember table for Jacobian matrices up to order ");
    put(file,k,1); put_line(file,"...");
    Create_Remember_Derivatives(jm,k,nd);
    monkeys := Monomial_Keys(natural32(k),nv(0));
    put_line(file,"The monomial keys : "); put(file,monkeys);
    put(file,"creating a remember table for deflation matrices up to order ");
    put(file,k,1); put_line(file,"...");
    evt := Create(natural32(k));
    put_line(file,"calling the numerical evaluation routine...");
    num_val := Eval(file,evt,nd,monkeys,k,nv,nq,R1,B,h,x);
    put_line(file,"calling the symbolic evaluation routine...");
    sym_val := Standard_Evaluate_Symbolic_Deflation(p,k,nv,nq,R1,B,h,x);
    check_sum := Difference_of_Matrices(num_val,sym_val);
    put(file,"Componentwise Difference Sum : ");
    put(file,check_sum,3); new_line(file);
    Show_Differences_in_Matrices(num_val,sym_val);
  end Compare_Symbolic_Numeric_Evaluation;

  procedure Read_Dimensions
             ( k,nv0,nq0 : in natural32;
               nv,nq,R1 : out Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Asks the user for the dimensions in the deflation procedure.

  -- ON ENTRY :
  --   k       number of deflation stages;
  --   nv0     number of variables in the original system;
  --   nq0     number of equations in the original system, nq0 >= nv0.

  -- REQUIRED : nv'range = 0..k, nq'range = 0..k, R1'range = 1..k.

  -- ON RETURN :
  --   nv      nv(i) is the number of variables at the i-th deflation;
  --   nq      nq(i) is the number of equations at the i-th deflation;
  --   R1      R1(i) is the number of multipliers at the i-th deflation,
  --           it is one plus the rank of the (i-1)-th deflation matrix.

  begin
    nv(0) := nv0;
    nq(0) := nq0;
    for i in 1..integer32(k) loop
      loop
        put("Give rank ("); put("<"); put(nv(i-1),1); put(")");
	put(" of deflation matrix "); put(i-1,1); put(" : ");
        R1(i) := 0; get(R1(i));
        R1(i) := R1(i) + 1;
        exit when (R1(i) <= nv(i-1));
        put("Rank must be smaller than "); put(nv(i-1),1);
        put_line(".  Please try again...");
      end loop;
      nv(i) := nv(i-1) + R1(i);
      nq(i) := 2*nq(i-1) + 1;
    end loop;
  end Read_Dimensions;

  procedure Standard_Random_Coefficients
              ( k : in natural32; output : in boolean;
                nv,R1 : in Standard_Natural_Vectors.Vector;
                B : out Standard_Complex_VecMats.VecMat;
                h : out Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Generates all random B's and h's used in deflating k times.

  -- ON ENTRY :
  --   k        number of deflation stages;
  --   nv       nv(i) is the number of variables in i-th system;
  --   R1       R1(i) is the number of multipliers in i-th stage.

  -- REQUIRED : nv'range = 0..k, R1'range = h'range = B'range = 1..k.

  -- ON RETURN :
  --   h        h(i) is used to scale the multipliers of the i-th stage;
  --   B        B(i) is used in the i-th deflation stage. 

  begin
    for i in 1..integer32(k) loop
      if output then
        put("h("); put(i,1); put(") is random ");
        put(R1(i),1); put_line("-vector");
        put("B("); put(i,1); put(") is random ");
        put(nv(i-1),1); put("-by-"); put(R1(i),1); put_line("-matrix");
      end if;
      declare
        hi : Standard_Complex_Vectors.Vector(1..integer32(R1(i)));
        Bi : Standard_Complex_Matrices.Matrix
               (1..integer32(nv(i-1)),1..integer32(R1(i)));
      begin
        for j in hi'range loop
          hi(j) := Random1;
        end loop;
        h(i) := new Standard_Complex_Vectors.Vector'(hi);
        for i1 in Bi'range(1) loop
          for j1 in Bi'range(2) loop
            Bi(i1,j1) := Random1;
          end loop;
        end loop;
        B(i) := new Standard_Complex_Matrices.Matrix'(Bi);
      end;
    end loop;
  end Standard_Random_Coefficients;

  procedure Multprec_Random_Coefficients
              ( size,k : in natural32; output : in boolean;
                nv,R1 : in Standard_Natural_Vectors.Vector;
                B : out Multprec_Complex_VecMats.VecMat;
                h : out Multprec_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Generates all random B's and h's used in deflating k times.

  -- ON ENTRY :
  --   k        number of deflation stages;
  --   nv       nv(i) is the number of variables in i-th system;
  --   R1       R1(i) is the number of multipliers in i-th stage.

  -- REQUIRED : nv'range = 0..k, R1'range = h'range = B'range = 1..k.

  -- ON RETURN :
  --   h        h(i) is used to scale the multipliers of the i-th stage;
  --   B        B(i) is used in the i-th deflation stage. 

  begin
    for i in 1..integer32(k) loop
      if output then
        put("h("); put(i,1); put(") is random ");
        put(R1(i),1); put_line("-vector");
        put("B("); put(i,1); put(") is random ");
        put(nv(i-1),1); put("-by-"); put(R1(i),1); put_line("-matrix");
      end if;
      declare
        hi : Multprec_Complex_Vectors.Vector(1..integer32(R1(i)));
        Bi : Multprec_Complex_Matrices.Matrix
               (1..integer32(nv(i-1)),1..integer32(R1(i)));
      begin
        for j in hi'range loop
          hi(j) := Random(size);
        end loop;
        h(i) := new Multprec_Complex_Vectors.Vector'(hi);
        for i1 in Bi'range(1) loop
          for j1 in Bi'range(2) loop
            Bi(i1,j1) := Random(size);
          end loop;
        end loop;
        B(i) := new Multprec_Complex_Matrices.Matrix'(Bi);
      end;
    end loop;
  end Multprec_Random_Coefficients;

  function Standard_Random_Evaluation_Data
               ( output : boolean; k : integer32; n : natural32;
                 R1 : Standard_Natural_Vectors.Vector )
               return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Generates random numbers for the input variables
  --   to evaluate the deflation matrices of order k.

  -- REQUIRED : R1'range = 1..k, vector on return has range 0..k.

  -- ON ENTRY :
  --   output    for writing intermediate output;
  --   k         order of the deflation;
  --   n         original number of variables;
  --   R1        R1(i) is number of multiplier variables in i-th stage.

  -- ON RETURN :
  --   array of range 0..k, the 0-th one is a random n-vector,
  --   while the i-th one is a random R1(i)-vector.

    x : Standard_Complex_VecVecs.VecVec(0..k);

  begin
    if output then
      put("generating "); put(n,1);
      put_line(" random values for the original values...");
    end if;
    declare -- will be x-values
      x0 : Standard_Complex_Vectors.Vector(1..integer32(n));
    begin
      for j in x0'range loop
        x0(j) := Random1;
      end loop;
      x(0) := new Standard_Complex_Vectors.Vector'(x0);
    end;
    for i in 1..k loop    -- generate values for i-th multipliers
      if output then
        put("generating "); put(R1(i),1);
        put(" random variables for multiplier vector "); put(i,1);
        put_line("...");
      end if;
      declare
        xi : Standard_Complex_Vectors.Vector(1..integer32(R1(i)));
      begin
        for j in xi'range loop
          xi(j) := Random1;
        end loop;
        x(i) := new Standard_Complex_Vectors.Vector'(xi);
      end;
    end loop;
    return x;
  end Standard_Random_Evaluation_Data;

  function Multprec_Random_Evaluation_Data
               ( output : boolean; size : natural32;
                 k : integer32; n : natural32;
                 R1 : Standard_Natural_Vectors.Vector )
               return Multprec_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Generates random numbers for the input variables
  --   to evaluate the deflation matrices of order k.

  -- REQUIRED : R1'range = 1..k, vector on return has range 0..k.

  -- ON ENTRY :
  --   output    for writing intermediate output;
  --   size      size of the multiprecision numbers;
  --   k         order of the deflation;
  --   n         original number of variables;
  --   R1        R1(i) is number of multiplier variables in i-th stage.

  -- ON RETURN :
  --   array of range 0..k, the 0-th one is a random n-vector,
  --   while the i-th one is a random R1(i)-vector.

    x : Multprec_Complex_VecVecs.VecVec(0..k);

  begin
    if output then
      put("generating "); put(n,1);
      put_line(" random values for the original values...");
    end if;
    declare -- will be x-values
      x0 : Multprec_Complex_Vectors.Vector(1..integer32(n));
    begin
      for j in x0'range loop
        x0(j) := Random(size);
      end loop;
      x(0) := new Multprec_Complex_Vectors.Vector'(x0);
    end;
    for i in 1..k loop   -- generate values for i-th multipliers
      if output then
        put("generating "); put(R1(i),1);
        put(" random variables for multiplier vector "); put(i,1);
        put_line("...");
      end if;
      declare
        xi : Multprec_Complex_Vectors.Vector(1..integer32(R1(i)));
      begin
        for j in xi'range loop
          xi(j) := Random(size);
        end loop;
        x(i) := new Multprec_Complex_Vectors.Vector'(xi);
      end;
    end loop;
    return x;
  end Multprec_Random_Evaluation_Data;

  procedure Standard_Symbolic_Deflation_Matrix
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32 ) is

  -- DESCRIPTION :
  --   Applies symbolic deflation to compute the k-th deflation matrix.

    use Standard_Complex_Polynomials;

    nv0 : constant natural32 := Number_of_Unknowns(p(p'first));
    nq0 : constant natural32 := natural32(p'last(1));
    nv,nq : Standard_Natural_Vectors.Vector(0..k);
    R1 : Standard_Natural_Vectors.Vector(1..k);
    B : Standard_Complex_VecMats.VecMat(1..k);
    h : Standard_Complex_VecVecs.VecVec(1..k);
    x : Standard_Complex_VecVecs.VecVec(0..k);

  begin
    Read_Dimensions(natural32(k),nv0,nq0,nv,nq,R1);
    put("Column numbers : "); put(nv); new_line;
    put("Number of rows : "); put(nq); new_line;
    put("Lambda numbers : "); put(R1); new_line;
    Standard_Random_Coefficients(natural32(k),true,nv,R1,B,h);
    x := Standard_Random_Evaluation_Data(true,k,nv(0),R1);
    declare
      eva : constant Standard_Complex_Matrices.Matrix
              (1..integer32(nq(k)),1..integer32(nv(k)))
          := Standard_Evaluate_Symbolic_Deflation(p,k,nv,nq,R1,B,h,x);
    begin
      Query_Matrix(eva);
    end;
  end Standard_Symbolic_Deflation_Matrix;

  procedure Standard_Evaluate_Deflation_Matrix
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                k : in integer32 ) is

    use Standard_Complex_Polynomials;

    nv0 : constant natural32 := Number_of_Unknowns(p(p'first));
    nq0 : constant natural32 := natural32(p'last(1));
    nv,nq : Standard_Natural_Vectors.Vector(0..k);
    R1 : Standard_Natural_Vectors.Vector(1..k);
    h : Standard_Complex_VecVecs.VecVec(1..k);
    B : Standard_Complex_VecMats.VecMat(1..k);
    x : Standard_Complex_VecVecs.VecVec(0..k);

  begin
    Read_Dimensions(natural32(k),nv0,nq0,nv,nq,R1);
    put("Column numbers : "); put(nv); new_line;
    put("Number of rows : "); put(nq); new_line;
    put("Lambda numbers : "); put(R1); new_line;
    Standard_Random_Coefficients(natural32(k),true,nv,R1,B,h);
    x := Standard_Random_Evaluation_Data(true,k,nv(0),R1);
    Standard_Interactive_Numeric_Evaluate_Deflation(p,k,nv,nq,R1,B,h,x);
  end Standard_Evaluate_Deflation_Matrix;

  procedure Multprec_Evaluate_Deflation_Matrix
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                size : in natural32; k : in integer32 ) is

    use Multprec_Complex_Polynomials;

    nv0 : constant natural32 := Number_of_Unknowns(p(p'first));
    nq0 : constant natural32 := natural32(p'last(1));
    nv,nq : Standard_Natural_Vectors.Vector(0..k);
    R1 : Standard_Natural_Vectors.Vector(1..k);
    h : Multprec_Complex_VecVecs.VecVec(1..k);
    B : Multprec_Complex_VecMats.VecMat(1..k);
    x : Multprec_Complex_VecVecs.VecVec(0..k);

  begin
    Read_Dimensions(natural32(k),nv0,nq0,nv,nq,R1);
    put("Column numbers : "); put(nv); new_line;
    put("Number of rows : "); put(nq); new_line;
    put("Lambda numbers : "); put(R1); new_line;
    Multprec_Random_Coefficients(size,natural32(k),true,nv,R1,B,h);
    x := Multprec_Random_Evaluation_Data(true,size,k,nv(0),R1);
    Multprec_Interactive_Numeric_Evaluate_Deflation(p,k,nv,nq,R1,B,h,x);
  end Multprec_Evaluate_Deflation_Matrix;

  generic
    with procedure Variables ( m : in Standard_Natural_Vectors.Vector );
    -- m(0) equals the number of variables in the original system
    -- m(i) is the rank of Jacobian matrix in the i-th deflation stage
  procedure Generate_Deflations ( n,k : in natural32 );

  -- DESCRIPTION :
  --   Generates all ways one can deflate a system with n variables
  --   successively k times.

  procedure Generate_Deflations ( n,k : in natural32 ) is

    accu : Standard_Natural_Vectors.Vector(0..integer32(k));

    procedure Enumerate ( i,s : in integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all possible choices for the i-th rank,
    --   for a system of s variables, when i <= k.
    --   s is the sum of the first i entries in accu.

    begin
      if i > integer32(k) then
        Variables(accu);
      else
        for j in 0..s-1 loop
          accu(i) := natural32(j);
          Enumerate(i+1,s+j+1);
        end loop;
      end if;
    end Enumerate;

  begin
    accu(0) := n;
    Enumerate(1,integer32(n));
  end Generate_Deflations;

  procedure Enumerate_Evaluations
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                k : in integer32 ) is

    use Standard_Complex_Polynomials;

    nv0 : constant natural32 := Number_of_Unknowns(p(p'first));
    nq0 : constant natural32 := natural32(p'last(1));
    cnt : natural32 := 0;
    sum : double_float := 0.0;
   -- ans : character;
    output : constant boolean := true;

    procedure Test_Deflation ( m : in Standard_Natural_Vectors.Vector ) is

      nv,nq : Standard_Natural_Vectors.Vector(0..k);
      R1 : Standard_Natural_Vectors.Vector(1..k);
      h : Standard_Complex_VecVecs.VecVec(1..k);
      B : Standard_Complex_VecMats.VecMat(1..k);
      x : Standard_Complex_VecVecs.VecVec(0..k);
      check_sum : double_float;

    begin
      put("Testing deflation for m ="); put(m,1);
      nv(0) := m(0);
      nq(0) := natural32(p'last(1));
      for i in 1..k loop
        R1(i) := m(i)+1;            -- #multipliers is rank + 1
        nv(i) := nv(i-1) + R1(i);   -- #variables in i-th stage
        nq(i) := 2*nq(i-1) + 1;     -- #equations in i-th stage
      end loop;
      if output then
        new_line;
        put("Column numbers : "); put(nv); new_line;
        put("Number of rows : "); put(nq); new_line;
        put("Lambda numbers : "); put(R1); new_line;
      end if;
      Standard_Random_Coefficients(natural32(k),output,nv,R1,B,h);
      x := Standard_Random_Evaluation_Data(output,k,nv(0),R1);
      if output then
        Compare_Symbolic_Numeric_Evaluation
          (standard_output,p,k,nv,nq,R1,B,h,x,check_sum);
      else
        Compare_Symbolic_Numeric_Evaluation(p,k,nv,nq,R1,B,h,x,check_sum);
      end if;
      put(" check sum = "); put(check_sum,3); new_line;
      sum := sum + check_sum;
      cnt := cnt + 1;
    end Test_Deflation;
    procedure Test_All is new Generate_Deflations(Test_Deflation);

  begin
   -- put("Do you wish intermediate output ? (y/n) ");
   -- Ask_Yes_or_No(ans);
   -- output := (ans = 'y');
    Test_All(nv0,natural32(k));
    put("Tested "); put(cnt,1); put_line(" cases.");
    put("Sum of all difference checks : "); put(sum,3); new_line;
  end Enumerate_Evaluations;

  procedure Write_Bytes ( n : in natural32 ) is

    n0 : constant natural32 := n mod 1024;
    n1 : natural32 := n/1024;
    n2 : natural32;

  begin
    if n1 > 1000 then
      n2 := n1/1024;
      n1 := n1 mod 1024;
      put(n2,1); put(",");
      put(n1,3); put(","); put(n0,3); put(" bytes");
    else
      put(n1,1); put(","); put(n0,3); put(" bytes");
    end if;
  end Write_Bytes;

-- MAIN TEST PROGRAMS :

  function Take_Order_from_Menu return character is

  -- DESCRIPTION :
  --   Shows the menu to the user and returns the choice made,
  --   as a character representing a number between 0 and 5.

    ans : character;

  begin
    new_line;
    put_line("Interactive test on deflation matrices ...");
    new_line;
    put_line("Choose one of the following :");
    put_line("  0. Show the complete tree unwinding the multipliers;");
    put_line("  1. Use remember table to display the unwinding tree;");
    put_line("  2. Write unwinding tree with dimension calculations;");
    put_line("  3. Enumerate derivative operators and test hashing;");
    put_line("  4. Enumerate derivatives of original Jacobian matrix;");
    put_line("  5. Test evaluation of symbolic deflation matrix;");
    put_line("  6. Test the running of a multi-loop;");
    put_line("  7. Test evaluation of numerical deflation matrix;");
    put_line("  8. Enumerate all possible numerical evaluations;");
    put_line("  9. Multiprecision evaluation of a deflation matrix.");
    put("Type 0, 1, 2, 3, 4, 5, 6, 7, 8, or 9 to make your choice : ");
    Ask_Alternative(ans,"0123456789");
    return ans;
  end Take_Order_from_Menu;

  procedure Show_Complete_Unwinding_Tree is

  -- DESCRIPTION :
  --   Asks the user for k, the deflation stage, and then displays
  --   the complete tree to unwind the multipliers.

    k : natural32 := 0;

  begin
    put("Give deflation stage k : "); get(k);
    new_line;
    Show_Unwinding_Tree(k);
  end Show_Complete_Unwinding_Tree;

  procedure Show_Directed_Acyclic_Graph is

  -- DESCRIPTION :
  --   When a stack is used to store every node only once, the tree to
  --   unwind the multipliers becomes a directed acyclic graph.

    use Standard_Evaluate_Deflation;
    use Standard_Evaluate_Deflation_io;

    k,cnt : natural32 := 0;
    evt : Link_to_Eval_Tree;

  begin
    put("Give deflation stage k : "); get(k);
    new_line;
    Show_Tree_with_Stack(k,cnt);
    new_line;
    put("Number of symbols in remember table : ");
    put(cnt,1); put_line(".");
    evt := Create(k);
    Write(evt.all);
    put("Number of nodes : ");
    put(Node_Count(evt.all),1); new_line;
    put("Number of different nodes : ");
    put(Different_Node_Count(evt.all),1); new_line;
  end Show_Directed_Acyclic_Graph;

  procedure Dimensions_of_Directed_Acyclic_Graph is

  -- DESCRIPTION :
  --   This procedure calculates the memory consumed for particular
  --   dimensions of an input system with a sequence of ranks.

    k,n,m,cnt_dm,cnt_jm,ssz : natural32 := 0;
    ans : character;

  begin
    put("Give deflation stage k : "); get(k);
    put("Give the number of variables : "); get(n);
    put("Give the number of equations : "); get(m);
    declare
      nv,nq : Standard_Natural_Vectors.Vector(0..integer32(k));
      R1 : Standard_Natural_Vectors.Vector(1..integer32(k));
    begin
      Read_Dimensions(k,n,m,nv,nq,R1);
      put("Do you wish to see the tree ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Show_Tree_and_Dimensions(k,nv,nq,R1,ssz,cnt_dm,cnt_jm);
       else Calculate_Dimensions(k,nv,nq,R1,ssz,cnt_dm,cnt_jm);
      end if;         
    end;
    put("#different operators : "); put(ssz,1); new_line;
    put("#numbers in deflation matrices : ");
    put(cnt_dm,1); new_line;
    put(" 1 complex = 8 bytes : "); Write_Bytes(8*cnt_dm);
    new_line;
    put("#numbers in Jacobian matrices : ");
    put(cnt_jm,1); new_line;
    put(" 1 complex = 8 bytes : "); Write_Bytes(8*cnt_jm);
    new_line;
  end Dimensions_of_Directed_Acyclic_Graph;

  procedure Test_Monomial_Hashing is

    k,n,cnt0,cnt1,cnt2 : natural32 := 0;

  begin
    put("Give degree of derivative : "); get(k);
    put("Give the number of variables : "); get(n);
    put_line("All distinct monomials : ");
    Enumerate_Distinct_Monomials(k,n,cnt0);
    put_line("Enumeration of all monomials :");
    Enumerate_All_Monomials(k,n,cnt1);
    put("The number of distinct monomials : ");
    put(cnt0,1); new_line;
    cnt2 := Monomial_Count(k,n);
    put("Counting the distinct monomials again : ");
    put(cnt2,1); new_line;
    put("The number of duplicate monomials : ");
    put(cnt1,1); new_line;
    put_line("The monomial codes :");
    declare
      monkeys : constant Standard_Natural64_VecVecs.VecVec(1..integer32(k))
              := Monomial_Keys(k,n);
    begin
      for i in monkeys'range loop
        put(monkeys(i).all); new_line;
      end loop;
      put_line("The indices to the monomials :");
      Test_Hash_Search(monkeys,k,n,cnt0);
    end;
  end Test_Monomial_Hashing;

  procedure Main is

    ans : character;
    size : natural32 := 0;
    k : integer32 := 0;
    stlp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mplp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    ans := Take_Order_from_Menu;
    new_line;
    case ans is
      when '0' => Show_Complete_Unwinding_Tree;
      when '1' => Show_Directed_Acyclic_Graph;
      when '2' => Dimensions_of_Directed_Acyclic_Graph;
      when '3' => Test_Monomial_Hashing;
      when '4' => get(stlp);
                  Enumerate_Jacobian_Matrices(stlp.all);
      when '5' => get(stlp);
                  new_line;
                  put("Give deflation stage k : "); get(k);
                  Standard_Symbolic_Deflation_Matrix(stlp.all,k);
      when '6' => Test_Multi_Loop;
      when '7' | '8'
               => get(stlp);
                  new_line;
                  put("Give deflation stage k : "); get(k);
                  if ans = '7'
                   then Standard_Evaluate_Deflation_Matrix(stlp.all,k);
                   else Enumerate_Evaluations(stlp.all,k);
                  end if;
      when '9' => get(mplp);
                  new_line;
                  put("Give the size of the numbers : "); get(size);
                  put("Give deflation stage k : "); get(k);
                  Multprec_Evaluate_Deflation_Matrix(mplp.all,size,k);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_defmat;
