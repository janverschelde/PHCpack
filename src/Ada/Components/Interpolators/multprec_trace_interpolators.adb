with unchecked_deallocation;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Complex_VecVecs;           use Multprec_Complex_VecVecs;
with Multprec_Complex_VecVecs_io;        use Multprec_Complex_VecVecs_io;
with Multprec_Complex_Matrices;          use Multprec_Complex_Matrices;
with Multprec_Complex_NesVecs;           use Multprec_Complex_NesVecs;
with Multprec_Complex_NesVecs_io;        use Multprec_Complex_NesVecs_io;
with Multprec_Complex_Poly_Systems;      use Multprec_Complex_Poly_Systems;
with Sample_Points;                      use Sample_Points;
with Rectangular_Sample_Grids;           use Rectangular_Sample_Grids;
with Multprec_Univariate_Interpolators;  use Multprec_Univariate_Interpolators;
with Multprec_Nvariate_Interpolators;    use Multprec_Nvariate_Interpolators;
with Multprec_Power_Traces;              use Multprec_Power_Traces;
with Hybrid_Durand_Kerner;               use Hybrid_Durand_Kerner;

package body Multprec_Trace_Interpolators is

-- DATA STRUCTURES :

  type Trace_Interpolator1_Rep ( n,d : integer32 ) is record
    rot : Vector(1..n);   -- normal to all parallel hyperplanes in grid
    trc : VecVec(1..d);   -- trc(i) contains coefficients of trace of deg i
  end record;

  type Newton_Form ( nv : natural32; n,d : integer32 ) is record -- nv = n
    x : VecVec(1..n);    -- grid is product of n vectors of range 0..d
    f : NesVec(nv,0,d);  -- divided differences for Newton interpolator
  end record;
  type Link_to_Newton_Form is access Newton_Form;
  type Newton_Forms is array ( integer32 range <> ) of Link_to_Newton_Form;

  type Trace_Interpolator_Rep ( k,n,d : integer32 ) is record
    case k is  -- k = #slices; n = #variables; d = degree
      when 1 => t1 : Trace_Interpolator1_Rep(n,d);
      when others => rot : VecVec(1..k);        -- normals to slices
                     trc : Newton_Forms(1..d);  -- low dim trace forms
    end case;
  end record;

-- AUXILIARIES TO THE CREATOR :

  function Power_Sum ( start,i,j : integer32;
                       m : Matrix ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the i-th power sum at the j-th slice, using the rotated
  --   and projected 2nd coordinate of the roots, starting with roots
  --   at the index start.

    res : Complex_Number := Create(integer(0));
    acc : Complex_Number;

  begin
    if i = 1 then
      for k in start..m'last(1) loop
        Add(res,m(k,j));
      end loop;
    else
      for k in start..m'last(1)-i+1 loop
        acc := Power_Sum(k+1,i-1,j,m);
        Mul(acc,m(k,j));
        Add(res,acc);
        Clear(acc);
      end loop;
    end if;
    return res;
  end Power_Sum;

  function Power_Sum ( file : file_type; start,i,j : integer32;
                       m : Matrix ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the i-th power sum at the j-th slice, using the rotated
  --   and projected 2nd coordinate of the roots, starting with roots
  --   at the index start.

    res : Complex_Number := Create(integer(0));
    acc : Complex_Number;

  begin
    if i = 1 then
      for k in start..m'last(1) loop
        put(file," + m("); put(file,k,1); put(file,",");
        put(file,j,1); put(file,")");
        Add(res,m(k,j));
      end loop;
      new_line(file);
    else
      for k in start..m'last(1)-i+1 loop
        put(file," + m("); put(file,k,1); put(file,",");
        put(file,j,1); put(file,")*");
        acc := Power_Sum(k+1,i-1,j,m);
        Mul(acc,m(k,j));
        Add(res,acc);
        Clear(acc);
      end loop;
      new_line(file);
    end if;
    return res;
  end Power_Sum;

  procedure Adjust_Signs ( m : in out Matrix ) is

  -- DESCRIPTION :
  --   Multiplies the i-th row of the matrix with (-1)^n.

    sign : boolean := (m'last(1) mod 2 = 1);

  begin
    for i in reverse m'range(1) loop
      if sign then
        for j in m'range(2) loop
          Min(m(i,j));
        end loop;
      end if;
      sign := not sign;
    end loop;
  end Adjust_Signs;

  function Power_Sums ( i : integer32; m : Matrix ) return Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..i with all power sums of degree i,
  --   over the columns of m, assuming the j-th columns contains the
  --   rotated and projected samples of slice j, for j in 0..i.

    res : Vector(0..i);

  begin
    for j in 0..i loop
      res(j) := Power_Sum(1,i,j,m);
    end loop;
    if i mod 2 = 1
     then Min(res);
    end if;
    return res;
  end Power_Sums;

  function Power_Sums ( m : Matrix ) return Matrix is

  -- DESCRIPTION :
  --   Given the rotated and projected 2nd coordinates of the roots,
  --   the matrix on return contains the power sums of the roots in m.
  --   If ps is the name of the matrix on return, then ps(i,j)
  --   contains the power sum of degree i evaluated at the j-th slice.

    res : Matrix(m'range(1),m'range(2));

  begin
    for i in res'range(1) loop
      for j in 0..i loop
        res(i,j) := Power_Sum(1,i,j,m);
      end loop;
    end loop;
    Adjust_Signs(res);
    return res;
  end Power_Sums;

  function Power_Sums ( file : file_type; m : Matrix ) return Matrix is

  -- DESCRIPTION :
  --   Given the rotated and projected 2nd coordinates of the roots,
  --   the matrix on return contains the power sums of the roots in m.
  --   If ps is the name of the matrix on return, then ps(i,j)
  --   contains the power sum of degree i evaluated at the j-th slice.

    res : Matrix(m'range(1),m'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Power_Sum(file,1,i,j,m);
      end loop;
    end loop;
    Adjust_Signs(res);
    return res;
  end Power_Sums;

  function Interpolate ( x : Vector; m : Matrix ) return VecVec is

  -- DESCRIPTION :
  --   Given the function values needed for the traces in the matrix m,
  --   the expanded coefficient vectors of the trace forms are returned.
  --   Note that the order of indexing is reversed: the i-th index
  --   of the vector on return contains trace d-i.

    d : constant integer32 := m'last(1);
    res : VecVec(1..d);

  begin
    for i in 1..d loop
      declare
        y,f : Vector(0..i);
      begin
        for j in 0..i loop
          y(j) := m(i,j);
        end loop;
        f := Create(x(0..i),y);
        res(i) := new Vector'(Expand(f,x(0..i)));
        Clear(f);
      end;
    end loop;
    return res;
  end Interpolate;

  function Interpolate ( file : file_type; x : Vector; m : Matrix )
                       return VecVec is

  -- DESCRIPTION :
  --   Given the function values needed for the traces in the matrix m,
  --   the expanded coefficient vectors of the trace forms are returned.
  --   Note that the order of indexing is reversed: the i-th index
  --   of the vector on return contains trace d-i.

    d : constant integer32 := m'last(1);
    res : VecVec(1..d);

  begin
    for i in 1..d loop
      declare
        y,f : Vector(0..i);
      begin
        for j in 0..i loop
          y(j) := m(i,j);
        end loop;
        f := Create(x(0..i),y);
        res(i) := new Vector'(Expand(f,x(0..i)));
        put_line(file,"The vector x : "); put_line(file,x(0..i));
        put_line(file,"The vector y : "); put_line(file,y);
        put_line(file,"The coefficients : ");
        put_line(file,res(i).all);
        put_line(file,"evaluating and comparing : ");
        for j in 0..i loop
          put(file,Evalc(res(i).all,x(j))); new_line(file);
          put(file,y(j)); new_line(file);
        end loop;
        Clear(f);
      end;
    end loop;
    return res;
  end Interpolate;

-- AUXILIARIES TO MULTIVARIATE CREATOR :

  function Rotate ( rot : VecVec; x : Vector ) return Vector is

  -- DESCRIPTION :
  --   Given the normals to the slices in the vector of vectors rot
  --   and a point x, the vector on return contains the inner product
  --   of all normals with x.

    res : Vector(rot'range);

  begin
    for i in res'range loop
      res(i) := rot(i)(x'range)*x;
    end loop;
    return res;
  end Rotate;

  function Grid_Points ( grid : Stacked_Sample_Grid ) return VecVec is

  -- DESCRIPTION :
  --   Returns the vectors whose product defines the interpolation grid.

    res : VecVec(1..grid.k);

    procedure Collect ( i : in integer32; sg : in Stacked_Sample_Grid ) is

    -- DESCRIPTION :
    --   Collects the points on the axes of the interpolation grid.

    begin
      res(i) := new Vector'(sg.pts);
      if sg.k > 1
       then Collect(i+1,sg.a(1).all);
      end if;
    end Collect;

  begin
    Collect(1,grid);
    return res;
  end Grid_Points;

  function Samples_on_Grid ( grid : Stacked_Sample_Grid; d : integer32 )
                           return NesVec is

  -- DESCRIPTION :
  --   Returns the (k+1)-st coordinate of every sample point in the grid.

    ind : constant natural32 := natural32(grid.k+1);
    res : NesVec(ind,0,d);

    procedure Extract_Samples
                ( sg : in Stacked_Sample_Grid; v : in out NesVec ) is

      tmp : Multprec_Sample_List;
      spt : Multprec_Sample;

    begin
      if sg.k = 1 then
        for i in v.w'range loop
          v.w(i) := new NesVec(1,1,d);
          tmp := sg.g(i);
          for j in v.w(i).v'range loop
            spt := Head_Of(tmp);
            Copy(Sample_Point(spt).v(integer32(ind)),v.w(i).v(j));
            tmp := Tail_Of(tmp);
          end loop;
        end loop;
      else
        for i in sg.a'range loop
          v.w(i) := new NesVec(v.n-1,v.a,v.b);
          Extract_Samples(sg.a(i).all,v.w(i).all);
        end loop;
      end if;
    end Extract_Samples;

  begin
    Extract_Samples(grid,res);
    return res;
  end Samples_on_Grid;

  function Power_Sum
             ( start,i : integer32; v : Vector ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the i-th power sum of the coordinates in v,
  --   starting in v at the index start.

    res : Complex_Number := Create(integer(0));
    acc : Complex_Number;

  begin
    if i = 1 then
      for k in start..v'last loop
        Add(res,v(k));
      end loop;
    else
      for k in start..v'last-i+1 loop
        acc := v(k)*Power_Sum(k+1,i-1,v);
        Add(res,acc);
        Clear(acc);
      end loop;
    end if;
    return res;
  end Power_Sum;

  function Power_Sums ( y : NesVec; d : integer32 ) return NesVec is

  -- DESCRIPTION :
  --   Given the nested vectors with the (k+1)-th sampled coordinate,
  --   the power sums of degree d are computed.
  --   The matrix on return has one level less than y.

    res : NesVec(y.n-1,0,d);
    minsgn : constant boolean := (d mod 2 = 1);

    procedure Recursive_Compute ( yv : in NesVec; v : in out NesVec ) is
    begin
      if v.n = 1 then
        for i in 0..d loop
          v.v(i) := Power_Sum(1,d,yv.w(i).v);
          if minsgn
           then Min(v.v(i));
          end if;
        end loop;
      else
        for i in 0..d loop
          v.w(i) := new NesVec(v.n-1,0,d);
          Recursive_Compute(yv.w(i).all,v.w(i).all);
        end loop;
      end if;
    end Recursive_Compute;

  begin
    Recursive_Compute(y,res);
    return res;
  end Power_Sums;

-- AUXILIARIES TO EXPAND :

  function Convert ( n : integer32; v1,v2 : Complex_Number ) return Poly is

  -- DESCRIPTION :
  --   Returns the polynomial v2*x1 - v1*x2, represented as a polynomial
  --   in n variables.

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    t.dg(1) := 1;
    t.cf := v2;
    res := Create(t);       -- res is the term v2*x1
    t.dg(1) := 0;
    t.dg(2) := 1;
    t.cf := v1;
    Sub(res,t);             -- res has become v2*x1 - v1*x2
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Convert;

  function Convert ( n : integer32; v : Vector ) return Poly is

  -- DESCRIPTION :
  --   Returns the polynomial v(1)*x1 + v(2)*x2 + .. + v(n)*xn.

    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    for i in 1..n loop
      t.dg(i) := 1;
      t.cf := v(i);
      Add(res,t);
      t.dg(i) := 0;
    end loop;
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Convert;

  function Substitute ( n : integer32; c : Vector; x : Poly ) return Poly is

  -- DESCRIPTION :
  --   Substitutes the linear equation for x into the polynomial
  --   represented by c(0) + c(1)*x + .. + c(d)*x^d, d = c'last,
  --   where x represents x(1)*x1 + x(2)*x2 + .. + x(n)*xn.
  --   The x1,x2,..,xn are the variables of the polynomial on return.
  --   The Horner scheme minimizes the polynomial multiplications.

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    t.cf := c(c'last);
    res := Create(t);
    for i in reverse 0..c'last-1 loop
      Mul(res,x);
      t.cf := c(i);
      Add(res,t);
    end loop;
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Substitute;

  function Leading_Coefficient ( p : Poly ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the leading coefficient of p.

    res : Complex_Number;

    procedure Assign_Coefficient ( t : in Term; continue : out boolean ) is
    begin
      res := t.cf;
      continue := false;
    end Assign_Coefficient;
    procedure Get_Coefficient is new Visiting_Iterator(Assign_Coefficient);

  begin
    Get_Coefficient(p);
    return res;
  end Leading_Coefficient;

  procedure Normalize ( p : in out Poly ) is

  -- DESCRIPTION :
  --   Makes the polynomial p monic by dividing p by the leading coefficient.

    lead : constant Complex_Number := Leading_Coefficient(p);
    one : Complex_Number := Create(integer(1));
    divisor : Complex_Number := one/lead;

  begin
    Mul(p,divisor);
    Clear(one);
    Clear(divisor);
  end Normalize;

  function Add_Last_Variable ( t : Term ) return Term is
 
  -- DESCRIPTION :
  --   Adds a place for a last variable, with power equal to zero.
 
    res : Term;
 
  begin
    Copy(t.cf,res.cf);
    res.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+1);
    res.dg(t.dg'range) := t.dg.all;
    res.dg(res.dg'last) := 0;
    return res;
  end Add_Last_Variable; 

  procedure Add_Last_Variable ( p : in out Poly ) is
 
  -- DESCRIPTION :
  --   Adds to every term space for a last variable, power equal to zero.
 
    res : Poly := Null_Poly;
 
    procedure Add_Term ( t : in Term; continue : out boolean ) is
 
      et : constant Term := Add_Last_Variable(t);
 
    begin
      Add(res,et);
      continue := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);
 
  begin
    Add_Terms(p);
    Clear(p); p := res;
  end Add_Last_Variable;           

  function Substitute ( n : integer32; hyp : VecVec; p : Poly ) return Poly is
 
  -- DESCRIPTION :
  --   Replaces the variables in p by the hyperplanes in hyp.
 
  -- REQUIRED : hyp'range = 1..Number_of_Unknowns(p), n refers to the last
  --   variable in p which will not be substituted.
 
    res : Poly := Null_Poly;
    px : Poly_Sys(hyp'range);
 
    procedure Substitute_Term ( t : in Term; cont : out boolean ) is
 
      et : Poly;
      fac : Term;
 
    begin
      fac.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
      Copy(t.cf,fac.cf);
      fac.dg(n) := t.dg(n);
      et := Create(fac);
      Clear(fac);
      for i in px'range loop
        for j in 1..t.dg(i) loop
          Mul(et,px(i));
        end loop;
      end loop;
      Add(res,et);
      Clear(et);
      cont := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Visiting_Iterator(Substitute_Term);
 
  begin
    for i in hyp'range loop
      px(i) := Convert(n,hyp(i).all);
    end loop;
    Substitute_Terms(p);
    Clear(px);
    return res;
  end Substitute;

-- AUXILIARIES FOR Create_on_Triangle :

  function Roots ( file : file_type;
                   c : Multprec_Complex_Vectors.Vector; size : natural32 )
                 return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Applies Durand-Kerner to compute the roots of the polynomial
  --   defined by the coefficients in c.  The range of c is 0..d,
  --   where d is the degree of the polynomial.

    d : constant integer32 := c'last;
    z : Multprec_Complex_Vectors.Vector(1..d) := Random_Vector(1,d,size);
    maxit : constant natural32 := 10*natural32(d)*size;
    dgts : constant natural32 := 4 - 8*(size+1);
    st_eps : constant double_float := 10.0**integer(dgts);
    mp_eps : Floating_Number := Create(st_eps);
    res : Multprec_Complex_Vectors.Vector(1..d);
    maxres,absres : Floating_Number;
    nb : natural32;
    fail : boolean;

    procedure Multprec_Write
                 ( step : in natural32;
                   z,res : in Multprec_Complex_Vectors.Vector ) is

      absres : Floating_Number;

    begin
      put(file,"Output after step  "); put(file,step,1); put_line(file," :");
      put_line(file,
      "--------------------------------------------------------------------");
      put_line(file,
      "|    APPROXIMATED ROOTS                                            |"); 
      put_line(file,
      "--------------------------------------------------------------------"); 
      for i in z'range loop
         put(file,z(i)); new_line(file);
      end loop;
      put_line(file,
      "--------------------------------------------------------------------");
      put_line(file,
      "|    RESIDUALS                                                     |");
      put_line(file,
      "--------------------------------------------------------------------");
      for i in z'range loop
        absres := AbsVal(res(i));
        put(file,absres); new_line(file);
        Clear(absres);
      end loop;
      put_line(file,
      "--------------------------------------------------------------------");
    end Multprec_Write;

    procedure hbdk is
      new Hybrid_Durand_Kerner.Reporting_Durand_Kerner(Multprec_Write);

  begin
   -- Silent_Durand_Kerner(c,z,res,maxit,mp_eps,nb,fail);
    hbdk(c,z,res,maxit,mp_eps,nb,fail);
    put_line(file,"The roots of the polynomial : "); put_line(file,z);
    put_line(file,"The residuals at the roots : "); put_line(file,res);
    maxres := AbsVal(res(1));
    for i in 2..d loop
      absres := Absval(res(i));
      if absres > maxres
       then Copy(absres,maxres);
      end if;
      Clear(absres);
    end loop;
    put(file,"Maximal residual at the roots : ");
    put(file,maxres,3); new_line(file);
    put(file,"Number of iterations : "); put(file,nb); new_line(file);
    Clear(mp_eps); Clear(res); Clear(maxres);
    return z;
  end Roots;

  function Roots ( c : Multprec_Complex_Vectors.Vector; size : natural32 )
                 return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Applies Durand-Kerner to compute the roots of the polynomial
  --   defined by the coefficients in c.  The range of c is 0..d,
  --   where d is the degree of the polynomial.

    d : constant integer32 := c'last;
    z : Multprec_Complex_Vectors.Vector(1..d) := Random_Vector(1,d,size);
    maxit : constant natural32 := 10*natural32(d)*size;
    dgts : constant natural32 := 4 - 8*(size+1);
    st_eps : constant double_float := 10.0**integer(dgts);
    mp_eps : Floating_Number := Create(st_eps);
    res : Multprec_Complex_Vectors.Vector(1..d);
    nb : natural32;
    fail : boolean;

  begin
    Silent_Durand_Kerner(c,z,res,maxit,mp_eps,nb,fail);
    Clear(mp_eps); Clear(res);
    return z;
  end Roots;

  procedure Newton_Samples 
                 ( file : in file_type; ya : in out Matrix;
                   x : in Vector; i : in integer32; size : in natural32;
                   trc : in VecVec ) is

  -- DESCRIPTION :
  --   Applies the Newton identities to add extra samples to the
  --   triangular grid of rotated and projected samples.

  -- ON ENTRY :
  --   file        for intermediate output and diagnostics;
  --   ya          matrix of projected sample points from triangular grid:
  --               in the i-th column we have values in rows 1..d-i+1;
  --   x           vector of abscisses in the grid;
  --   i           current column in ya to add new samples to;
  --   size        size of the multi-precision numbers;
  --   trc         coefficients of traces trc(j), for j = 1..i-1.

  -- ON RETURN :
  --   ya          updated matrix with projected sample points.

    d : constant integer32 := ya'last(1);

  begin
    if i = 2 then
      ya(d,2) := Evalc(trc(1).all,x(2));  -- traces come with sign
      Min(ya(d,2));
      for k in ya'first(1)..d-1 loop
        Sub(ya(d,2),ya(k,2));
      end loop;
    else
      declare
        valtrc : Vector(1..i-1);  -- values at j-th trace, j=1,..,i-1
        valpow : Vector(1..i-1);  -- values of power sums
        polcff : Vector(0..i-1);  -- coefficients of polynomial
        ind : integer32;
        pow : Complex_Number;
      begin
        for j in valtrc'range loop
          valtrc(j) := Evalc(trc(j).all,x(i));
        end loop;
        valpow := Traces_to_Power_Sums(valtrc);
        for j in valpow'range loop  -- modify the power sums
          for k in 1..d-i+1 loop    -- use values ya(k,i)
            pow := ya(k,i)**integer(j);
            Sub(valpow(j),pow);
            Clear(pow);
          end loop;
        end loop;
        Clear(valtrc);
        valtrc := Power_Sums_to_Traces(valpow);   -- modified traces
        polcff(i-1) := Create(integer(1));                 -- monic polynomial
        for j in valtrc'range loop    -- fill in coefficients of poly
          Copy(valtrc(j),polcff(i-1-j));
        end loop;
        Clear(valpow);
        valpow := Roots(file,polcff,size);        -- missing values
        ind := 1;
        for k in d-i+2..d loop                    -- fill in i-1 values
          Copy(valpow(ind),ya(k,i));
          ind := ind+1;
        end loop;
        Clear(valtrc); Clear(valpow); Clear(polcff);
      end;
    end if;
   -- put_line(file,"The grid of projected samples with the added ones :");
   -- put(file,ya,3);
  end Newton_Samples;

  procedure Newton_Samples 
                 ( ya : in out Matrix;
                   x : in Vector; i : in integer32; size : in natural32;
                   trc : in VecVec ) is

  -- DESCRIPTION :
  --   Applies the Newton identities to add extra samples to the
  --   triangular grid of rotated and projected samples.

  -- ON ENTRY :
  --   ya          matrix of projected sample points from triangular grid:
  --               in the i-th column we have values in rows 1..d-i+1;
  --   x           vector of abscisses in the grid;
  --   i           current column in ya to add new samples to;
  --   size        size of the multi-precision numbers;
  --   trc         coefficients of traces trc(j), for j = 1..i-1.

  -- ON RETURN :
  --   ya          updated matrix with projected sample points.

    d : constant integer32 := ya'last(1);

  begin
    if i = 2 then
      ya(d,2) := Evalc(trc(1).all,x(2));  -- traces come with sign
      Min(ya(d,2));
      for k in ya'first(1)..d-1 loop
        Sub(ya(d,2),ya(k,2));
      end loop;
    else
      declare
        valtrc : Vector(1..i-1);  -- values at j-th trace, j=1,..,i-1
        valpow : Vector(1..i-1);  -- values of power sums
        polcff : Vector(0..i-1);  -- coefficients of polynomial
        ind : integer32;
        pow : Complex_Number;
      begin
        for j in valtrc'range loop
          valtrc(j) := Evalc(trc(j).all,x(i));
        end loop;
        valpow := Traces_to_Power_Sums(valtrc);
        for j in valpow'range loop  -- modify the power sums
          for k in 1..d-i+1 loop    -- use values ya(k,i)
            pow := ya(k,i)**integer(j);
            Sub(valpow(j),pow);
            Clear(pow);
          end loop;
        end loop;
        Clear(valtrc);
        valtrc := Power_Sums_to_Traces(valpow);   -- modified traces
        polcff(i-1) := Create(integer(1));                 -- monic polynomial
        for j in valtrc'range loop    -- fill in coefficients of poly
          Copy(valtrc(j),polcff(i-1-j));
        end loop;
        Clear(valpow);
        valpow := Roots(polcff,size);             -- missing values
        ind := 1;
        for k in d-i+2..d loop                    -- fill in i-1 values
          Copy(valpow(ind),ya(k,i));
          ind := ind+1;
        end loop;
        Clear(valtrc); Clear(valpow); Clear(polcff);
      end;
    end if;
  end Newton_Samples;

-- CREATORS :

  function Create ( grid : Array_of_Multprec_Sample_Lists )
                  return Trace_Interpolator1_Rep is

    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := grid'last;
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    res : Trace_Interpolator1_Rep(n,d);
    ya,ps : Matrix(1..d,0..d);
    x : Vector(0..d) := Abscisses(grid,natural32(d));

  begin
    res.rot := hyp(hyp'first)(1..n);
    ya := Rotate_Samples(natural32(d),natural32(d),res.rot,grid);
    ps := Power_Sums(ya);
    res.trc := Interpolate(x,ps);
    Clear(ya); Clear(ps); Clear(x);
    return res;
  end Create;

  function Create ( grid : Array_of_Multprec_Sample_Lists )
                  return Trace_Interpolator1 is

    res : Trace_Interpolator1;

  begin
    res := new Trace_Interpolator1_Rep'(Create(grid));
    return res;
  end Create;

  function Create ( file : file_type; grid : Array_of_Multprec_Sample_Lists )
                  return Trace_Interpolator1_Rep is

    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := grid'last;
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    res : Trace_Interpolator1_Rep(n,d);
    ya,ps : Matrix(1..d,0..d);
    x : Vector(0..d) := Abscisses(grid,natural32(d));

  begin
    res.rot := hyp(hyp'first)(1..n);
    ya := Rotate_Samples(natural32(d),natural32(d),res.rot,grid);
    ps := Power_Sums(file,ya);
    res.trc := Interpolate(file,x,ps);
    Clear(ya); Clear(ps); Clear(x);
    return res;
  end Create;

  function Create ( file : file_type; grid : Array_of_Multprec_Sample_Lists )
                  return Trace_Interpolator1 is

    res : Trace_Interpolator1;

  begin
    res := new Trace_Interpolator1_Rep'(Create(file,grid));
    return res;
  end Create;

  function Create ( grid : Array_of_Multprec_Sample_Lists; i : integer32 )
                  return Vector is

    res : Vector(0..i);
    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := integer32(Length_Of(grid(0)));
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    rot : constant Vector := hyp(hyp'first)(1..n);
    x : Vector(0..i) := Abscisses(grid,natural32(i));
    ya : Matrix(1..d,0..i);
    ps,f : Vector(0..i);

  begin
    ya := Rotate_Samples(natural32(d),natural32(i),rot,grid);
    ps := Power_Sums(i,ya);
    f := Create(x,ps);
    res := Expand(f,x);
    Clear(ya); Clear(ps); Clear(f); Clear(x);
    return res;
  end Create;

  function Create_on_Triangle
                  ( file : file_type; grid : Array_of_Multprec_Sample_Lists;
                    size : natural32 )
                  return Trace_Interpolator1 is

    res : Trace_Interpolator1;
    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := integer32(Length_Of(grid(0)));
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    res_rep : Trace_Interpolator1_Rep(n,d);
    x : constant Vector(0..d) := Abscisses(grid,natural32(d));
    ya : Matrix(1..d,0..d);
    ps,f : Vector(0..d);

  begin
    res_rep.rot := hyp(hyp'first)(1..n);
    ya := Rotate_Samples(natural32(d),natural32(d),res_rep.rot,grid);
   -- put_line(file,"The rotated samples :"); put(file,ya,3);
    ps(0..1) := Power_Sums(1,ya);
    f(0..1) := Create(x(0..1),ps(0..1));
    res_rep.trc(1) := new Vector'(Expand(f(0..1),x(0..1)));
    Clear(ps(0..1)); Clear(f(0..1));
    for i in 2..d loop
       Newton_Samples(file,ya,x,i,size,res_rep.trc);
       ps(0..i) := Power_Sums(i,ya);
       f(0..i) := Create(x(0..i),ps(0..i));
       res_rep.trc(i) := new Vector'(Expand(f(0..i),x(0..i)));
       Clear(ps(0..i)); Clear(f(0..i));
    end loop;
    res := new Trace_Interpolator1_Rep'(res_rep);
    return res;
  end Create_on_Triangle;

  function Create_on_Triangle
                  ( grid : Array_of_Multprec_Sample_Lists; size : natural32 )
                  return Trace_Interpolator1 is

    res : Trace_Interpolator1;
    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := integer32(Length_Of(grid(0)));
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    res_rep : Trace_Interpolator1_Rep(n,d);
    x : constant Vector(0..d) := Abscisses(grid,natural32(d));
    ya : Matrix(1..d,0..d);
    ps,f : Vector(0..d);

  begin
    res_rep.rot := hyp(hyp'first)(1..n);
    ya := Rotate_Samples(natural32(d),natural32(d),res_rep.rot,grid);
    ps(0..1) := Power_Sums(1,ya);
    f(0..1) := Create(x(0..1),ps(0..1));
    res_rep.trc(1) := new Vector'(Expand(f(0..1),x(0..1)));
    Clear(ps(0..1)); Clear(f(0..1));
    for i in 2..d loop
       Newton_Samples(ya,x,i,size,res_rep.trc);
       ps(0..i) := Power_Sums(i,ya);
       f(0..i) := Create(x(0..i),ps(0..i));
       res_rep.trc(i) := new Vector'(Expand(f(0..i),x(0..i)));
       Clear(ps(0..i)); Clear(f(0..i));
    end loop;
    res := new Trace_Interpolator1_Rep'(res_rep);
    return res;
  end Create_on_Triangle;

  function Create ( grid : Stacked_Sample_Grid; d : integer32 )
                  return Trace_Interpolator_Rep is

  -- ASSUMPTION : grid.k > 1.

    res : Trace_Interpolator_Rep(grid.k,integer32(grid.n),d);
    x : constant VecVec(1..grid.k) := Grid_Points(grid);
    y : constant NesVec(natural32(grid.k)+1,0,d) := Samples_on_Grid(grid,d);

  begin
    res.rot := grid.hyp;
    for i in 1..d loop
      declare
        ps : constant NesVec(natural32(grid.k),0,i) := Power_Sums(y,i);
        nf : Newton_Form(natural32(grid.k),grid.k,i);
      begin
        nf.x := x;
        nf.f := Create(natural32(grid.k),natural32(i),x,ps);
        res.trc(i) := new Newton_Form'(nf);
      end;
    end loop;
    return res;
  end Create;

  function Create ( file : file_type;
                    grid : Stacked_Sample_Grid; d : integer32 )
                  return Trace_Interpolator_Rep is

  -- ASSUMPTION : grid.k > 1.

    res : Trace_Interpolator_Rep(grid.k,integer32(grid.n),d);
    x : constant VecVec(1..grid.k) := Grid_Points(grid);
    y : constant NesVec(natural32(grid.k)+1,0,d) := Samples_on_Grid(grid,d);

  begin
    put_line(file,"The vectors whose product defines the grid : ");
    put_line(file,x);
    put_line(file,"The multi-dimensional matrix of extracted samples : ");
    put(file,y);
    res.rot := grid.hyp;
    for i in 1..d loop
      put(file,"Compute trace of degree "); put(file,i,1);
      put_line(file," :");
      declare
        ps : constant NesVec(natural32(grid.k),0,i) := Power_Sums(y,i);
        nf : Newton_Form(natural32(grid.k),grid.k,i);
        maxerr : Floating_Number;
      begin
        put_line(file,"The power sums :"); put(file,ps);
        nf.x := x;
        nf.f := Create_on_Square(natural32(grid.k),natural32(i),x,ps);
        put_line(file,"The generalized divided differences : ");
        put(file,nf.f);
        Eval_Grid(file,nf.x,ps,nf.f,maxerr);
        res.trc(i) := new Newton_Form'(nf);
      end;
    end loop;
    return res;
  end Create;

  function Create ( grid : Stacked_Sample_Grid; d : integer32 )
                  return Trace_Interpolator is

    res : Trace_Interpolator;
    res_rep : Trace_Interpolator_Rep(grid.k,integer32(grid.n),d);

  begin
    if grid.k = 1
     then res_rep.t1 := Create(grid.g.all);
     else res_rep := Create(grid,d);
    end if;
    res := new Trace_Interpolator_Rep'(res_rep);
    return res;
  end Create;

  function Create ( file : file_type;
                    grid : Stacked_Sample_Grid; d : integer32 )
                  return Trace_Interpolator is

    res : Trace_Interpolator;
    res_rep : Trace_Interpolator_Rep(grid.k,integer32(grid.n),d);

  begin
    if grid.k = 1
     then res_rep.t1 := Create(grid.g.all);
     else res_rep := Create(file,grid,d);
    end if;
    res := new Trace_Interpolator_Rep'(res_rep);
    return res;
  end Create;

  function Expand ( t : Trace_Interpolator1_Rep ) return Poly is

    res : Poly := Null_Poly;
    n : constant integer32 := t.n-1;
    px : Poly := Convert(n,t.rot);
    py : Poly := Convert(n,t.rot(1),t.rot(2));
    trace : Poly;
    tt : Term;

  begin
    tt.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    tt.cf := Create(integer(1));
    res := Create(tt);
    for i in 1..t.d loop
      Mul(res,py);
      trace := Substitute(n,t.trc(i).all,px);
      Add(res,trace);
      Clear(trace);
    end loop;
    Normalize(res);
    Clear(px); Clear(py); Clear(tt);                
    return res;
  end Expand;

  function Expand ( t : Trace_Interpolator_Rep ) return Poly is

    res,tp,subres : Poly;
    fac : Term;
    nv : integer32;
 
  begin
    if t.k = 1 then
      res := Expand(t.t1);
    else
      res := Null_Poly;
      fac.dg := new Standard_Natural_Vectors.Vector'(1..t.k+1 => 0);
      fac.dg(t.k+1) := natural32(t.d);
      fac.cf := Create(integer(1));
      res := Create(fac);
      for i in 1..t.d loop
        tp := Expand(t.trc(i).f,t.trc(i).x);
        Add_Last_Variable(tp);
        fac.dg(t.k+1) := natural32(t.d-i);
        Mul(tp,fac);
        Add(res,tp);
        Clear(tp);
      end loop;
      Clear(fac);                              
      nv := t.k+1; --t.n - t.k;
      subres := Substitute(nv,t.rot,res);
      Clear(res); res := subres;
      Normalize(res);
    end if;
    return res;
  end Expand;

  function Expand ( t : Trace_Interpolator1 ) return Poly is
  begin
    if t = null
     then return Null_Poly;
     else return Expand(t.all);
    end if;
  end Expand;
 
  function Expand ( t : Trace_Interpolator ) return Poly is
  begin
    if t = null
     then return Null_Poly;
     else return Expand(t.all);
    end if;
  end Expand;

-- SELECTORS :

  function Degree ( t : Trace_Interpolator1 ) return integer32 is
  begin
    if t = null
     then return -1;
     else return t.d;
    end if;
  end Degree;

  function Dimension ( t : Trace_Interpolator1 ) return integer32 is
  begin
    if t = null
     then return -1;
     else return t.n;
    end if;
  end Dimension;

  function Trace ( t : Trace_Interpolator1; i : integer32 ) return Vector is
  begin
    return t.trc(i).all;
  end Trace;

-- EVALUATORS :

  function Eval ( t : Trace_Interpolator1_Rep; x : Vector )
                return Complex_Number is

    res,acc : Complex_Number;
    pt : constant Vector(1..2) := Rotate_and_Project(t.rot,x);

  begin
    res := Create(integer(1));                   -- monic
    for i in 1..t.d loop
      Mul(res,pt(2));
      acc := Evalc(t.trc(i).all,pt(1));
      Add(res,acc);
      Clear(acc);
    end loop;
    return res;
  end Eval;

  function Eval ( t : Trace_Interpolator1; x : Vector )
                return Complex_Number is

  begin
    if t = null
     then return Create(integer(0));
     else return Eval(t.all,x);
    end if;
  end Eval;

  procedure Eval_Trace ( t : in Vector; d : in integer32;
                         sps : in Multprec_Sample_List;
                         val,eva : out Complex_Number ) is

    hyp : constant VecVec := Hyperplane_Sections(Head_Of(sps));
    tmp : Multprec_Sample_List := sps;
    len : constant integer32 := integer32(Length_Of(sps));
    y : Matrix(1..len,1..1);
    x : Complex_Number; -- := -hyp(1)(0);
    n : constant integer32 := Number_of_Variables(Head_Of(sps));

  begin
    x := hyp(1)(1..n)*Sample_Point(Head_Of(tmp)).v;
    val := Evalc(t,x);
    for i in 1..len loop
      y(i,1) := Rotate_and_Project2(hyp(1).all,Sample_Point(Head_Of(tmp)).v);
      tmp := Tail_Of(tmp);
    end loop;
    eva := Power_Sum(1,d,1,y);
    if d mod 2 = 1
     then Min(eva);
    end if;
    Clear(y);
  end Eval_Trace;

  function Eval ( t : Trace_Interpolator_Rep; x : Vector )
                return Complex_Number is

    res : Complex_Number;

  begin
    if t.k = 1 then
      res := Eval(t.t1,x);
    else
      res := Create(integer(1));
      declare
        pt : Vector(t.rot'range) := Rotate(t.rot,x);
        eva : Complex_Number;
      begin
        for i in 1..t.d loop
          eva := Eval0(t.trc(i).f,t.trc(i).x,pt);
          Mul(res,x(t.k+1));
          Add(res,eva);
          Clear(eva);
        end loop;
        Clear(pt);
      end;
    end if;
    return res;
  end Eval;

  function Eval ( file : file_type; t : Trace_Interpolator_Rep;
                  x : Vector ) return Complex_Number is

    res : Complex_Number;

  begin
    if t.k = 1 then
      res := Eval(t.t1,x);
    else
      res := Create(integer(1));
      declare
        pt : constant Vector(t.rot'range) := Rotate(t.rot,x);
        eva : Complex_Number;
      begin
        put_line(file,"Eval at the rotated point : ");
        put_line(file,pt);
        put(file,"x("); put(file,t.k+1,1); put(file,") : ");
        put(file,x(t.k+1)); new_line(file);
        for i in 1..t.d loop
          eva := Eval0(t.trc(i).f,t.trc(i).x,pt);
          put(file,"Trc "); put(file,i,1); put(file," Val : ");
          put(file,eva); new_line(file);
          Mul(res,x(t.k+1));
          Add(res,eva);
          Clear(eva);
        end loop;
      end;
    end if;
    return res;
  end Eval;

  function Eval ( t : Trace_Interpolator; x : Vector ) return Complex_Number is
  begin
    if t = null
     then return Create(integer(0));
     else return Eval(t.all,x);
    end if;
  end Eval;

  function Eval ( file : file_type;
                  t : Trace_Interpolator; x : Vector ) return Complex_Number is
  begin
    if t = null
     then return Create(integer(0));
     else return Eval(t.all,x);
    end if;
  end Eval;

-- DIAGNOSTICS :

  function Errors ( t : Trace_Interpolator1;
                    grid : Array_of_Multprec_Sample_Lists )
                  return Multprec_Floating_Matrices.Matrix is

    res : Multprec_Floating_Matrices.Matrix
            (grid'range,1..integer32(Length_Of(grid(grid'first))));
    tmp : Multprec_Sample_List;

  begin
    for i in grid'range loop
      tmp := grid(i);
      for j in res'range(2) loop
        Clear(res(i,j));
        res(i,j) := AbsVal(Eval(t,Sample_Point(Head_Of(tmp)).v));
        tmp := Tail_Of(tmp);
        if Is_Null(tmp) then
          for k in j+1..res'last(2) loop
            res(i,k) := Create(integer(0));
          end loop;
          exit;
        end if;
      end loop;
    end loop;
    return res;
  end Errors;

  procedure Write_Errors ( file : in file_type; t : in Trace_Interpolator1;
                           grid : in Array_of_Multprec_Sample_Lists ) is

    tmp : Multprec_Sample_List;

  begin
    for i in grid'range loop
      tmp := grid(i);
      for j in 1..integer32(Length_Of(grid(i))) loop
        put(file,"("); put(file,i,1); put(file,",");
        put(file,j,1); put(file,") : ");
        put(file,Eval(t,Sample_Point(Head_Of(tmp)).v));
        new_line(file);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Write_Errors;

  function Maximal_Error ( residuals : Multprec_Floating_Matrices.Matrix )
                         return Floating_Number is

    max : Floating_Number;

  begin
    Copy(residuals(residuals'first(1),residuals'first(2)),max);
    for i in residuals'range(1) loop
      for j in residuals'range(2) loop
        if residuals(i,j) > max
         then Copy(residuals(i,j),max);
        end if;
      end loop;
    end loop;
    return max;
  end Maximal_Error;

  function Maximal_Error ( t : Trace_Interpolator1;
                           grid : Array_of_Multprec_Sample_Lists )
                         return Floating_Number is
  begin
    return Maximal_Error(Errors(t,grid));
  end Maximal_Error;

  procedure Write_Errors
              ( file : in file_type; t : in Trace_Interpolator;
                grid : in Stacked_Sample_Grid;
                maxerr : out Floating_Number ) is

    res : Floating_Number := Create(integer32(-1));
    tmp : Multprec_Sample_List;
    eva : Complex_Number;
    abseva : Floating_Number;

  begin
    put(file,"Evaluating trace interpolation in stacked grid at level ");
    put(file,grid.k,1); put_line(file," :");
    if grid.k = 1 then
      for i in grid.g'range loop
        tmp := grid.g(i);
        while not Is_Null(tmp) loop
          eva := Eval(file,t,Sample_Point(Head_Of(tmp)).v);
          put(file,eva); new_line(file);
          abseva := AbsVal(eva);
          if abseva > res
           then Copy(abseva,res);
          end if;
          Clear(abseva); Clear(eva);
          tmp := Tail_Of(tmp);
        end loop;
      end loop;
    else
      for i in grid.a'range loop
        Write_Errors(file,t,grid.a(i).all,abseva);
        if abseva > res
         then Copy(abseva,res);
        end if;
        Clear(abseva);
      end loop;
    end if;
    maxerr := res;
  end Write_Errors;

  function Maximal_Error
             ( t : Trace_Interpolator; grid : Stacked_Sample_Grid )
             return Floating_Number is

    res : Floating_Number := Create(integer(-1));
    tmp : Multprec_Sample_List;
    eva : Complex_Number;
    abseva : Floating_Number;

  begin
    if grid.k = 1 then
      for i in grid.g'range loop
        tmp := grid.g(i);
        while not Is_Null(tmp) loop
          eva := Eval(t,Sample_Point(Head_Of(tmp)).v);
          abseva := AbsVal(eva);
          if abseva > res
           then Copy(abseva,res);
          end if;
          Clear(abseva); Clear(eva);
          tmp := Tail_Of(tmp);
        end loop;
      end loop;
    else
      for i in grid.a'range loop
        abseva := Maximal_Error(t,grid.a(i).all);
        if abseva > res
         then Copy(abseva,res);
        end if;
        Clear(abseva);
      end loop;
    end if;
    return res;
  end Maximal_Error;

-- DESTRUCTORS :

  procedure Clear ( t : in out Trace_Interpolator1 ) is

    procedure free is new unchecked_deallocation
      (Trace_Interpolator1_Rep,Trace_Interpolator1);

  begin
    if t /= null
     then Multprec_Complex_VecVecs.Clear(t.trc);
          free(t);
    end if;
  end Clear;

  procedure Clear ( nf : in out Newton_Form ) is
  begin
    Clear(nf.x);
    Clear(nf.f);
  end Clear;

  procedure Clear ( nf : in out Link_to_Newton_Form ) is

    procedure free is new unchecked_deallocation
      (Newton_Form,Link_to_Newton_Form);

  begin
    if nf /= null
     then Clear(nf.all);
          free(nf);
    end if;
  end Clear;

  procedure Clear ( t : in out Trace_Interpolator ) is

    procedure free is new unchecked_deallocation
      (Trace_Interpolator_Rep,Trace_Interpolator);

  begin
    if t /= null then
      if t.k = 1 then
        Multprec_Complex_VecVecs.Clear(t.t1.trc);
      else
        Clear(t.rot);
        for i in t.trc'range loop
          Clear(t.trc(i));
        end loop;
      end if;
      free(t);
    end if;
  end Clear;

end Multprec_Trace_Interpolators;
