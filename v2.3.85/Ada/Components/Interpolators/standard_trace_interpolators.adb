with unchecked_deallocation;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
--with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Natural_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_NesVecs;           use Standard_Complex_NesVecs;
with Standard_Complex_NesVecs_io;        use Standard_Complex_NesVecs_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Sample_Points;                      use Sample_Points;
with Rectangular_Sample_Grids;           use Rectangular_Sample_Grids;
with Standard_Univariate_Interpolators;  use Standard_Univariate_Interpolators;
with Standard_Power_Traces;              use Standard_Power_Traces;
with Standard_Durand_Kerner;             use Standard_Durand_Kerner;

package body Standard_Trace_Interpolators is

-- DATA STRUCTURES :

  type Trace_Interpolator1_Rep ( n,d : integer32 ) is record
    rot : Vector(1..n);   -- normal to all parallel hyperplanes in grid
    trc : VecVec(1..d);   -- trc(i) contains coefficients of trace of deg i
  end record;

  type Trace_Interpolator_Rep ( k,n,d : integer32 ) is record
    case k is  -- k = #slices; n = dimension of ambient space; d = degree
      when 1 => t1 : Trace_Interpolator1_Rep(n,d);
      when others => rot : VecVec(1..k);        -- normals to slices
                     trc : Newton_Forms(1..d);  -- low dim trace forms
    end case;
  end record;

-- AUXILIARIES TO THE CREATOR :

  function Power_Sum ( start,i,j : integer32; m : Matrix )
                     return Complex_Number is

  -- DESCRIPTION :
  --   Returns the i-th power sum at the j-th slice, using the rotated
  --   and projected 2nd coordinate of the roots, starting with roots
  --   at the index start.

    res : Complex_Number := Create(integer(0));

  begin
    if i = 1 then
      for k in start..m'last(1) loop
        res := res + m(k,j);
      end loop;
    else
      for k in start..m'last(1)-i+1 loop
        res := res + m(k,j)*Power_Sum(k+1,i-1,j,m);
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

  begin
    if i = 1 then
      for k in start..m'last(1) loop
        put(file," + m("); put(file,k,1); put(file,",");
        put(file,j,1); put(file,")");
        res := res + m(k,j);
      end loop;
      new_line(file);
    else
      for k in start..m'last(1)-i+1 loop
        put(file," + m("); put(file,k,1); put(file,",");
        put(file,j,1); put(file,")*");
        res := res + m(k,j)*Power_Sum(k+1,i-1,j,m);
      end loop;
      new_line(file);
    end if;
    return res;
  end Power_Sum;

  procedure Adjust_Signs ( m : in out Matrix ) is

  -- DESCRIPTION :
  --   Multiplies the i-th row of the matrix with (-1)^n.

    min : boolean := (m'last(1) mod 2 = 1);

  begin
    for i in reverse m'range(1) loop
      if min
       then for j in m'range(2) loop
              m(i,j) := -m(i,j);
            end loop;
      end if;
      min := not min;
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
      end;
    end loop;
    return res;
  end Interpolate;

  function Interpolate ( file : file_type; x : Vector; m : Matrix )
                       return VecVec is

  -- DESCRIPTION :
  --   Given the function values needed for the traces in the matrix m,
  --   the expanded coefficient vectors of the trace forms are returned.

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
      end;
    end loop;
    return res;
  end Interpolate;

-- AUXILIARIES TO MULTIVARIATE CREATOR :

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
  --   Returns a random projection of the samples in the grid.
  --   Returns the (k+1)-st coordinate of every sample point in the grid,
  --   set up to interpolate the i-th trace of a degree d component.

    ind : constant integer32 := grid.k+1;
    res : NesVec(natural32(ind),0,d);
   -- provec : Vector(1..grid.n) := Random_Vector(1,grid.n);

    procedure Extract_Samples 
                ( sg : in Stacked_Sample_Grid; v : in out NesVec ) is

      tmp : Standard_Sample_List;
      spt : Standard_Sample;

    begin
      if sg.k = 1 then
        for j1 in v.w'range loop
          v.w(j1) := new NesVec(1,1,d);
          tmp := sg.g(j1);
          for j2 in v.w(j1).v'range loop
            spt := Head_Of(tmp);
            v.w(j1).v(j2) := Sample_Point(spt).v(ind);
           -- declare
           --   pt : constant Vector := Sample_Point(spt).v;
           -- begin
           --   v.w(j1).v(j2) := pt(provec'range)*provec;
           -- end;
            tmp := Tail_Of(tmp);
          end loop;
        end loop;
      else
        for j1 in sg.a'range loop
          v.w(j1) := new NesVec(v.n-1,v.a,v.b);
          Extract_Samples(sg.a(j1).all,v.w(j1).all);
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

  begin
    if i = 1 then
      for k in start..v'last loop
        res := res + v(k);
      end loop;
    else
      for k in start..v'last-i+1 loop
        res := res + v(k)*Power_Sum(k+1,i-1,v);
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
    min : constant boolean := (d mod 2 = 1);

    procedure Recursive_Compute ( yv : in NesVec; v : in out NesVec ) is
    begin
      if v.n = 1 then
        for i in 0..d loop
          v.v(i) := Power_Sum(1,d,yv.w(i).v);
          if min
           then v.v(i) := -v.v(i);
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
    Clear(t);
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
    Clear(t);
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
    Clear(t);
    return res;
  end Substitute;

  function Leading_Coefficient ( p : Poly; tol : double_float )
                               return Complex_Number is

  -- DESCRIPTION :
  --   Returns the leading coefficient of p, under the condition
  --   that it be higher than tol.

    res : Complex_Number;

    procedure Assign_Coefficient ( t : in Term; continue : out boolean ) is
    begin
      if AbsVal(t.cf) > tol
       then res := t.cf;
            continue := false;
       else continue := true;
      end if;
    end Assign_Coefficient;
    procedure Get_Coefficient is new Visiting_Iterator(Assign_Coefficient);

  begin
    Get_Coefficient(p);
    return res;
  end Leading_Coefficient;

  procedure Normalize ( p : in out Poly ) is

  -- DESCRIPTION :
  --   Makes the polynomial p monic by dividing p by the leading coefficient.

    tol : constant double_float := 1.0E-10;
    lead : constant Complex_Number := Leading_Coefficient(p,tol);
    divisor : constant Complex_Number := Create(1.0)/lead;

  begin
    Mul(p,divisor);
  end Normalize;

  function Add_Last_Variable ( t : Term ) return Term is

  -- DESCRIPTION :
  --   Adds a place for a last variable, with power equal to zero.

    res : Term;

  begin
    res.cf := t.cf;
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
      fac.cf := t.cf;
      fac.dg(n) := t.dg(n);
      et := Create(fac);
      Clear(fac.dg);
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
   -- put_line("The polynomial before substitution : ");
   -- put_line(p);
    for i in hyp'range loop
      px(i) := Convert(n,hyp(i).all);
     -- put("px("); put(i,1); put_line(") : ");
     -- put_line(px(i));
    end loop;
    Substitute_Terms(p);
    Clear(px);
    return res;
  end Substitute;

-- AUXILIARIES FOR Create_on_Triangle :

  function Roots ( file : file_type;
                   c : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Applies Durand-Kerner to compute the roots of the polynomial
  --   defined by the coefficients in c.  The range of c is 0..d,
  --   where d is the degree of the polynomial.

    d : constant integer32 := c'last;
    z : Standard_Complex_Vectors.Vector(1..d) := Random_Vector(1,d) ;
    maxit : constant natural32 := 10*natural32(d);
    eps : constant double_float := 1.0E-14;
    res : Standard_Complex_Vectors.Vector(1..d);
    maxres : double_float;
    nb : natural32;
    fail : boolean;

  begin
    Silent_Durand_Kerner(c,z,res,maxit,eps,nb,fail);
    put_line(file,"The roots of the polynomial : "); put_line(file,z);
    put_line(file,"The residuals at the roots : "); put_line(file,res);
    maxres := AbsVal(res(1));
    for i in 2..d loop
      if AbsVal(res(i)) > maxres
       then maxres := AbsVal(res(i));
      end if;
    end loop;
    put(file,"Maximal residual at the roots : ");
    put(file,maxres,3); new_line(file);
    return z;
  end Roots;

  function Roots ( c : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Applies Durand-Kerner to compute the roots of the polynomial
  --   defined by the coefficients in c.  The range of c is 0..d,
  --   where d is the degree of the polynomial.

    d : constant integer32 := c'last;
    z : Standard_Complex_Vectors.Vector(1..d) := Random_Vector(1,d) ;
    maxit : constant natural32 := 10*natural32(d);
    eps : constant double_float := 1.0E-14;
    res : Standard_Complex_Vectors.Vector(1..d);
    nb : natural32;
    fail : boolean;

  begin
    Silent_Durand_Kerner(c,z,res,maxit,eps,nb,fail);
    return z;
  end Roots;

  procedure Newton_Samples 
                 ( file : in file_type; ya : in out Matrix;
                   x : in Vector; i : in integer32; trc : in VecVec ) is

  -- DESCRIPTION :
  --   Applies the Newton identities to add extra samples to the
  --   triangular grid of rotated and projected samples.

  -- ON ENTRY :
  --   file        for intermediate output and diagnostics;
  --   ya          matrix of projected sample points from triangular grid:
  --               in the i-th column we have values in rows 1..d-i+1;
  --   x           vector of abscisses in the grid;
  --   i           current column in ya to add new samples to;
  --   trc         coefficients of traces trc(j), for j = 1..i-1.

  -- ON RETURN :
  --   ya          updated matrix with projected sample points.

    d : constant integer32 := ya'last(1);

  begin
    if i = 2 then
      ya(d,2) := -Evalc(trc(1).all,x(2));  -- traces come with sign
      for k in ya'first(1)..d-1 loop
        ya(d,2) := ya(d,2) - ya(k,2);
      end loop;
    else
      declare
        valtrc : Vector(1..i-1);  -- values at j-th trace, j=1,..,i-1
        valpow : Vector(1..i-1);  -- values of power sums
        polcff : Vector(0..i-1);  -- coefficients of polynomial
        ind : integer32;
      begin
        for j in valtrc'range loop
          valtrc(j) := Evalc(trc(j).all,x(i));
        end loop;
        valpow := Traces_to_Power_Sums(valtrc);
        for j in valpow'range loop  -- modify the power sums
          for k in 1..d-i+1 loop    -- use values ya(k,i)
            valpow(j) := valpow(j) - ya(k,i)**integer(j);
          end loop;
        end loop;
        valtrc := Power_Sums_to_Traces(valpow);   -- modified traces
        polcff(i-1) := Create(1.0);               -- monic polynomial
        for j in valtrc'range loop    -- fill in coefficients of poly
          polcff(i-1-j) := valtrc(j);
        end loop;
        valpow := Roots(file,polcff);             -- missing values
        ind := 1;
        for k in d-i+2..d loop                    -- fill in i-1 values
          ya(k,i) := valpow(ind);
          ind := ind+1;
        end loop;
      end;
    end if;
   -- put_line(file,"The grid of projected samples with the added ones :");
   -- put(file,ya,3);
  end Newton_Samples;

  procedure Newton_Samples 
                 ( ya : in out Matrix;
                   x : in Vector; i : in integer32; trc : in VecVec ) is

  -- DESCRIPTION :
  --   Applies the Newton identities to add extra samples to the
  --   triangular grid of rotated and projected samples.

  -- ON ENTRY :
  --   ya          matrix of projected sample points from triangular grid:
  --               in the i-th column we have values in rows 1..d-i+1;
  --   x           vector of abscisses in the grid;
  --   i           current column in ya to add new samples to;
  --   trc         coefficients of traces trc(j), for j = 1..i-1.

  -- ON RETURN :
  --   ya          updated matrix with projected sample points.

    d : constant integer32 := ya'last(1);

  begin
    if i = 2 then
      ya(d,2) := -Evalc(trc(1).all,x(2));  -- traces come with sign
      for k in ya'first(1)..d-1 loop
        ya(d,2) := ya(d,2) - ya(k,2);
      end loop;
    else
      declare
        valtrc : Vector(1..i-1);  -- values at j-th trace, j=1,..,i-1
        valpow : Vector(1..i-1);  -- values of power sums
        polcff : Vector(0..i-1);  -- coefficients of polynomial
        ind : integer32;
      begin
        for j in valtrc'range loop
          valtrc(j) := Evalc(trc(j).all,x(i));
        end loop;
        valpow := Traces_to_Power_Sums(valtrc);
        for j in valpow'range loop  -- modify the power sums
          for k in 1..d-i+1 loop    -- use values ya(k,i)
            valpow(j) := valpow(j) - ya(k,i)**integer(j);
          end loop;
        end loop;
        valtrc := Power_Sums_to_Traces(valpow);   -- modified traces
        polcff(i-1) := Create(1.0);               -- monic polynomial
        for j in valtrc'range loop    -- fill in coefficients of poly
          polcff(i-1-j) := valtrc(j);
        end loop;
        valpow := Roots(polcff);             -- missing values
        ind := 1;
        for k in d-i+2..d loop                    -- fill in i-1 values
          ya(k,i) := valpow(ind);
          ind := ind+1;
        end loop;
      end;
    end if;
  end Newton_Samples;

-- TEST PROCEDURE :

  procedure Test_Residuals ( file : in file_type; x : in Vector;
                             ya,ps : in Matrix; traces : in VecVec ) is

  -- DESCRIPTION :
  --   While only the first two slices are used in constructing the
  --   first trace, all slices should vanish on the trace.

  -- ON ENTRY :
  --   file       for output of the computations;
  --   x          vector with the abscisses, x'range = 0..degree;
  --   ya         rotated and projected 2nd coordinate of samples;
  --   ps         power sums;
  --   traces     coefficients of the traces.

    sum : Complex_Number;

  begin
    put_line(file,"Evaluating trace 1 in the slices : ");
    for i in x'range loop
      put(file,"At slice "); put(file,i,1); put_line(file," :");
      sum := Evalc(traces(1).all,x(i));
      put(file,"Calculated : "); put(file,sum); new_line(file);
      put(file,"Sampled    : "); put(file,ps(1,i)); new_line(file);
    end loop;
  end Test_Residuals;

-- CREATORS :

  function Create ( grid : Array_of_Standard_Sample_Lists )
                  return Trace_Interpolator1_Rep is

    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := grid'last;
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    res : Trace_Interpolator1_Rep(n,d);
    ya,ps : Matrix(1..d,0..d);
    x : constant Vector(0..d) := Abscisses(grid,natural32(d));

  begin
    res.rot := hyp(hyp'first)(1..n);
    ya := Rotate_Samples(natural32(d),natural32(d),res.rot,grid);
    ps := Power_Sums(ya);
    res.trc := Interpolate(x,ps);
    return res;
  end Create;

  function Create ( grid : Array_of_Standard_Sample_Lists )
                  return Trace_Interpolator1 is

    res : Trace_Interpolator1;

  begin
    res := new Trace_Interpolator1_Rep'(Create(grid));
    return res;
  end Create;

  function Create ( file : file_type; grid : Array_of_Standard_Sample_Lists )
                  return Trace_Interpolator1 is

    res : Trace_Interpolator1;
    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := grid'last;
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    res_rep : Trace_Interpolator1_Rep(n,d);
    ya,ps : Matrix(1..d,0..d);
    x : constant Vector(0..d) := Abscisses(grid,natural32(d));

  begin
    res_rep.rot := hyp(hyp'first)(1..n);
    ya := Rotate_Samples(natural32(d),natural32(d),res_rep.rot,grid);
    ps := Power_Sums(file,ya);
    res_rep.trc := Interpolate(file,x,ps);
    Test_Residuals(file,x,ya,ps,res_rep.trc);
    res := new Trace_Interpolator1_Rep'(res_rep);
    return res;
  end Create;

  function Create ( grid : Array_of_Standard_Sample_Lists; i : integer32 )
                  return Vector is

    n : constant integer32 := Number_of_Variables(Head_Of(grid(0)));
    d : constant integer32 := integer32(Length_Of(grid(0)));
    hyp : constant VecVec := Hyperplane_Sections(Head_Of(grid(0)));
    rot : constant Vector := hyp(hyp'first)(1..n);
    x : constant Vector(0..i) := Abscisses(grid,natural32(i));
    ya : Matrix(1..d,0..i);
    ps,f : Vector(0..i);

  begin
    ya := Rotate_Samples(natural32(d),natural32(i),rot,grid);
    ps := Power_Sums(i,ya);
    f := Create(x,ps);
    return Expand(f,x);
  end Create;

  function Create_on_Triangle 
                  ( file : file_type; grid : Array_of_Standard_Sample_Lists )
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
    for i in 2..d loop
       Newton_Samples(file,ya,x,i,res_rep.trc);
       ps(0..i) := Power_Sums(i,ya);
       f(0..i) := Create(x(0..i),ps(0..i));
       res_rep.trc(i) := new Vector'(Expand(f(0..i),x(0..i)));
    end loop;
    res := new Trace_Interpolator1_Rep'(res_rep);
    return res;
  end Create_on_Triangle;

  function Create_on_Triangle
                  ( grid : Array_of_Standard_Sample_Lists ) 
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
    for i in 2..d loop
       Newton_Samples(ya,x,i,res_rep.trc);
       ps(0..i) := Power_Sums(i,ya);
       f(0..i) := Create(x(0..i),ps(0..i));
       res_rep.trc(i) := new Vector'(Expand(f(0..i),x(0..i)));
    end loop;
    res := new Trace_Interpolator1_Rep'(res_rep);
    return res;
  end Create_on_Triangle;

  function Create ( grid : Stacked_Sample_Grid; d : integer32 )
                  return Trace_Interpolator_Rep is

  -- ASSUMPTION : grid.k > 1.

    res : Trace_Interpolator_Rep(grid.k,integer32(grid.n),d);
    x : constant VecVec(1..grid.k) := Grid_Points(grid);
    y : NesVec(natural32(grid.k)+1,0,d) := Samples_on_Grid(grid,d);

  begin
    res.rot := grid.hyp;
    for i in 1..d loop
      declare
        ps : NesVec(natural32(grid.k),0,i) := Power_Sums(y,i);
        nf : Newton_Form(natural32(grid.k),grid.k,i);
      begin
        nf.x := x;
        nf.f := Create(natural32(grid.k),natural32(i),x,ps);
        res.trc(i) := new Newton_Form'(nf);
        Clear(ps);
      end;
    end loop;
    Clear(y);
    return res;
  end Create;

  function Create ( file : file_type;
                    grid : Stacked_Sample_Grid; d : integer32 )
                  return Trace_Interpolator_Rep is

  -- ASSUMPTION : grid.k > 1.

    res : Trace_Interpolator_Rep(grid.k,integer32(grid.n),d);
    x : constant VecVec(1..grid.k) := Grid_Points(grid);
    y : NesVec(natural32(grid.k)+1,0,d) := Samples_on_Grid(grid,d);

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
        ps : NesVec(natural32(grid.k),0,i) := Power_Sums(y,i);
        nf : Newton_Form(natural32(grid.k),grid.k,i);
        maxerr : double_float;
      begin
        put_line(file,"The power sums :"); put(file,ps);
        nf.x := x;
        nf.f := Create_on_Square(natural32(grid.k),natural32(i),x,ps);
        put_line(file,"The generalized divided differences : ");
        put(file,nf.f);
        Eval_Grid(file,nf.x,ps,nf.f,maxerr);
        res.trc(i) := new Newton_Form'(nf);
        Clear(ps);
      end;
    end loop;
    Clear(y);
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

  function Create ( grid : Stacked_Sample_Grid; d,i : integer32 )
                  return Newton_Form is

    res : Newton_Form(natural32(grid.k),grid.k,i);
    y : NesVec(natural32(grid.k)+1,0,d); -- := Samples_on_Grid(grid,i);
    ps : NesVec(natural32(grid.k),0,d); -- := Power_Sums(y,i);

  begin
    y := Samples_on_Grid(grid,d);
    ps := Power_Sums(y,i);
    res.x := Grid_Points(grid);
    res.f := Create(natural32(grid.k),natural32(i),res.x,ps);
    Clear(y); Clear(ps);
    return res;
  end Create;

  function Create ( file : file_type;
                    grid : Stacked_Sample_Grid; d,i : integer32 )
                  return Newton_Form is

    res : Newton_Form(natural32(grid.k),grid.k,i);
    y : NesVec(natural32(grid.k)+1,0,d); -- := Samples_on_Grid(grid,i);
    ps : NesVec(natural32(grid.k),0,d); -- := Power_Sums(y,i);

  begin
    y := Samples_on_Grid(grid,d);
    put_line(file,"The sampled points : "); put(file,y);
    ps := Power_Sums(y,i);
    put_line(file,"The power sums : "); put(file,ps);
    res.x := Grid_Points(grid);
    put_line(file,"The grid points : "); put_line(file,res.x);
    res.f := Create_on_Square(natural32(grid.k),natural32(i),res.x,ps);
    put_line(file,"The divided differences : "); put(file,res.f);
    return res;
  end Create;

  function Expand ( t : Trace_Interpolator1_Rep ) return Poly is

  -- IMPORTANT NOTICE :
  --   This expansion routine will not work for polynomials like x*y.

    res : Poly := Null_Poly;
    n : constant integer32 := t.n-1;
    px : Poly := Convert(n,t.rot);
    py : Poly := Convert(n,t.rot(1),t.rot(2));
    trace : Poly;
    tt : Term;

  begin
   -- put_line("px = "); put_line(px);
   -- put_line("py = "); put_line(py);
    tt.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    tt.cf := Create(1.0);
    res := Create(tt);
    for i in 1..t.d loop
      Mul(res,py);
      trace := Substitute(n,t.trc(i).all,px);
      Add(res,trace);
      Clear(trace);
    end loop;
   -- put_line("Before normalization : "); put_line(res);
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
      fac.cf := Create(1.0);
      res := Create(fac);
      for i in 1..t.d loop
        tp := Expand(t.trc(i).f,t.trc(i).x);
        Add_Last_Variable(tp);
        fac.dg(t.k+1) := natural32(t.d-i);
        Mul(tp,fac);
        Add(res,tp);
        Clear(tp);
      end loop;
      Clear(fac.dg);
      nv := t.k+1; --t.n - t.k;
      -- put_line("The polynomial before substitution : ");
      -- put_line(res);
      subres := Substitute(nv,t.rot,res);
      Clear(res); res := subres;
      -- put_line("The polynomial before normalization : ");
      -- put_line(res);
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

  function Trace ( t : Trace_Interpolator; i : integer32 )
                 return Newton_Form is
  begin
    return t.trc(i).all;
  end Trace;

-- EVALUATORS :

  function Eval ( t : Trace_Interpolator1_Rep; x : Vector )
                return Complex_Number is

    res : Complex_Number := Create(1.0);   -- trace polynomial is monic
    pt : constant Vector(1..2) := Rotate_and_Project(t.rot,x);

  begin
    for i in 1..t.d loop
      res := res*pt(2) + Evalc(t.trc(i).all,pt(1));
    end loop;
    return res;
  end Eval;

  function Eval ( t : Trace_Interpolator1; x : Vector )
                return Complex_Number is
  begin
    if t = null
     then return Create(0.0);
     else return Eval(t.all,x);
    end if;
  end Eval;

  procedure Eval_Trace ( t : in Vector; d : in integer32; 
                         sps : in Standard_Sample_List;
                         val,eva : out Complex_Number ) is

    hyp : constant VecVec := Hyperplane_Sections(Head_Of(sps));
    tmp : Standard_Sample_List := sps;
    len : constant integer32 := integer32(Length_Of(sps));
    y : Matrix(1..len,1..1);
    x : constant Complex_Number := -hyp(1)(0);

  begin
    for i in 1..len loop
      y(i,1) := Rotate_and_Project2(hyp(1).all,Sample_Point(Head_Of(tmp)).v);
      tmp := Tail_Of(tmp);
    end loop;
    val := Evalc(t,x);
    eva := Power_Sum(1,d,1,y);
    if d mod 2 = 1
     then eva := -eva;
    end if;
  end Eval_Trace;

  procedure Eval_Trace ( t : in Newton_Form; d : in integer32;
                         sps : in Standard_Sample_List;
                         val,eva : out Complex_Number ) is

    hyp : constant VecVec := Hyperplane_Sections(Head_Of(sps));
    tmp : Standard_Sample_List := sps;
    x : constant Vector(hyp'range) := Rotate(hyp,Sample_Point(Head_Of(sps)).v);
    len : constant integer32 := integer32(Length_Of(sps));
    y : Matrix(1..len,1..1);

  begin
    for i in 1..len loop
      y(i,1) := Sample_Point(Head_Of(tmp)).v(t.n+1);
      tmp := Tail_Of(tmp);
    end loop;
    val := Eval0(t.f,t.x,x);
    eva := Power_Sum(1,d,1,y);
    if d mod 2 = 1
     then eva := -eva;
    end if;
  end Eval_Trace;

  function Eval ( t : Trace_Interpolator_Rep; x : Vector )
                return Complex_Number is

    res : Complex_Number;

  begin
    if t.k = 1 then
      res := Eval(t.t1,x);
    else
      res := Create(1.0);
      declare
        pt : constant Vector := Rotate(t.rot,x);
      begin
        for i in 1..t.d loop
          res := res*x(t.k+1) + Eval0(t.trc(i).f,t.trc(i).x,pt);
        end loop;
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
      res := Create(1.0);
      declare
        pt : constant Vector := Rotate(t.rot,x);
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
          res := res*x(t.k+1) + eva;
        end loop;
      end;
    end if;
    return res;
  end Eval;

  function Eval ( t : Trace_Interpolator; x : Vector ) return Complex_Number is
  begin
    if t = null
     then return Create(0.0);
     else return Eval(t.all,x);
    end if;
  end Eval;

  function Eval ( file : file_type;
                  t : Trace_Interpolator; x : Vector ) return Complex_Number is
  begin
    if t = null
     then return Create(0.0);
     else return Eval(file,t.all,x);
    end if;
  end Eval;

-- DIAGNOSTICS :

  function Errors ( t : Trace_Interpolator1;
                    grid : Array_of_Standard_Sample_Lists )
                  return Standard_Floating_Matrices.Matrix is

    res : Standard_Floating_Matrices.Matrix
            (grid'range,1..integer32(Length_Of(grid(grid'first))));
    tmp : Standard_Sample_List;

  begin
    for i in grid'range loop
      tmp := grid(i);
      for j in res'range(2) loop
        res(i,j) := AbsVal(Eval(t,Sample_Point(Head_Of(tmp)).v));
        tmp := Tail_Of(tmp);
        if Is_Null(tmp) then
          for k in j+1..res'last(2) loop
            res(i,k) := 0.0;
          end loop;
          exit;
        end if;
      end loop;
    end loop;
    return res;
  end Errors;

  procedure Write_Errors ( file : in file_type; t : in Trace_Interpolator1;
                           grid : in Array_of_Standard_Sample_Lists ) is

    tmp : Standard_Sample_List;

  begin
    for i in grid'range loop
      tmp := grid(i);
      for j in 1..Length_Of(grid(i)) loop
        put(file,"("); put(file,i,1); put(file,",");
        put(file,j,1); put(file,") : ");
        put(file,Eval(t,Sample_Point(Head_Of(tmp)).v));
        new_line(file);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Write_Errors;

  function Maximal_Error ( residuals : Standard_Floating_Matrices.Matrix )
                         return double_float is

    max : double_float := residuals(residuals'first(1),residuals'first(2));

  begin
    for i in residuals'range(1) loop
      for j in residuals'range(2) loop
        if residuals(i,j) > max
         then max := residuals(i,j);
        end if;
      end loop;
    end loop;
    return max;
  end Maximal_Error;

  function Maximal_Error
             ( t : Trace_Interpolator1;
               grid : Array_of_Standard_Sample_Lists ) return double_float is
  begin
    return Maximal_Error(Errors(t,grid));
  end Maximal_Error;

  function Maximal_Error ( t : Vector; d : integer32;
                           grid : Array_of_Standard_Sample_Lists ) 
                         return double_float is

    res : double_float := 0.0;
    val,eva : Complex_Number;
    abseva : double_float;

  begin
    for i in grid'range loop
      Eval_Trace(t,d,grid(i),val,eva);
      abseva := AbsVal(val-eva);
      if abseva > res
       then res := abseva;
      end if;
    end loop;
    return res;
  end Maximal_Error;

  procedure Write_Errors
              ( file : in file_type; t : in Trace_Interpolator;
                grid : in Stacked_Sample_Grid; maxerr : out double_float ) is

    res : double_float := -1.0;
    tmp : Standard_Sample_List;
    eva : Complex_Number;
    abseva : double_float;

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
           then res := abseva;
          end if;
          tmp := Tail_Of(tmp);
        end loop;
      end loop;
    else
      for i in grid.a'range loop
        Write_Errors(file,t,grid.a(i).all,abseva);
        if abseva > res
         then res := abseva;
        end if;
      end loop;
    end if;
    maxerr := res;
  end Write_Errors;

  function Maximal_Error
             ( t : Trace_Interpolator; grid : Stacked_Sample_Grid )
             return double_float is

    res : double_float := -1.0;
    tmp : Standard_Sample_List;
    eva : Complex_Number;
    abseva : double_float;

  begin
    if grid.k = 1 then
      for i in grid.g'range loop
        tmp := grid.g(i);
        while not Is_Null(tmp) loop
          eva := Eval(t,Sample_Point(Head_Of(tmp)).v);
          abseva := AbsVal(eva);
          if abseva > res
           then res := abseva;
          end if;
          tmp := Tail_Of(tmp);
        end loop;
      end loop;
    else 
      for i in grid.a'range loop
        abseva := Maximal_Error(t,grid.a(i).all);
        if abseva > res
         then res := abseva;
        end if;
      end loop;
    end if;
    return res;
  end Maximal_Error;

  function Maximal_Error
             ( file : file_type;
               t : Newton_Form; d,i : integer32; grid : Stacked_Sample_Grid )
             return double_float is

  -- ALGORITHM :
  --   Run through the grid of samples and evaluate the trace
  --   at the bottom level when all grid values are specified.
  --   So for every list we only do one evaluation of the trace,
  --   predicting the power sum of degree d, which we then compare
  --   with the power sum computed from the actual samples.

    res : double_float := 0.0;
   -- x : constant VecVec(1..grid.k) := Grid_Points(grid);
    y : constant NesVec(natural32(grid.k)+1,0,d) := Samples_on_Grid(grid,d);
    ps : constant NesVec(natural32(grid.k),0,d) := Power_Sums(y,i);

    procedure Eval_Trace_on_Grid
                ( sg : in Stacked_Sample_Grid; v : in NesVec ) is

    -- DESCRIPTION :
    --   sg is the sample grid and v the structure of power sums.

      spt : Standard_Sample;

    begin
      if sg.k = 1 then                     -- this must imply that v.n = 1
        for j in v.v'range loop
          spt := Head_Of(sg.g(j));
          declare
            hyp : constant VecVec := Hyperplane_Sections(spt);
            pt : constant Vector(hyp'range) := Rotate(hyp,Sample_Point(spt).v);
            val : constant Complex_Number := v.v(j);
            eva : constant Complex_Number := Eval0(t.f,t.x,pt);
            abseva : constant double_float := AbsVal(val-eva);
          begin
            put(file,"Power Sum at Samples : "); put(file,val);
            new_line(file);
            put(file,"The Evaluated Trace  : "); put(file,eva);
            new_line(file);
            put(file,"-> Absolute value of difference : ");
            put(file,abseva); new_line(file);
            if abseva > res
             then res := abseva;
            end if;
          end;
        end loop;
      else -- NOTE : for the i-th trace, only power sums up to slice i 
           --   exist in ps, therefore we stop the testing at j = i.
        for j in sg.a'first..i loop
          Eval_Trace_on_Grid(sg.a(j).all,v.w(j).all);
        end loop;
      end if;
    end Eval_Trace_on_Grid;

  begin
    Eval_Trace_on_Grid(grid,ps);
    return res;
  end Maximal_Error;

-- DESTRUCTORS :

  procedure Clear ( t : in out Trace_Interpolator1 ) is

    procedure free is new unchecked_deallocation
      (Trace_Interpolator1_Rep,Trace_Interpolator1);

  begin
    if t /= null
     then Standard_Complex_VecVecs.Clear(t.trc); free(t);
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
     then Clear(nf.all); free(nf);
    end if;
  end Clear;

  procedure Clear ( t : in out Trace_Interpolator ) is

    procedure free is new unchecked_deallocation
      (Trace_Interpolator_Rep,Trace_Interpolator);

  begin
    if t /= null then
      if t.k = 1 then
        Standard_Complex_VecVecs.Clear(t.t1.trc);
      else
        Clear(t.rot);
        for i in t.trc'range loop
          Clear(t.trc(i));
        end loop;
      end if;
      free(t);
    end if;
  end Clear;

end Standard_Trace_Interpolators;
