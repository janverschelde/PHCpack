with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Sample_Point_Lists_io;              use Sample_Point_Lists_io;
with Standard_Lined_Hypersurfaces;       use Standard_Lined_Hypersurfaces;
with Hypersurface_Sample_Grids;          use Hypersurface_Sample_Grids;
with Standard_Stacked_Sample_Grids;      use Standard_Stacked_Sample_Grids;
with Standard_Trace_Interpolators;       use Standard_Trace_Interpolators;
with Standard_Divided_Differences;       use Standard_Divided_Differences;
with Monodromy_Partitions;               use Monodromy_Partitions;
with Monodromy_Polynomial_Breakup;       use Monodromy_Polynomial_Breakup;
with Combinatorial_Factorization;        use Combinatorial_Factorization;
with Factored_Witness_Vectors;           use Factored_Witness_Vectors;
with Certify_Factor_with_Trace;          use Certify_Factor_with_Trace;

package body Multivariate_Factorization is

-- ORGANIZATION :
--   There are three stages in the factorization :
--     1) calculation of the generic points;
--     2) monodromy breakup;
--     3) call for the interpolation routines.
--   In addition, there is the complication of multiple factors.

  procedure Subfactor_with_Multiplicities
               ( n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector;
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in Standard_Complex_Poly_Systems.Poly_Sys;
                 eva_rdp : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicy mu, without output.

    tol : constant double_float := 1.0E-8;
    w0 : constant Standard_Complex_Vectors.Vector
       := Select_Multiple_Factors(m,t,mu);
    wp : constant Standard_Natural_Vectors.Vector := Positions(t,w0,tol);
    subfactors : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    if k > mu then
      subfactors := Init_Factors(natural32(wp'last));
      Monodromy_Breakup(rdp(integer32(mu)),eva_rdp(integer32(mu)),
                        b,v,w0,subfactors);
      if Number_of_Factors(subfactors.all) /= wp'length then
        Assign_Legend(subfactors,wp);
        Merge(factors,subfactors);
      end if;
    end if;
  end Subfactor_with_Multiplicities;

  procedure Subfactor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in Standard_Complex_Poly_Systems.Poly_Sys;
                 eva_rdp : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu,
  --   with intermediate output.

    tol : constant double_float := 1.0E-8;
    w0 : constant Standard_Complex_Vectors.Vector
       := Select_Multiple_Factors(m,t,mu);
    wp : constant Standard_Natural_Vectors.Vector := Positions(t,w0,tol);
    subfactors : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    put(file,"Finding factors of multiplicity ");
    put(file,mu,1); put(file," with sum of degrees = ");
    put(file,k,1); put_line(file," ...");
    put(file,"The positions of those points : ");
    put(file,wp); new_line(file);
    if mu = k then
      put_line(file," -> only one factor.");
    else
      subfactors := Init_Factors(natural32(wp'last));
      Monodromy_Breakup(file,rdp(integer32(mu)),eva_rdp(integer32(mu)),
                        b,v,w0,output,subfactors);
      if Number_of_Factors(subfactors.all) = natural32(wp'length) then
        put_line(file,"No monodromy groupings found.");
      else
        Assign_Legend(subfactors,wp);
        put_line(file,"The factors after legend assignment :");
        Write_Factors(file,subfactors.all);
        Merge(factors,subfactors);
        put_line(file,"The factors after merging : ");
        Write_Factors(file,factors.all);
      end if;
    end if;
  end Subfactor_with_Multiplicities;

  procedure Sub_Trace_Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 mu,k : in natural32;
                 rdp : in Standard_Complex_Poly_Systems.Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu, without output.

    tol : constant double_float := 1.0E-8;
    w0 : constant Standard_Complex_Vectors.Vector
       := Select_Multiple_Factors(m,t,mu);
    wp : constant Standard_Natural_Vectors.Vector := Positions(t,w0,tol);
    subfactors : Standard_Natural_VecVecs.Link_to_VecVec;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    if k > mu then
      Hypersurface_Sample_Grids.Initialize(rdp(integer32(mu)));
      grid := Parallel_Sample1(b,v,w0,2);
      subfactors
        := new Standard_Natural_VecVecs.VecVec'(Factor(w0'length,grid));
      if Number_of_Factors(subfactors.all) /= natural32(wp'length) then
        Assign_Legend(subfactors,wp);
        Merge(factors,subfactors);
      end if;
    end if;
  end Sub_Trace_Factor_with_Multiplicities;

  procedure Sub_Trace_Factor_with_Multiplicities
               ( file : in file_type;
                 n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in Standard_Complex_Poly_Systems.Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu,
  --   with intermediate output.

    tol : constant double_float := 1.0E-8;
    w0 : constant Standard_Complex_Vectors.Vector
       := Select_Multiple_Factors(m,t,mu);
    wp : constant Standard_Natural_Vectors.Vector := Positions(t,w0,tol);
    subfactors : Standard_Natural_VecVecs.Link_to_VecVec;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    put(file,"Finding factors of multiplicity ");
    put(file,mu,1); put(file," with sum of degrees = ");
    put(file,k,1); put_line(file," ...");
    put(file,"The positions of those points : ");
    put(file,wp); new_line(file);
    if mu = k then
      put_line(file," -> only one factor.");
    else
      Hypersurface_Sample_Grids.Initialize(rdp(integer32(mu)));
      grid := Parallel_Sample1(file,false,b,v,w0,2);
      subfactors
        := new Standard_Natural_VecVecs.VecVec'(Factor(file,w0'length,grid));
      if Number_of_Factors(subfactors.all) /= natural32(wp'length) then
        Assign_Legend(subfactors,wp);
        put_line(file,"The factors after legend assignment :");
        Write_Factors(file,subfactors.all);
        Merge(factors,subfactors);
        put_line(file,"The factors after merging : ");
        Write_Factors(file,factors.all);
      end if;
    end if;
  end Sub_Trace_Factor_with_Multiplicities;

  procedure Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out Standard_Complex_Vectors.Link_to_Vector;     
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Performs the factorization when multiplicities > 1, 
  --   without any intermediate output.

  -- ON ENTRY :
  --   n         number of variables of the polynomial;
  --   p         a polynomial in n variables;
  --   ep        nested Horner form of the polynomial;
  --   b         offset vector of a random affine line b + t*v;
  --   v         direction of a random affine line b + t*v;
  --   t         values for t on the random affine line b + t*v;
  --   m         multiplicities sorted in ascending order of t-values;
  --   rdp       contains random derivatives of p,
  --             rdp'last = highest multiplicity in m.

  -- ON RETURN :
  --   factors   labels for generic points grouped together on same factor;
  --   wp        witness points with duplicates removed;
  --   mw        mw(i) is the multiplicity of wp(i).

    tol : constant double_float := 1.0E-8;
    k : natural32;
    eva_rdp : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(rdp'range);
    t1 : constant Standard_Complex_Vectors.Vector
       := Remove_Duplicates(t,tol);
    m1 : constant Standard_Natural_Vectors.Vector
       := Remove_Duplicates(t,tol,m);

  begin
    eva_rdp(1) := ep;
    for i in 2..rdp'last loop
      eva_rdp(i) := Standard_Complex_Poly_Functions.Create(rdp(i));
    end loop;
    factors := Init_Factors(natural32(m1'last));
    for i in 1..rdp'last loop
      k := Countmu(m,natural32(i));
      if k > 0 then
        Subfactor_with_Multiplicities
          (n,d,p,ep,b,v,t1,m1,natural32(i),k,rdp.all,eva_rdp,factors); 
      end if;
    end loop;
    for i in 2..rdp'last loop
      Standard_Complex_Poly_Functions.Clear(eva_rdp(i));
    end loop;
    wp := new Standard_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Factor_with_Multiplicities;

  procedure Factor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out Standard_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Performs the factorization when multiplicities > 1,
  --   with intermediate output to file.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   output    if output needed during continuation;
  --   n         number of variables of the polynomial;
  --   p         a polynomial in n variables;
  --   ep        nested Horner form of the polynomial;
  --   b         offset vector of a random affine line b + t*v;
  --   v         direction of a random affine line b + t*v;
  --   t         values for t on the random affine line b + t*v;
  --   m         multiplicities sorted in ascending order of t-values;
  --   rdp       contains random derivatives of p,
  --             rdp'last = highest multiplicity in m.

  -- ON RETURN :
  --   factors   labels for generic points grouped together on same factor;
  --   wp        witness points with duplicates removed;
  --   mw        mw(i) is the multiplicity of wp(i).

    tol : constant double_float := 1.0E-8;
    k : natural32;
    eva_rdp : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(rdp'range);
    t1 : constant Standard_Complex_Vectors.Vector
       := Remove_Duplicates(t,tol);
    m1 : constant Standard_Natural_Vectors.Vector
       := Remove_Duplicates(t,tol,m);

  begin
    eva_rdp(1) := ep;
    for i in 2..rdp'last loop
      eva_rdp(i) := Standard_Complex_Poly_Functions.Create(rdp(i));
    end loop;
    put(file,"Sorted multiplicities : "); put(file,m);
    put(file," with max = "); put(file,rdp'last,1); new_line(file);
    put_line(file,"The original list of witness points:");
    put_line(file,t);
    put_line(file,"witness point with duplicates removed : ");
    put_line(file,t1);
    put(file,"with corresponding multiplicities : ");
    put(file,m1); new_line(file);
    factors := Init_Factors(natural32(m1'last));
    for i in 1..natural32(rdp'last) loop
      k := Countmu(m,i);
      if k = 0 then
        put(file,"There is no factor with multiplicity ");
        put(file,i,1); put_line(file,".");
      else
        Subfactor_with_Multiplicities
          (file,output,n,d,p,ep,b,v,t1,m1,i,k,rdp.all,eva_rdp,factors); 
      end if;
    end loop;
    for i in 2..rdp'last loop
      Standard_Complex_Poly_Functions.Clear(eva_rdp(i));
    end loop;
    wp := new Standard_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Factor_with_Multiplicities;

  procedure Trace_Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out Standard_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Applies the combinatorial factorization with linear traces
  --   to the case with witness points of higher multiplicity,
  --   without any intermediate output.

  -- ON ENTRY :
  --   n         number of variables of the polynomial;
  --   p         a polynomial in n variables;
  --   b         offset vector of a random affine line b + t*v;
  --   v         direction of a random affine line b + t*v;
  --   t         values for t on the random affine line b + t*v;
  --   m         multiplicities sorted in ascending order of t-values;
  --   rdp       contains random derivatives of p,
  --             rdp'last = highest multiplicity in m.

  -- ON RETURN :
  --   factors   labels for generic points grouped together on same factor;
  --   wp        witness points with duplicates removed;
  --   mw        mw(i) is the multiplicity of wp(i).

    tol : constant double_float := 1.0E-8;
    k : natural32;
    t1 : constant Standard_Complex_Vectors.Vector
       := Remove_Duplicates(t,tol);
    m1 : constant Standard_Natural_Vectors.Vector
       := Remove_Duplicates(t,tol,m);

  begin
    factors := Init_Factors(natural32(m1'last));
    for i in 1..natural32(rdp'last) loop
      k := Countmu(m,i);
      if k > 0 then
        Sub_Trace_Factor_with_Multiplicities
          (n,d,p,b,v,t1,m1,i,k,rdp.all,factors); 
      end if;
    end loop;
    wp := new Standard_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Trace_Factor_with_Multiplicities;

  procedure Trace_Factor_with_Multiplicities
               ( file : in file_type; n,d : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 b,v,t : in Standard_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out Standard_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Applies the combinatorial factorization with linear traces
  --   to the case with witness points of higher multiplicity,
  --   with intermediate output to file.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   n         number of variables of the polynomial;
  --   p         a polynomial in n variables;
  --   ep        nested Horner form of the polynomial;
  --   b         offset vector of a random affine line b + t*v;
  --   v         direction of a random affine line b + t*v;
  --   t         values for t on the random affine line b + t*v;
  --   m         multiplicities sorted in ascending order of t-values;
  --   rdp       contains random derivatives of p,
  --             rdp'last = highest multiplicity in m.

  -- ON RETURN :
  --   factors   labels for generic points grouped together on same factor;
  --   wp        witness points with duplicates removed;
  --   mw        mw(i) is the multiplicity of wp(i).

    tol : constant double_float := 1.0E-8;
    k : natural32;
    t1 : constant Standard_Complex_Vectors.Vector
       := Remove_Duplicates(t,tol);
    m1 : constant Standard_Natural_Vectors.Vector
       := Remove_Duplicates(t,tol,m);

  begin
    put(file,"Sorted multiplicities : "); put(file,m);
    put(file," with max = "); put(file,rdp'last,1); new_line(file);
    put_line(file,"The original list of witness points:");
    put_line(file,t);
    put_line(file,"witness point with duplicates removed : ");
    put_line(file,t1);
    put(file,"with corresponding multiplicities : ");
    put(file,m1); new_line(file);
    factors := Init_Factors(natural32(m1'last));
    for i in 1..natural32(rdp'last) loop
      k := Countmu(m,i);
      if k = 0 then
        put(file,"There is no factor with multiplicity ");
        put(file,i,1); put_line(file,".");
      else
        Sub_Trace_Factor_with_Multiplicities
          (file,n,d,p,b,v,t1,m1,i,k,rdp.all,factors); 
      end if;
    end loop;
    wp := new Standard_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Trace_Factor_with_Multiplicities;

-- AUXILIARY OPERATIONS FOR INTERPOLATION :

  function Interpolate_Factor
              ( p : Standard_Complex_Polynomials.Poly;
                b,v,w : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Finds a polynomial interpolating through a factor of p,
  --   without any intermediate output.

  -- NOTICE : 
  --   For n = 2, the Newton-Taylor form is used, while for n > 2,
  --   we rely on traces to find the expanded form.

    use Standard_Complex_Numbers;

    res : Standard_Complex_Polynomials.Poly;
    n : constant integer32 := b'length;
    d : constant integer32 := w'length;

  begin
    Hypersurface_Sample_Grids.Initialize(p);
    if n = 2 then
      declare
        grid : Array_of_Standard_Sample_Lists(0..d)
             := Parallel_Sample(b,v,w,d);
        q1 : Newton_Interpolator1 := Create(grid,Create(1.0));
        eq : Newton_Form_Evaluator1 := Create(q1);
      begin
        res := Expand(eq);
        Deep_Clear(grid);
        Clear(q1); Clear(eq);
      end;
    else
      declare
        grid : constant Stacked_Sample_Grid(n-1,d) := Full_Sample(b,v,w);
        t : constant Trace_Interpolator := Create(grid,d);
      begin
        res := Expand(t);
      end;
    end if;
    Hypersurface_Sample_Grids.Clear;
    Normalize(res);
    return res;
  end Interpolate_Factor;

  function Interpolate_Factor
              ( file : file_type;
                p : Standard_Complex_Polynomials.Poly;
                b,v,w : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Finds a polynomial interpolating through a factor of p,
  --   with intermediate output.

    use Standard_Complex_Numbers;

    res : Standard_Complex_Polynomials.Poly;
    n : constant integer32 := b'length;
    d : constant integer32 := w'length;
    max_err : double_float;

  begin
    Hypersurface_Sample_Grids.Initialize(p);
    if n = 2 then
      declare
        grid : Array_of_Standard_Sample_Lists(0..d)
             := Parallel_Sample(file,false,b,v,w,d);
        q1 : Newton_Interpolator1 := Create(grid,Create(1.0));
        eq : Newton_Form_Evaluator1 := Create(q1);
      begin
        max_err := Maximal_Error(q1,grid);
        put(file,"The maximal error on the grid : ");
        put(file,max_err,3); new_line(file);
        res := Expand(eq);
        Deep_Clear(grid);
        Clear(q1); Clear(eq);
      end;
    else
      new_line(file);
      put(file,"Computing ");
      put(file,Full_Grid_Size(natural32(n),natural32(d)),1);
      put_line(file," samples...");
      flush(file);
      declare
        grid : constant Stacked_Sample_Grid(n-1,d)
             := Full_Sample(file,b,v,w);
        t : constant Trace_Interpolator := Create(file,grid,d);
      begin
        max_err := Maximal_Error(t,grid);
        put(file,"The maximal error on the grid : ");
        put(file,max_err,3); new_line(file);
        res := Expand(t);
      end;
    end if;
    Hypersurface_Sample_Grids.Clear;
    Normalize(res);
    return res;
  end Interpolate_Factor;

-- TARGET ROUTINES :

  function Leading_Coefficient
              ( p : Standard_Complex_Polynomials.Poly;
                tol : double_float )
              return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the first coefficient in p larger than tol.

    use Standard_Complex_Numbers;
    res : Complex_Number;

    procedure Scan_Term
                ( t : in Standard_Complex_Polynomials.Term;
                  continue : out boolean ) is
    begin
      if AbsVal(t.cf) > tol then
        res := t.cf;
        continue := false;
      else
        continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Leading_Coefficient;

  function Leading_Coefficient
              ( p : DoblDobl_Complex_Polynomials.Poly;
                tol : double_float )
              return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the first coefficient in p larger than tol.

    use DoblDobl_Complex_Numbers;
    res : Complex_Number;

    procedure Scan_Term
                ( t : in DoblDobl_Complex_Polynomials.Term;
                  continue : out boolean ) is

      val : constant double_double := AbsVal(t.cf);
 
    begin
      if val > tol then
        res := t.cf;
        continue := false;
      else
        continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Leading_Coefficient;

  function Leading_Coefficient
              ( p : QuadDobl_Complex_Polynomials.Poly;
                tol : double_float )
              return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the first coefficient in p larger than tol.

    use QuadDobl_Complex_Numbers;
    res : Complex_Number;

    procedure Scan_Term
                ( t : in QuadDobl_Complex_Polynomials.Term;
                  continue : out boolean ) is

      val : constant quad_double := AbsVal(t.cf);
 
    begin
      if val > tol then
        res := t.cf;
        continue := false;
      else
        continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Leading_Coefficient;

  procedure Normalize ( p : in out Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Numbers;

    tol : constant double_float := 1.0E-10;
    leadcff : constant Complex_Number := Leading_Coefficient(p,tol);

    procedure Normalize_Term
                ( t : in out Standard_Complex_Polynomials.Term;
                  continue : out boolean ) is
    begin
      t.cf := t.cf/leadcff;
      continue := true;
    end Normalize_Term;
    procedure Normalize_Terms is
      new Standard_Complex_Polynomials.Changing_Iterator(Normalize_Term);

  begin
    Normalize_Terms(p);
  end Normalize;

  procedure Normalize ( p : in out DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Numbers;

    tol : constant double_float := 1.0E-10;
    leadcff : constant Complex_Number := Leading_Coefficient(p,tol);

    procedure Normalize_Term
                ( t : in out DoblDobl_Complex_Polynomials.Term;
                  continue : out boolean ) is
    begin
      t.cf := t.cf/leadcff;
      continue := true;
    end Normalize_Term;
    procedure Normalize_Terms is
      new DoblDobl_Complex_Polynomials.Changing_Iterator(Normalize_Term);

  begin
    Normalize_Terms(p);
  end Normalize;

  procedure Normalize ( p : in out QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Numbers;

    tol : constant double_float := 1.0E-10;
    leadcff : constant Complex_Number := Leading_Coefficient(p,tol);

    procedure Normalize_Term
                ( t : in out QuadDobl_Complex_Polynomials.Term;
                  continue : out boolean ) is
    begin
      t.cf := t.cf/leadcff;
      continue := true;
    end Normalize_Term;
    procedure Normalize_Terms is
      new QuadDobl_Complex_Polynomials.Changing_Iterator(Normalize_Term);

  begin
    Normalize_Terms(p);
  end Normalize;

  procedure Normalize ( p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Normalize(p(i));
    end loop;
  end Normalize;

  procedure Normalize ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Normalize(p(i));
    end loop;
  end Normalize;

  procedure Normalize ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Normalize(p(i));
    end loop;
  end Normalize;

  procedure Factor ( p : in Standard_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;
                     b,v : out Standard_Complex_Vectors.Vector;
                     wp : out Standard_Complex_Vectors.Link_to_Vector;
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                     rad,dst : out Standard_Floating_Vectors.Vector;
                     fail : out boolean ) is

    genpts : Standard_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    ep : Standard_Complex_Poly_Functions.Eval_Poly
       := Standard_Complex_Poly_Functions.Create(p);
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 20*d;

  begin
    for i in 1..2 loop
      b := Random_Vector(1,integer32(n)); v := Random_Vector(1,integer32(n));
      Generic_Points(p,ep,d,b,v,eps,maxit,genpts,fail,m,rdp,rad,dst);
      exit when not fail;
    end loop;
    if not fail then
      if rdp'last = 1 then
        factors := Init_Factors(d);
        Monodromy_Breakup(p,ep,b,v,genpts,factors);
        wp := new Standard_Complex_Vectors.Vector'(genpts);
        mw := new Standard_Natural_Vectors.Vector'(m);
      else
        Normalize(rdp.all);
        Sort(m,genpts);
        Factor_with_Multiplicities
          (n,d,p,ep,b,v,genpts,m,rdp,factors,wp,mw);
      end if;
    else
      factors := Init_Factors(d);
    end if;
    Remove_Empty_Entries(factors);
    mf := new Standard_Natural_Vectors.Vector'
                (Multiplicity_of_Factors(factors.all,mw.all));
  end Factor;

  procedure Factor ( p : in Standard_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;
                     b,v : out Standard_Complex_Vectors.Vector;
                     wp : out Standard_Complex_Vectors.Link_to_Vector;
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                     fail : out boolean ) is

    rad,dst : Standard_Floating_Vectors.Vector(1..integer32(d));

  begin
    Factor(p,n,d,factors,mf,b,v,wp,mw,rdp,rad,dst,fail);
  end Factor;

  procedure Factor ( file : in file_type; output : in boolean;
                     p : in Standard_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;     
                     b,v : out Standard_Complex_Vectors.Vector;
                     wp : out Standard_Complex_Vectors.Link_to_Vector;     
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                     rad,dst : out Standard_Floating_Vectors.Vector;
                     fail : out boolean ) is

    genpts : Standard_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    ep : Standard_Complex_Poly_Functions.Eval_Poly
       := Standard_Complex_Poly_Functions.Create(p);
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 20*d;

  begin
    for i in 1..2 loop
      b := Random_Vector(1,integer32(n)); v := Random_Vector(1,integer32(n));
      Generic_Points(file,p,ep,d,b,v,eps,maxit,genpts,fail,m,rdp,rad,dst);
      exit when not fail;
    end loop;
    if not fail then
      if rdp'last = 1 then
        factors := Init_Factors(d);
        Monodromy_Breakup(file,p,ep,b,v,genpts,output,factors);
        wp := new Standard_Complex_Vectors.Vector'(genpts);
        mw := new Standard_Natural_Vectors.Vector'(m);
      else
        Normalize(rdp.all);
        Sort(m,genpts);
        Factor_with_Multiplicities
          (file,output,n,d,p,ep,b,v,genpts,m,rdp,factors,wp,mw);
        put(file,"The multiplicity of the factors :");
        put(file,Multiplicity_of_Factors(factors.all,mw.all));
        new_line(file);
      end if;
    else
      factors := Init_Factors(d);
    end if;
    if fail then
      put_line(file,"Failed to converge");
    else
      put_line(file,"The factors : ");
      Write_Factors(file,factors.all,mw.all);
    end if;
    Remove_Empty_Entries(factors);
    mf := new Standard_Natural_Vectors.Vector'
                (Multiplicity_of_Factors(factors.all,mw.all));
    Standard_Complex_Poly_Functions.Clear(ep);
  end Factor;

  procedure Factor ( file : in file_type; output : in boolean;
                     p : in Standard_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;     
                     b,v : out Standard_Complex_Vectors.Vector;
                     wp : out Standard_Complex_Vectors.Link_to_Vector;     
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                     fail : out boolean ) is

    rad,dst : Standard_Floating_Vectors.Vector(1..integer32(d));

  begin
    Factor(file,output,p,n,d,factors,mf,b,v,wp,mw,rdp,rad,dst,fail);
  end Factor;

  procedure Certify ( p : in Standard_Complex_Polynomials.Poly;
                      b,v,wp : in Standard_Complex_Vectors.Vector;
                      mw : in Standard_Natural_Vectors.Vector;
                      f : in Standard_Natural_VecVecs.VecVec;
                      mf : in Standard_Natural_Vectors.Vector;
                      rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                      maxdif : out double_float ) is

    res : double_float;

  begin
    maxdif := 0.0;
    for i in f'range loop
      if mf(i) = 1
       then Certify_Factor(p,b,v,Select_Points(wp,f(i).all),res);
       else Certify_Factor(rdp(integer32(mf(i))),b,v,
                           Select_Points(wp,f(i).all),res);
      end if;
      if res > maxdif
       then maxdif := res;
      end if;
    end loop;
  end Certify;

  procedure Certify ( file : in file_type;
                      p : in Standard_Complex_Polynomials.Poly;
                      b,v,wp : in Standard_Complex_Vectors.Vector;
                      mw : in Standard_Natural_Vectors.Vector;
                      f : in Standard_Natural_VecVecs.VecVec;
                      mf : in Standard_Natural_Vectors.Vector;
                      rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                      maxdif : out double_float ) is

    res : double_float;

  begin
    maxdif := 0.0;
    for i in f'range loop
      put(file,"Certification of factor "); put(file,i,1);
      put_line(file," :");
      if mf(i) = 1
       then Certify_Factor(file,p,b,v,Select_Points(wp,f(i).all),res);
       else Certify_Factor(file,rdp(integer32(mf(i))),b,v,
                           Select_Points(wp,f(i).all),res);
      end if;
      if res > maxdif
       then maxdif := res;
      end if;
    end loop;
  end Certify;

  procedure Interpolate
              ( p : in Standard_Complex_Polynomials.Poly;
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors
                  : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    factors := new Standard_Complex_Poly_Systems.Poly_Sys(f'range);
    for i in f'range loop
      if mf(i) = 1 then
        factors(i) := Interpolate_Factor(p,b,v,Select_Points(wp,f(i).all));
      else
        factors(i)
          := Interpolate_Factor(rdp(integer32(mf(i))),b,v,
                                Select_Points(wp,f(i).all));
      end if;
    end loop;
  end Interpolate;

  procedure Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly;
                b,v,wp : in Standard_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors
                  : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    factors := new Standard_Complex_Poly_Systems.Poly_Sys(f'range);
    for i in f'range loop
      put(file,"Interpolation of factor "); put(file,i,1);
      put_line(file," :");
      if mf(i) = 1 then
        factors(i) 
          := Interpolate_Factor(file,p,b,v,Select_Points(wp,f(i).all));
      else
        factors(i)
          := Interpolate_Factor
               (file,rdp(integer32(mf(i))),b,v,Select_Points(wp,f(i).all));
      end if;
    end loop;
  end Interpolate;

  function Multiply ( factors : Standard_Complex_Poly_Systems.Poly_Sys;
                      mu : Standard_Natural_Vectors.Vector )
                    return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly;

  begin
    Copy(factors(factors'first),res);
    for i in 2..mu(mu'first) loop
      Mul(res,factors(factors'first));
    end loop;
    for i in factors'first+1..factors'last loop
      for j in 1..mu(i) loop
        Mul(res,factors(i));
      end loop;
    end loop;
    return res;
  end Multiply;

  procedure Trace_Factor
              ( p : in Standard_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    ep : Standard_Complex_Poly_Functions.Eval_Poly
       := Standard_Complex_Poly_Functions.Create(p);
    genpts : Standard_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 20*d;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    for i in 1..2 loop
      b := Random_Vector(1,integer32(n)); v := Random_Vector(1,integer32(n));
      Generic_Points(p,ep,d,b,v,eps,maxit,genpts,fail,m,rdp,rad,dst);
      exit when not fail;
    end loop;
    if not fail then
      if rdp'last = 1 then
        Hypersurface_Sample_Grids.Initialize(p);
        grid := Parallel_Sample1(b,v,genpts,2);
        wp := new Standard_Complex_Vectors.Vector'(genpts);
        mw := new Standard_Natural_Vectors.Vector'(m);
        factors := new Standard_Natural_VecVecs.VecVec'(Factor(d,grid));
      else
        Normalize(rdp.all);
        Sort(m,genpts);
        Trace_Factor_with_Multiplicities
          (n,d,p,b,v,genpts,m,rdp,factors,wp,mw);
      end if;
    else 
      factors := Init_Factors(d);
    end if;
    Remove_Empty_Entries(factors);
    mf := new Standard_Natural_Vectors.Vector'
                (Multiplicity_of_Factors(factors.all,mw.all));
  end Trace_Factor;

  procedure Trace_Factor
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out Standard_Complex_Vectors.Vector;
                wp : out Standard_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    ep : Standard_Complex_Poly_Functions.Eval_Poly
       := Standard_Complex_Poly_Functions.Create(p);
    genpts : Standard_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 20*d;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    for i in 1..2 loop
      b := Random_Vector(1,integer32(n)); v := Random_Vector(1,integer32(n));
      Generic_Points(p,ep,d,b,v,eps,maxit,genpts,fail,m,rdp,rad,dst);
      exit when not fail;
    end loop;
    if not fail then
      if rdp'last = 1 then
        Hypersurface_Sample_Grids.Initialize(p);
        grid := Parallel_Sample1(file,false,b,v,genpts,2);
        wp := new Standard_Complex_Vectors.Vector'(genpts);
        mw := new Standard_Natural_Vectors.Vector'(m);
        factors := new Standard_Natural_VecVecs.VecVec'(Factor(file,d,grid));
      else
        Normalize(rdp.all);
        Sort(m,genpts);
        Trace_Factor_with_Multiplicities
          (file,n,d,p,b,v,genpts,m,rdp,factors,wp,mw);
        put(file,"The multiplicity of the factors :");
        put(file,Multiplicity_of_Factors(factors.all,mw.all));
        new_line(file);
      end if;
    else
      factors := Init_Factors(d);
    end if;
    Remove_Empty_Entries(factors);
    mf := new Standard_Natural_Vectors.Vector'
                (Multiplicity_of_Factors(factors.all,mw.all));
  end Trace_Factor;

end Multivariate_Factorization; 
