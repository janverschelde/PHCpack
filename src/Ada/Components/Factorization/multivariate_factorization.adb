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
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with DoblDobl_Random_Vectors;            use DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_SysFun;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Sample_Point_Lists_io;              use Sample_Point_Lists_io;
with DoblDobl_Sample_Lists;              use DoblDobl_Sample_Lists;
with QuadDobl_Sample_Lists;              use QuadDobl_Sample_Lists;
with Standard_Lined_Hypersurfaces;       use Standard_Lined_Hypersurfaces;
with DoblDobl_Lined_Hypersurfaces;       use DoblDobl_Lined_Hypersurfaces;
with QuadDobl_Lined_Hypersurfaces;       use QuadDobl_Lined_Hypersurfaces;
with Hypersurface_Sample_Grids;          use Hypersurface_Sample_Grids;
with DoblDobl_Gridded_Hypersurfaces;     use DoblDobl_Gridded_Hypersurfaces;
with QuadDobl_Gridded_Hypersurfaces;     use QuadDobl_Gridded_Hypersurfaces;
with Standard_Stacked_Sample_Grids;      use Standard_Stacked_Sample_Grids;
with Standard_Trace_Interpolators;       use Standard_Trace_Interpolators;
with Standard_Divided_Differences;       use Standard_Divided_Differences;
with Monodromy_Partitions;               use Monodromy_Partitions;
with Monodromy_Polynomial_Breakup;       use Monodromy_Polynomial_Breakup;
with Combinatorial_Factorization;        use Combinatorial_Factorization;
with Factored_Witness_Vectors;           use Factored_Witness_Vectors;
with Certify_Factor_with_Trace;          use Certify_Factor_with_Trace;
with Interpolate_Multivariate_Factor;    use Interpolate_Multivariate_Factor;

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
               ( n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector;
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 eva_rdp : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicy mu, without output.

    tol : constant double_float := 1.0E-8;
    w0 : constant DoblDobl_Complex_Vectors.Vector
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
               ( n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector;
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 eva_rdp : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicy mu, without output.

    tol : constant double_float := 1.0E-8;
    w0 : constant QuadDobl_Complex_Vectors.Vector
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

  procedure Subfactor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 eva_rdp : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu,
  --   with intermediate output.

    tol : constant double_float := 1.0E-8;
    w0 : constant DoblDobl_Complex_Vectors.Vector
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

  procedure Subfactor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 eva_rdp : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu,
  --   with intermediate output.

    tol : constant double_float := 1.0E-8;
    w0 : constant QuadDobl_Complex_Vectors.Vector
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
               ( n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 mu,k : in natural32;
                 rdp : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu, without output.

    tol : constant double_float := 1.0E-8;
    w0 : constant DoblDobl_Complex_Vectors.Vector
       := Select_Multiple_Factors(m,t,mu);
    wp : constant Standard_Natural_Vectors.Vector := Positions(t,w0,tol);
    subfactors : Standard_Natural_VecVecs.Link_to_VecVec;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    if k > mu then
      DoblDobl_Gridded_Hypersurfaces.Initialize(rdp(integer32(mu)));
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
               ( n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 mu,k : in natural32;
                 rdp : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu, without output.

    tol : constant double_float := 1.0E-8;
    w0 : constant QuadDobl_Complex_Vectors.Vector
       := Select_Multiple_Factors(m,t,mu);
    wp : constant Standard_Natural_Vectors.Vector := Positions(t,w0,tol);
    subfactors : Standard_Natural_VecVecs.Link_to_VecVec;
    grid : Array_of_QuadDobl_Sample_Lists(0..2);

  begin
    if k > mu then
      QuadDobl_Gridded_Hypersurfaces.Initialize(rdp(integer32(mu)));
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

  procedure Sub_Trace_Factor_with_Multiplicities
               ( file : in file_type;
                 n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu,
  --   with intermediate output.

    tol : constant double_float := 1.0E-8;
    w0 : constant DoblDobl_Complex_Vectors.Vector
       := Select_Multiple_Factors(m,t,mu);
    wp : constant Standard_Natural_Vectors.Vector := Positions(t,w0,tol);
    subfactors : Standard_Natural_VecVecs.Link_to_VecVec;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    put(file,"Finding factors of multiplicity ");
    put(file,mu,1); put(file," with sum of degrees = ");
    put(file,k,1); put_line(file," ...");
    put(file,"The positions of those points : ");
    put(file,wp); new_line(file);
    if mu = k then
      put_line(file," -> only one factor.");
    else
      DoblDobl_Gridded_Hypersurfaces.Initialize(rdp(integer32(mu)));
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

  procedure Sub_Trace_Factor_with_Multiplicities
               ( file : in file_type;
                 n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector; mu,k : in natural32;
                 rdp : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 factors : in Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Subroutine to factor one component of multiplicity mu,
  --   with intermediate output.

    tol : constant double_float := 1.0E-8;
    w0 : constant QuadDobl_Complex_Vectors.Vector
       := Select_Multiple_Factors(m,t,mu);
    wp : constant Standard_Natural_Vectors.Vector := Positions(t,w0,tol);
    subfactors : Standard_Natural_VecVecs.Link_to_VecVec;
    grid : Array_of_QuadDobl_Sample_Lists(0..2);

  begin
    put(file,"Finding factors of multiplicity ");
    put(file,mu,1); put(file," with sum of degrees = ");
    put(file,k,1); put_line(file," ...");
    put(file,"The positions of those points : ");
    put(file,wp); new_line(file);
    if mu = k then
      put_line(file," -> only one factor.");
    else
      QuadDobl_Gridded_Hypersurfaces.Initialize(rdp(integer32(mu)));
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

-- FACTORING WITH GIVEN GENERIC POINTS :

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
               ( n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out DoblDobl_Complex_Vectors.Link_to_Vector;     
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

    tol : constant double_float := 1.0E-8;
    k : natural32;
    eva_rdp : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(rdp'range);
    t1 : constant DoblDobl_Complex_Vectors.Vector
       := Remove_Duplicates(t,tol);
    m1 : constant Standard_Natural_Vectors.Vector
       := Remove_Duplicates(t,tol,m);

  begin
    eva_rdp(1) := ep;
    for i in 2..rdp'last loop
      eva_rdp(i) := DoblDobl_Complex_Poly_Functions.Create(rdp(i));
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
      DoblDobl_Complex_Poly_Functions.Clear(eva_rdp(i));
    end loop;
    wp := new DoblDobl_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Factor_with_Multiplicities;

  procedure Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out QuadDobl_Complex_Vectors.Link_to_Vector;     
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

    tol : constant double_float := 1.0E-8;
    k : natural32;
    eva_rdp : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(rdp'range);
    t1 : constant QuadDobl_Complex_Vectors.Vector
       := Remove_Duplicates(t,tol);
    m1 : constant Standard_Natural_Vectors.Vector
       := Remove_Duplicates(t,tol,m);

  begin
    eva_rdp(1) := ep;
    for i in 2..rdp'last loop
      eva_rdp(i) := QuadDobl_Complex_Poly_Functions.Create(rdp(i));
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
      QuadDobl_Complex_Poly_Functions.Clear(eva_rdp(i));
    end loop;
    wp := new QuadDobl_Complex_Vectors.Vector'(t1);
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

  procedure Factor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

    tol : constant double_float := 1.0E-8;
    k : natural32;
    eva_rdp : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(rdp'range);
    t1 : constant DoblDobl_Complex_Vectors.Vector
       := Remove_Duplicates(t,tol);
    m1 : constant Standard_Natural_Vectors.Vector
       := Remove_Duplicates(t,tol,m);

  begin
    eva_rdp(1) := ep;
    for i in 2..rdp'last loop
      eva_rdp(i) := DoblDobl_Complex_Poly_Functions.Create(rdp(i));
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
      DoblDobl_Complex_Poly_Functions.Clear(eva_rdp(i));
    end loop;
    wp := new DoblDobl_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Factor_with_Multiplicities;

  procedure Factor_with_Multiplicities
               ( file : in file_type; output : in boolean;
                 n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

    tol : constant double_float := 1.0E-8;
    k : natural32;
    eva_rdp : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(rdp'range);
    t1 : constant QuadDobl_Complex_Vectors.Vector
       := Remove_Duplicates(t,tol);
    m1 : constant Standard_Natural_Vectors.Vector
       := Remove_Duplicates(t,tol,m);

  begin
    eva_rdp(1) := ep;
    for i in 2..rdp'last loop
      eva_rdp(i) := QuadDobl_Complex_Poly_Functions.Create(rdp(i));
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
      QuadDobl_Complex_Poly_Functions.Clear(eva_rdp(i));
    end loop;
    wp := new QuadDobl_Complex_Vectors.Vector'(t1);
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
               ( n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

    tol : constant double_float := 1.0E-8;
    k : natural32;
    t1 : constant DoblDobl_Complex_Vectors.Vector
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
    wp := new DoblDobl_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Trace_Factor_with_Multiplicities;

  procedure Trace_Factor_with_Multiplicities
               ( n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

    tol : constant double_float := 1.0E-8;
    k : natural32;
    t1 : constant QuadDobl_Complex_Vectors.Vector
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
    wp := new QuadDobl_Complex_Vectors.Vector'(t1);
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

  procedure Trace_Factor_with_Multiplicities
               ( file : in file_type; n,d : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,t : in DoblDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

    tol : constant double_float := 1.0E-8;
    k : natural32;
    t1 : constant DoblDobl_Complex_Vectors.Vector
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
    wp := new DoblDobl_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Trace_Factor_with_Multiplicities;

  procedure Trace_Factor_with_Multiplicities
               ( file : in file_type; n,d : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,t : in QuadDobl_Complex_Vectors.Vector; 
                 m : in Standard_Natural_Vectors.Vector;
                 rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 factors : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                 mw : out Standard_Natural_Vectors.Link_to_Vector ) is

    tol : constant double_float := 1.0E-8;
    k : natural32;
    t1 : constant QuadDobl_Complex_Vectors.Vector
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
    wp := new QuadDobl_Complex_Vectors.Vector'(t1);
    mw := new Standard_Natural_Vectors.Vector'(m1);
  end Trace_Factor_with_Multiplicities;

-- FACTORING WITHOUT GIVEN GENERIC POINTS :

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

  procedure Factor ( p : in DoblDobl_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;
                     b,v : out DoblDobl_Complex_Vectors.Vector;
                     wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                     rad,dst : out Standard_Floating_Vectors.Vector;
                     fail : out boolean ) is

    genpts : DoblDobl_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    ep : DoblDobl_Complex_Poly_Functions.Eval_Poly
       := DoblDobl_Complex_Poly_Functions.Create(p);
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
        wp := new DoblDobl_Complex_Vectors.Vector'(genpts);
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

  procedure Factor ( p : in QuadDobl_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;
                     b,v : out QuadDobl_Complex_Vectors.Vector;
                     wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                     rad,dst : out Standard_Floating_Vectors.Vector;
                     fail : out boolean ) is

    genpts : QuadDobl_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    ep : QuadDobl_Complex_Poly_Functions.Eval_Poly
       := QuadDobl_Complex_Poly_Functions.Create(p);
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
        wp := new QuadDobl_Complex_Vectors.Vector'(genpts);
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

  procedure Factor ( p : in DoblDobl_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;
                     b,v : out DoblDobl_Complex_Vectors.Vector;
                     wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                     fail : out boolean ) is

    rad,dst : Standard_Floating_Vectors.Vector(1..integer32(d));

  begin
    Factor(p,n,d,factors,mf,b,v,wp,mw,rdp,rad,dst,fail);
  end Factor;

  procedure Factor ( p : in QuadDobl_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;
                     b,v : out QuadDobl_Complex_Vectors.Vector;
                     wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
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
                     p : in DoblDobl_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;     
                     b,v : out DoblDobl_Complex_Vectors.Vector;
                     wp : out DoblDobl_Complex_Vectors.Link_to_Vector;     
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                     rad,dst : out Standard_Floating_Vectors.Vector;
                     fail : out boolean ) is

    genpts : DoblDobl_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    ep : DoblDobl_Complex_Poly_Functions.Eval_Poly
       := DoblDobl_Complex_Poly_Functions.Create(p);
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
        wp := new DoblDobl_Complex_Vectors.Vector'(genpts);
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
    DoblDobl_Complex_Poly_Functions.Clear(ep);
  end Factor;

  procedure Factor ( file : in file_type; output : in boolean;
                     p : in QuadDobl_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;     
                     b,v : out QuadDobl_Complex_Vectors.Vector;
                     wp : out QuadDobl_Complex_Vectors.Link_to_Vector;     
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                     rad,dst : out Standard_Floating_Vectors.Vector;
                     fail : out boolean ) is

    genpts : QuadDobl_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    ep : QuadDobl_Complex_Poly_Functions.Eval_Poly
       := QuadDobl_Complex_Poly_Functions.Create(p);
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
        wp := new QuadDobl_Complex_Vectors.Vector'(genpts);
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
    QuadDobl_Complex_Poly_Functions.Clear(ep);
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

  procedure Factor ( file : in file_type; output : in boolean;
                     p : in DoblDobl_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;     
                     b,v : out DoblDobl_Complex_Vectors.Vector;
                     wp : out DoblDobl_Complex_Vectors.Link_to_Vector;     
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                     fail : out boolean ) is

    rad,dst : Standard_Floating_Vectors.Vector(1..integer32(d));

  begin
    Factor(file,output,p,n,d,factors,mf,b,v,wp,mw,rdp,rad,dst,fail);
  end Factor;

  procedure Factor ( file : in file_type; output : in boolean;
                     p : in QuadDobl_Complex_Polynomials.Poly;
                     n,d : in natural32;
                     factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                     mf : out Standard_Natural_Vectors.Link_to_Vector;     
                     b,v : out QuadDobl_Complex_Vectors.Vector;
                     wp : out QuadDobl_Complex_Vectors.Link_to_Vector;     
                     mw : out Standard_Natural_Vectors.Link_to_Vector;
                     rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
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

  procedure Certify ( p : in DoblDobl_Complex_Polynomials.Poly;
                      b,v,wp : in DoblDobl_Complex_Vectors.Vector;
                      mw : in Standard_Natural_Vectors.Vector;
                      f : in Standard_Natural_VecVecs.VecVec;
                      mf : in Standard_Natural_Vectors.Vector;
                      rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
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

  procedure Certify ( p : in QuadDobl_Complex_Polynomials.Poly;
                      b,v,wp : in QuadDobl_Complex_Vectors.Vector;
                      mw : in Standard_Natural_Vectors.Vector;
                      f : in Standard_Natural_VecVecs.VecVec;
                      mf : in Standard_Natural_Vectors.Vector;
                      rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
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

  procedure Certify ( file : in file_type;
                      p : in DoblDobl_Complex_Polynomials.Poly;
                      b,v,wp : in DoblDobl_Complex_Vectors.Vector;
                      mw : in Standard_Natural_Vectors.Vector;
                      f : in Standard_Natural_VecVecs.VecVec;
                      mf : in Standard_Natural_Vectors.Vector;
                      rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
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

  procedure Certify ( file : in file_type;
                      p : in QuadDobl_Complex_Polynomials.Poly;
                      b,v,wp : in QuadDobl_Complex_Vectors.Vector;
                      mw : in Standard_Natural_Vectors.Vector;
                      f : in Standard_Natural_VecVecs.VecVec;
                      mf : in Standard_Natural_Vectors.Vector;
                      rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
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
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                b,v,wp : in DoblDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors
                  : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    factors := new DoblDobl_Complex_Poly_Systems.Poly_Sys(f'range);
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
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                b,v,wp : in QuadDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors
                  : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    factors := new QuadDobl_Complex_Poly_Systems.Poly_Sys(f'range);
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

  procedure Interpolate
              ( file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly;
                b,v,wp : in DoblDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors
                  : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    factors := new DoblDobl_Complex_Poly_Systems.Poly_Sys(f'range);
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

  procedure Interpolate
              ( file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly;
                b,v,wp : in QuadDobl_Complex_Vectors.Vector;
                mw : in Standard_Natural_Vectors.Vector;
                f : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Natural_Vectors.Vector;
                rdp : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                factors
                  : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    factors := new QuadDobl_Complex_Poly_Systems.Poly_Sys(f'range);
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

  function Multiply ( factors : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                      mu : Standard_Natural_Vectors.Vector )
                    return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;

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

  function Multiply ( factors : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                      mu : Standard_Natural_Vectors.Vector )
                    return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

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
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out DoblDobl_Complex_Vectors.Vector;
                wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    ep : DoblDobl_Complex_Poly_Functions.Eval_Poly
       := DoblDobl_Complex_Poly_Functions.Create(p);
    genpts : DoblDobl_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 20*d;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    for i in 1..2 loop
      b := Random_Vector(1,integer32(n)); v := Random_Vector(1,integer32(n));
      Generic_Points(p,ep,d,b,v,eps,maxit,genpts,fail,m,rdp,rad,dst);
      exit when not fail;
    end loop;
    if not fail then
      if rdp'last = 1 then
        DoblDobl_Gridded_Hypersurfaces.Initialize(p);
        grid := Parallel_Sample1(b,v,genpts,2);
        wp := new DoblDobl_Complex_Vectors.Vector'(genpts);
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
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out QuadDobl_Complex_Vectors.Vector;
                wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    ep : QuadDobl_Complex_Poly_Functions.Eval_Poly
       := QuadDobl_Complex_Poly_Functions.Create(p);
    genpts : QuadDobl_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 20*d;
    grid : Array_of_QuadDobl_Sample_Lists(0..2);

  begin
    for i in 1..2 loop
      b := Random_Vector(1,integer32(n)); v := Random_Vector(1,integer32(n));
      Generic_Points(p,ep,d,b,v,eps,maxit,genpts,fail,m,rdp,rad,dst);
      exit when not fail;
    end loop;
    if not fail then
      if rdp'last = 1 then
        QuadDobl_Gridded_Hypersurfaces.Initialize(p);
        grid := Parallel_Sample1(b,v,genpts,2);
        wp := new QuadDobl_Complex_Vectors.Vector'(genpts);
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

  procedure Trace_Factor
              ( file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out DoblDobl_Complex_Vectors.Vector;
                wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    ep : DoblDobl_Complex_Poly_Functions.Eval_Poly
       := DoblDobl_Complex_Poly_Functions.Create(p);
    genpts : DoblDobl_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 20*d;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    for i in 1..2 loop
      b := Random_Vector(1,integer32(n)); v := Random_Vector(1,integer32(n));
      Generic_Points(p,ep,d,b,v,eps,maxit,genpts,fail,m,rdp,rad,dst);
      exit when not fail;
    end loop;
    if not fail then
      if rdp'last = 1 then
        DoblDobl_Gridded_Hypersurfaces.Initialize(p);
        grid := Parallel_Sample1(file,false,b,v,genpts,2);
        wp := new DoblDobl_Complex_Vectors.Vector'(genpts);
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

  procedure Trace_Factor
              ( file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly;
                n,d : in natural32;
                factors : out Standard_Natural_VecVecs.Link_to_VecVec;
                mf : out Standard_Natural_Vectors.Link_to_Vector;
                b,v : out QuadDobl_Complex_Vectors.Vector;
                wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                mw : out Standard_Natural_Vectors.Link_to_Vector;
                rdp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                rad,dst : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    ep : QuadDobl_Complex_Poly_Functions.Eval_Poly
       := QuadDobl_Complex_Poly_Functions.Create(p);
    genpts : QuadDobl_Complex_Vectors.Vector(1..integer32(d));
    m : Standard_Natural_Vectors.Vector(1..integer32(d));
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 20*d;
    grid : Array_of_QuadDobl_Sample_Lists(0..2);

  begin
    for i in 1..2 loop
      b := Random_Vector(1,integer32(n)); v := Random_Vector(1,integer32(n));
      Generic_Points(p,ep,d,b,v,eps,maxit,genpts,fail,m,rdp,rad,dst);
      exit when not fail;
    end loop;
    if not fail then
      if rdp'last = 1 then
        QuadDobl_Gridded_Hypersurfaces.Initialize(p);
        grid := Parallel_Sample1(file,false,b,v,genpts,2);
        wp := new QuadDobl_Complex_Vectors.Vector'(genpts);
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
