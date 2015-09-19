with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Natural_Vectors;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Random_Polynomials;       use Standard_Random_Polynomials;
with Standard_Stacked_Sample_Grids;     use Standard_Stacked_Sample_Grids;
with Hypersurface_Sample_Grids;         use Hypersurface_Sample_Grids;
with Standard_Lined_Hypersurfaces;      use Standard_Lined_Hypersurfaces;
with Standard_Trace_Interpolators;      use Standard_Trace_Interpolators;

procedure ts_trapol is

-- DESCRIPTION :
--   Test on trace form of the interpolators for one polynomial.

  procedure Read_Polynomial ( n : out natural32; p : out Poly ) is

  -- DESCRIPTION :
  --   Allows the user to enter a polynomial or reads from file.

    ans : character;
    file : file_type;

  begin
    n := 0;
    new_line;
    put("Is the polynomial on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the file for the polynomial...");
      Read_Name_and_Open_File(file);
      get(file,n);
      Symbol_Table.Init(n);
      get(file,p);
    else
      new_line;
      put("Give the number of variables : ");
      get(n);
      Symbol_Table.Init(n);
      put("Give your polynomial : ");
      get(p);
      skip_line;
    end if;
    new_line;
    put("Your polynomial is "); put(p); new_line;
  end Read_Polynomial;

  procedure Generate_Polynomial ( n : out natural32; p : out Poly ) is

  -- DESCRIPTION :
  --   Generates a polynomial in n variables with complex coefficients.

    d,m : natural32 := 0;

  begin
    n := 0;
    new_line;
    put("Give the number of variables : "); get(n);
    put("Give the degree bound : "); get(d);
    put("Give the number of terms : "); get(m);
    p := Random_Sparse_Poly(n,d,m,0);
    put_line("A random polynomial : ");
    put_line(p);
    skip_line;
  end Generate_Polynomial;

  procedure Witness_Points
              ( file : in file_type;
                n,d : in integer32; p : in Poly; ep : in Eval_Poly;
                b,v,w : out Standard_Complex_Vectors.Vector;
                mu : out Standard_Natural_Vectors.Vector;
                rdp : out Link_to_Poly_Sys; fail : out boolean ) is

  -- DESCRIPTION :
  --   Computes witness points for the given polynomial,
  --   writing intermediate results to a file.

    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 10*natural32(d);

  begin
    b := Random_Vector(1,n);
    v := Random_Vector(1,n);
    Generic_Points(file,p,ep,natural32(d),b,v,eps,maxit,w,fail,mu,rdp);
  end Witness_Points;

  procedure Interpolate ( file : in file_type;
                          n : in integer32; p : in Poly ) is

    d : constant integer32 := Degree(p);
    ep : constant Eval_Poly := Create(p);
    grid : Stacked_Sample_Grid(n-1,d);
    fail : boolean;
    b,v : Standard_Complex_Vectors.Vector(1..n);
    w : Standard_Complex_Vectors.Vector(1..d);
    m : Standard_Natural_Vectors.Vector(1..d);
    t : Trace_Interpolator;
    rdp : Link_to_Poly_Sys;
    max_err : double_float;
    ip : Poly;

  begin
    put(file,"p = "); put(file,p); new_line(file);
    Witness_Points(file,n,d,p,ep,b,v,w,m,rdp,fail);
    Hypersurface_Sample_Grids.Initialize(p);
    grid := Full_Sample(file,b,v,w);
    Hypersurface_Sample_Grids.Clear;
    max_err := Maximal_Error(grid);
    put(file,"Maximal error on grid : ");
    put(file,max_err,3); new_line(file);
    Write_Grid_Values(file,grid);
    t := Create(grid,d);
    Write_Errors(file,t,grid,max_err);
    put(file,"Maximal residual of interpolator : ");
    put(file,max_err,3); new_line(file);
    ip := Expand(t);
    put_line(file,"The expanded interpolating polynomial : ");
    put_line(file,ip);
    new_line;
    put_line("The expanded interpolating polynomial : ");
    put_line(ip);
  end Interpolate;

  procedure Main is

    n : natural32;
    p : Poly;
    ans : character;
    file : file_type;
    first_time : boolean := true;

  begin
    new_line;
    put_line("Interpolating Multivariate Polynomials with Traces");
    loop
      new_line;
      put_line("Choose one of the following options : ");
      put_line("  0. Exit this program.");
      put_line("  1. Interpolate a randomly generated complex polynomial.");
      put_line("  2. Give your own polynomial to interpolate.");
      put("Type 0, 1, or 2 to choose : ");
      Ask_Alternative(ans,"012");
      case ans is
        when '0' => return;
        when '1' => Generate_Polynomial(n,p);
        when '2' => Read_Polynomial(n,p);
        when others => null;
      end case;
      new_line;
      if first_time then
        put_line("Reading the name of the output file...");
        Read_Name_and_Create_File(file);
        first_time := false;
      end if;
      Interpolate(file,integer32(n),p);
    end loop;
  end Main;

begin
  Main;
end ts_trapol;
