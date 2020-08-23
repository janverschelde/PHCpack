with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Standard_Complex_Matrices;
with Multprec_Complex_Matrices;
with Standard_Random_Vectors;
with Standard_Random_Matrices;
with Multprec_Random_Vectors;
with Multprec_Random_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_to_Multprec_Convertors;   use Standard_to_Multprec_Convertors;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions;
with Standard_Deflate_Singularities;
with Multprec_Deflate_Singularities;
with Standard_Deflation_Trees_io;       use Standard_Deflation_Trees_io;
with Multprec_Deflation_Trees_io;       use Multprec_Deflation_Trees_io;
with Standard_Deflation_Methods;
with DoblDobl_Deflation_Methods;
with QuadDobl_Deflation_Methods;
with Multprec_Deflation_Methods;

package body Drivers_to_Deflate_Singularities is

  procedure Read_Tolerance ( tol : in out double_float ) is

    ans : character;

  begin
    new_line;
    loop
      put("Current tolerance for numerical rank is ");
      put(tol,3); new_line;
      put("Do you wish to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      put("Give new tolerance : "); get(tol);
    end loop;
  end Read_Tolerance;

  function Test_Standard_Deflate
              ( p : Standard_Complex_Poly_Systems.Poly_Sys; k,m : natural32 )
              return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system p, at deflation number k.

    use Standard_Complex_Vectors,Standard_Random_Vectors;
    use Standard_Complex_Matrices,Standard_Random_Matrices;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Polynomials_io;
    use Standard_Complex_Poly_Systems;
    use Standard_Deflate_Singularities;

    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    ne : constant integer32 := p'last-p'first+1;
    cf : constant Vector(1..integer32(m)) := Random_Vector(1,integer32(m));
    a : constant Matrix(1..nv,1..integer32(m))
      := Random_Matrix(natural32(nv),m);
    dp : constant Poly_Sys(p'first..p'last+ne+1) := Deflate(p,m,a,cf);

  begin
    Add_Multiplier_Symbols(k,m);
    put(natural32(dp'last),natural32(nv)+m,dp(p'range));
    for i in p'last+1..dp'last loop
      put_line(dp(i));
    end loop;
    return dp;
  end Test_Standard_Deflate;

  function Test_Multprec_Deflate
              ( p : Multprec_Complex_Poly_Systems.Poly_Sys;
                k,m,size : natural32 )
              return Multprec_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system p, at deflation number k.

    use Multprec_Complex_Vectors,Multprec_Random_Vectors;
    use Multprec_Complex_Matrices,Multprec_Random_Matrices;
    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Polynomials_io;
    use Multprec_Complex_Poly_Systems;
    use Multprec_Deflate_Singularities;

    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    ne : constant integer32 := p'last-p'first+1;
    cf : constant Vector(1..integer32(m))
       := Random_Vector(1,integer32(m),size);
    a : constant Matrix(1..nv,1..integer32(m))
      := Random_Matrix(natural32(nv),m,size);
    dp : constant Poly_Sys(p'first..p'last+ne+1) := Deflate(p,m,a,cf);

  begin
    Add_Multiplier_Symbols(k,m);
    put(dp(p'range));
    for i in p'last+1..dp'last loop
      put_line(dp(i));
    end loop;
    return dp;
  end Test_Multprec_Deflate;

  procedure Multiple_Standard_Deflations
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    ans : character;
    wp : Link_to_Poly_Sys := new Poly_Sys'(p);
    k : natural32 := 1;
    m : natural32 := 0;

  begin
    loop
      put("Give number of multipliers : "); get(m);
      declare
        dp : Link_to_Poly_Sys
           := new Poly_Sys'(Test_Standard_Deflate(wp.all,k,m));
      begin
        put("Do you want to deflate once more ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        Clear(wp);
        wp := new Poly_Sys(dp'range);
        for i in dp'range loop
          Copy(dp(i),wp(i));
        end loop;
        put(file,"DEFLATED SYSTEM #"); put(file,k,1); put_line(file," :");
        Write_Deflated_System(file,p,dp.all);
        Clear(dp);
        k := k + 1;
      end;
    end loop;
  end Multiple_Standard_Deflations;

  procedure Multiple_Multprec_Deflations
               ( file : in file_type;
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 size : in natural32 ) is

    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;

    ans : character;
    wp : Link_to_Poly_Sys := new Poly_Sys'(p);
    k : natural32 := 1;
    m : natural32 := 0;

  begin
    loop
      put("Give number of multipliers : "); get(m);
      declare
        dp : Link_to_Poly_Sys
           := new Poly_Sys'(Test_Multprec_Deflate(wp.all,k,m,size));
      begin
        put("Do you want to deflate once more ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        Clear(wp);
        wp := new Poly_Sys(dp'range);
        for i in dp'range loop
          Copy(dp(i),wp(i));
        end loop;
        put(file,"DEFLATED SYSTEM #"); put(file,k,1); put_line(file," :");
        Write_Deflated_System(file,p,dp.all);
        Clear(dp);
        k := k + 1;
      end;
    end loop;
  end Multiple_Multprec_Deflations;

  procedure Set_Default_Parameters
              ( symbolic,output : out boolean;
                maxitr,maxdef,nbdgts : out natural32;
	        tolerr,tolres,tolrnk : out double_float ) is
  begin
    symbolic := false;
    output := false;
    maxitr := 3;
    maxdef := 3;
    nbdgts := 15;
    tolerr := 1.0E-12;
    tolres := 1.0E-12;
    tolrnk := 1.0E-6;
  end Set_Default_Parameters;

  procedure Display_Parameters
              ( file : in file_type;
                symbolic,output : in boolean;
                maxitr,maxdef,nbdgts : in natural32;
	        tolerr,tolres,tolrnk : in double_float ) is
  begin
    put_line(file,"MENU with current settings for the deflation method : ");
    if symbolic then
      put(file,"  1. Symbolic evaluation of Jacobian matrix : ");
    else
      put(file,"  1. Algorithmic evaluation of derivatives  : ");
    end if;
    put_line(file," yes");
    put(file,"  2. Output during the iterations           : ");
    if output
     then put_line(file," yes");
     else put_line(file," no");
    end if;
    put(file,"  3. Maximal #iterations allowed per root   :  ");
    put(file,maxitr,1); new_line(file);
    put(file,"  4. Maximal #deflations allowed per root   :  ");
    put(file,maxdef,1); new_line(file);
    put(file,"  5. #decimal places in working precision   :  ");
    put(file,nbdgts,1); new_line(file);
    put(file,"  6. Tolerance for error on the root        :");
    put(file,tolerr,3); new_line(file);
    put(file,"  7. Tolerance on the residual              :");
    put(file,tolres,3); new_line(file);
    put(file,"  8. Tolerance to decide the numerical rank :");
    put(file,tolrnk,3); new_line(file);
  end Display_Parameters;

  procedure Determine_Parameters
              ( inprecision : in natural32; 
                symbolic,output : out boolean;
                maxitr,maxdef,nbdgts : out natural32;
                tolerr,tolres,tolrnk : out double_float ) is

  -- DESCRIPTION :
  --   Allows the user to change the default values for the parameters.
  --   For double double and quad double precision, the value for
  --   inprecision should respectively be 31 and 63.

    ans : character;

  begin
    Set_Default_Parameters
      (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    if inprecision > 15
     then nbdgts := inprecision;
    end if;
    loop
      Display_Parameters(standard_output,symbolic,output,
                         maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
      put("Type 1, 2, 3, 4, 5, 6, 7, or 8 to change a value, or 0 to exit : ");
      Ask_Alternative(ans,"012345678");
      case ans is
        when '1' => if symbolic then
                      put("Do you want algorithmic evaluation (y/n) ? ");
                      Ask_Yes_or_No(ans);
                      symbolic := not (ans = 'y');
                    else
                      put("Do you want symbolic evaluation (y/n) ? ");
                      Ask_Yes_or_No(ans);
                      symbolic := (ans = 'y');
                    end if;
        when '2' => put("Do you want output during the iterations ? (y/n) ");
                    Ask_Yes_or_No(ans);
                    output := (ans = 'y');
        when '3' => put("Give the maximal #iterations/root : "); get(maxitr);
        when '4' => put("Give the maximal #deflations/root : "); get(maxdef);
        when '5' => put("Give the number of decimal places : "); get(nbdgts);
        when '6' => put("Give the tolerance for the error on the root : ");
                    get(tolerr);
        when '7' => put("Give the tolerance on the residual : "); get(tolres);
        when '8' => put("Give the tolerance to decide the numerical rank : ");
                    get(tolrnk);
        when others => null;
      end case;
      exit when (ans = '0');
    end loop;
  end Determine_Parameters;

  procedure Deflate_Singularities
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    symbolic,output : boolean;
    maxitr,maxdef,nbdgts : natural32;
    tolerr,tolres,tolrnk : double_float;

  begin
    if vrblvl > 0 then
      put("-> in drivers_to_deflate_singularities.");
      put_line("Deflate_Singularities 1 ...");
    end if;
    Set_Default_Parameters
      (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    Standard_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (p,sols,maxitr,maxdef,tolerr,tolres,tolrnk);
  end Deflate_Singularities;

  procedure Deflate_Singularities
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    symbolic,output : boolean;
    maxitr,maxdef,nbdgts : natural32;
    tolerr,tolres,tolrnk : double_float;

  begin
    if vrblvl > 0 then
      put("-> in drivers_to_deflate_singularities.");
      put_line("Deflate_Singularities 2 ...");
    end if;
    Set_Default_Parameters
      (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    DoblDobl_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (p,sols,maxitr,maxdef,tolerr,tolres,tolrnk);
  end Deflate_Singularities;

  procedure Deflate_Singularities
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    symbolic,output : boolean;
    maxitr,maxdef,nbdgts : natural32;
    tolerr,tolres,tolrnk : double_float;

  begin
    if vrblvl > 0 then
      put("-> in drivers_to_deflate_singularities.");
      put_line("Deflate_Singularities 3 ...");
    end if;
    Set_Default_Parameters
      (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    QuadDobl_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (p,sols,maxitr,maxdef,tolerr,tolres,tolrnk);
  end Deflate_Singularities;

  procedure Deflate_Singularities
              ( file : in file_type; outfilename : in string;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    timer : Timing_Widget;
    symbolic,output : boolean;
    maxitr,maxdef,nbdgts : natural32;
    tolerr,tolres,tolrnk : double_float;

  begin
    if vrblvl > 0 then
      put("-> in drivers_to_deflate_singularities.");
      put_line("Deflate_Singularities 4 ...");
    end if;
    new_line;
    Determine_Parameters
      (15,symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    new_line;
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
    if nbdgts <= 15 then
      if symbolic then
        Standard_Deflation_Methods.Symbolic_Deflation_and_Clustering
          (file,outfilename,p,sols,
           output,maxitr,maxdef,tolerr,tolres,tolrnk);
      else
        Standard_Deflation_Methods.Algorithmic_Deflation_and_Clustering
          (file,outfilename,p,sols,
           output,maxitr,maxdef,tolerr,tolres,tolrnk);
        new_line(file);
        put_line(file,"THE SOLUTIONS after deflation :");
        put(file,Standard_Complex_Solutions.Length_Of(sols),
                 natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    else
      declare
        size : constant natural32 := Decimal_to_Size(nbdgts);
        mp : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);
        msols : Multprec_Complex_Solutions.Solution_List
              := Multprec_Complex_Solutions.Create(sols);
      begin
        mp := Convert(p);
        Set_Size(mp,size);
        Multprec_Complex_Solutions.Set_Size(msols,size);
        if symbolic then
          Multprec_Deflation_Methods.Symbolic_Deflation_and_Clustering
            (file,outfilename,mp,msols,
             output,maxitr,maxdef,size,tolerr,tolres,tolrnk);
        else
          Multprec_Deflation_Methods.Interactive_Algorithmic_Deflation
            (file,mp,size,sols,tolrnk);
        end if;
      end;
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Deflating Isolated Singularities");
  end Deflate_Singularities;

  procedure Deflate_Singularities
              ( file : in file_type; outfilename : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    timer : Timing_Widget;
    symbolic,output : boolean;
    maxitr,maxdef,nbdgts : natural32;
    tolerr,tolres,tolrnk : double_float;

  begin
    if vrblvl > 0 then
      put("-> in drivers_to_deflate_singularities.");
      put_line("Deflate_Singularities 5 ...");
    end if;
    new_line;
    Determine_Parameters
      (31,symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    new_line;
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
    if symbolic then
      DoblDobl_Deflation_Methods.Symbolic_Deflation_and_Clustering
        (file,outfilename,p,sols,
         output,maxitr,maxdef,tolerr,tolres,tolrnk);
    else
      DoblDobl_Deflation_Methods.Algorithmic_Deflation_and_Clustering
        (file,outfilename,p,sols,
         output,maxitr,maxdef,tolerr,tolres,tolrnk);
      new_line(file);
      put_line(file,"THE SOLUTIONS after deflation :");
      put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
               natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Deflating Isolated Singularities");
  end Deflate_Singularities;

  procedure Deflate_Singularities
              ( file : in file_type; outfilename : in string;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    timer : Timing_Widget;
    symbolic,output : boolean;
    maxitr,maxdef,nbdgts : natural32;
    tolerr,tolres,tolrnk : double_float;

  begin
    if vrblvl > 0 then
      put("-> in drivers_to_deflate_singularities.");
      put_line("Deflate_Singularities 6 ...");
    end if;
    new_line;
    Determine_Parameters
      (63,symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    new_line;
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
    if symbolic then
      QuadDobl_Deflation_Methods.Symbolic_Deflation_and_Clustering
        (file,outfilename,p,sols,
         output,maxitr,maxdef,tolerr,tolres,tolrnk);
    else
      QuadDobl_Deflation_Methods.Algorithmic_Deflation_and_Clustering
        (file,outfilename,p,sols,
         output,maxitr,maxdef,tolerr,tolres,tolrnk);
      new_line(file);
      put_line(file,"THE SOLUTIONS after deflation :");
      put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
               natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Deflating Isolated Singularities");
  end Deflate_Singularities;

end Drivers_to_Deflate_Singularities;
