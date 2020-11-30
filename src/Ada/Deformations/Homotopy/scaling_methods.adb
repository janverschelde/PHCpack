with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Multprec_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Scaling;
with DoblDobl_Scaling;
with QuadDobl_Scaling;
with Multprec_Scaling;

package body Scaling_Methods is

  procedure Display_Info is

    i : array(1..12) of string(1..65);

  begin
    i( 1):="By scaling the coefficients are transformed so that they  do  not";
    i( 2):="have extreme values.  The purpose is to avoid numerical problems.";
    i( 3):="  Equation scaling means that every polynomial is divided by  its";
    i( 4):="average coefficient.                                             ";
    i( 5):="  Variable scaling uses transformations like z  = (10^c)*x.   The";
    i( 6):="transformation  is  such  that  real  solutions remain real.  The";
    i( 7):="inverse  of  the  condition  number  of the linear system that is";
    i( 8):="solved to set up this transformation gives an indication  on  the";
    i( 9):="condition of the original polynomial system.                     ";
    i(10):="  Solution scaling transforms the solutions of  a  scaled  system";
    i(11):="back  into  the  original  coordinate  system.   Note that in the";
    i(12):="original coordinates, the solutions can be ill-conditioned.      ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Display_Info;

  procedure Equation_Scaling
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

    timer : Timing_Widget;

  begin
    put_line(file,"EQUATION SCALING :");
    tstart(timer);
    Standard_Scaling.Scale(p);
    tstop(timer);
    new_line(file); print_times(file,timer,"Equation Scaling"); new_line(file);
  end Equation_Scaling;

  procedure Equation_Scaling
              ( file : in file_type;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    timer : Timing_Widget;

  begin
    put_line(file,"EQUATION SCALING :");
    tstart(timer);
    DoblDobl_Scaling.Scale(p);
    tstop(timer);
    new_line(file); print_times(file,timer,"Equation Scaling"); new_line(file);
  end Equation_Scaling;

  procedure Equation_Scaling
              ( file : in file_type;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    timer : Timing_Widget;

  begin
    put_line(file,"EQUATION SCALING :");
    tstart(timer);
    QuadDobl_Scaling.Scale(p);
    tstop(timer);
    new_line(file); print_times(file,timer,"Equation Scaling"); new_line(file);
  end Equation_Scaling;

  procedure Equation_Scaling
              ( file : in file_type;
                p : in out Multprec_Complex_Poly_Systems.Poly_Sys ) is

    timer : Timing_Widget;

  begin
    put_line(file,"EQUATION SCALING :");
    tstart(timer);
    Multprec_Scaling.Scale(p);
    tstop(timer);
    new_line(file); print_times(file,timer,"Equation Scaling"); new_line(file);
  end Equation_Scaling;

  procedure Variable_Scaling
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                basis : out natural32;
                scvc : out Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Polynomials;

    yn : character;
    timer : Timing_Widget;
    rcond : double_float;
    bas : constant natural32 := 10;
    dim : constant integer32
        := p'length + integer32(Number_of_Unknowns(p(p'first)));
    scalecoeff : Standard_Complex_Vectors.Vector(1..dim);

  begin
    put_line(file,"EQUATION AND VARIABLE SCALING :");
    put("  Reducing the variability of coefficients ? (y/n) ");
    Ask_Yes_or_No(yn);
    tstart(timer);
    if yn = 'y' then
      put_line(file,"  Reducing the variability of coefficients.");
      Standard_Scaling.scale(p,bas,true,rcond,scalecoeff);
    else
      put_line(file,"  No reduction of variability of coefficients.");
      Standard_Scaling.scale(p,bas,false,rcond,scalecoeff);
    end if;
    tstop(timer);
    put("  The inverse condition is");
    put(rcond,3); put_line(".");
    put(file,"  The inverse condition is");
    put(file,rcond,3); put_line(file,".");
    basis := bas;
    scvc := new Standard_Complex_Vectors.Vector'(scalecoeff);
    new_line(file); print_times(file,timer,"Variable Scaling"); new_line(file);
  end Variable_Scaling;

  procedure Variable_Scaling
              ( file : in file_type;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                basis : out natural32;
                scvc : out DoblDobl_Complex_Vectors.Link_to_Vector ) is

    use DoblDobl_Complex_Polynomials;

    yn : character;
    timer : Timing_Widget;
    rcond : double_double;
    bas : constant natural32 := 10;
    dim : constant integer32
        := p'length + integer32(Number_of_Unknowns(p(p'first)));
    scalecoeff : DoblDobl_Complex_Vectors.Vector(1..dim);

  begin
    put_line(file,"EQUATION AND VARIABLE SCALING :");
    put("  Reducing the variability of coefficients ? (y/n) ");
    Ask_Yes_or_No(yn);
    tstart(timer);
    if yn = 'y' then
      put_line(file,"  Reducing the variability of coefficients.");
      DoblDobl_Scaling.scale(p,bas,true,rcond,scalecoeff);
    else
      put_line(file,"  No reduction of variability of coefficients.");
      DoblDobl_Scaling.scale(p,bas,false,rcond,scalecoeff);
    end if;
    tstop(timer);
    put("  The inverse condition is ");
    put(rcond,3); put_line(".");
    put(file,"  The inverse condition is ");
    put(file,rcond,3); put_line(file,".");
    basis := bas;
    scvc := new DoblDobl_Complex_Vectors.Vector'(scalecoeff);
    new_line(file); print_times(file,timer,"Variable Scaling"); new_line(file);
  end Variable_Scaling;

  procedure Variable_Scaling
              ( file : in file_type;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                basis : out natural32;
                scvc : out QuadDobl_Complex_Vectors.Link_to_Vector ) is

    use QuadDobl_Complex_Polynomials;

    yn : character;
    timer : Timing_Widget;
    rcond : quad_double;
    bas : constant natural32 := 10;
    dim : constant integer32
        := p'length + integer32(Number_of_Unknowns(p(p'first)));
    scalecoeff : QuadDobl_Complex_Vectors.Vector(1..dim);

  begin
    put_line(file,"EQUATION AND VARIABLE SCALING :");
    put("  Reducing the variability of coefficients ? (y/n) ");
    Ask_Yes_or_No(yn);
    tstart(timer);
    if yn = 'y' then
      put_line(file,"  Reducing the variability of coefficients.");
      QuadDobl_Scaling.scale(p,bas,true,rcond,scalecoeff);
    else
      put_line(file,"  No reduction of variability of coefficients.");
      QuadDobl_Scaling.scale(p,bas,false,rcond,scalecoeff);
    end if;
    tstop(timer);
    put("  The inverse condition is ");
    put(rcond,3); put_line(".");
    put(file,"  The inverse condition is ");
    put(file,rcond,3); put_line(file,".");
    basis := bas;
    scvc := new QuadDobl_Complex_Vectors.Vector'(scalecoeff);
    new_line(file); print_times(file,timer,"Variable Scaling"); new_line(file);
  end Variable_Scaling;

  procedure Variable_Scaling
              ( file : in file_type;
                p : in out Multprec_Complex_Poly_Systems.Poly_Sys;
                basis : out natural32;
                scvc : out Multprec_Complex_Vectors.Link_to_Vector ) is

    use Multprec_Complex_Polynomials;

    yn : character;
    timer : Timing_Widget;
    rcond : Floating_Number;
    bas : constant natural32 := 10;
    dim : constant integer32
        := p'length + integer32(Number_of_Unknowns(p(p'first)));
    scalecoeff : Multprec_Complex_Vectors.Vector(1..dim);

  begin
    put_line(file,"EQUATION AND VARIABLE SCALING :");
    put("  Reducing the variability of coefficients ? (y/n) ");
    Ask_Yes_or_No(yn);
    tstart(timer);
    if yn = 'y' then
      put_line(file,"  Reducing the variability of coefficients.");
      Multprec_Scaling.scale(p,bas,true,rcond,scalecoeff);
    else
      put_line(file,"  No reduction of variability of coefficients.");
      Multprec_Scaling.scale(p,bas,false,rcond,scalecoeff);
    end if;
    tstop(timer);
    put("  The inverse condition is ");
    put(rcond,3); put_line(".");
    put(file,"  The inverse condition is ");
    put(file,rcond,3); put_line(file,".");
    basis := bas;
    scvc := new Multprec_Complex_Vectors.Vector'(scalecoeff);
    new_line(file); print_times(file,timer,"Variable Scaling"); new_line(file);
  end Variable_Scaling;

  procedure Write_Results
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scvc : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Polynomials;

  begin
    new_line(file);
    put_line(file,"THE SCALED SYSTEM :");
    new_line(file);
    put(file,p'length,Number_of_Unknowns(p(p'first)),p);
    new_line(file);
    if basis /= 0 then
      new_line(file);
      put_line(file,"SCALING COEFFICIENTS :");
      new_line(file);
      put(file,basis,1); new_line(file);
      put_line(file,scvc);
    end if;
  end Write_Results;

  procedure Write_Results
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scvc : in DoblDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    new_line(file);
    put_line(file,"THE SCALED SYSTEM :");
    new_line(file);
    put(file,p);
    new_line(file);
    if basis /= 0 then
      new_line(file);
      put_line(file,"SCALING COEFFICIENTS :");
      new_line(file);
      put(file,basis,1); new_line(file);
      put_line(file,scvc);
    end if;
  end Write_Results;

  procedure Write_Results
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scvc : in QuadDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    new_line(file);
    put_line(file,"THE SCALED SYSTEM :");
    new_line(file);
    put(file,p);
    new_line(file);
    if basis /= 0 then
      new_line(file);
      put_line(file,"SCALING COEFFICIENTS :");
      new_line(file);
      put(file,basis,1); new_line(file);
      put_line(file,scvc);
    end if;
  end Write_Results;

  procedure Write_Results
              ( file : in file_type;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scvc : in Multprec_Complex_Vectors.Link_to_Vector ) is
  begin
    new_line(file);
    put_line(file,"THE SCALED SYSTEM :");
    new_line(file);
    put(file,p);
    new_line(file);
    if basis /= 0 then
      new_line(file);
      put_line(file,"SCALING COEFFICIENTS :");
      new_line(file);
      put(file,basis,1); new_line(file);
      put_line(file,scvc);
    end if;
  end Write_Results;

  procedure Main ( file : in file_type;
                   p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                   basis : out natural32;
                   scvc : out Standard_Complex_Vectors.Link_to_Vector;
                   verbose : in integer32 := 0 ) is

    ans : character;
    res_scvc : Standard_Complex_Vectors.Link_to_Vector;
    bas : natural32 := 0;

  begin
    if verbose > 0
     then put_line("-> in scaling_methods.Main 1 ...");
    end if;
    loop
      new_line;
      put_line("MENU for Scaling Polynomial Systems :");
      put_line("  0 : No Scaling       : leave the menu                    ");
      put_line("  1 : Equation Scaling : divide by average coefficient     ");
      put_line("  2 : Variable Scaling : change of variables, z = (10^c)*x ");
      put("Type 0, 1, or 2 to select scaling, or i for info : ");
      Ask_Alternative(ans,"012i");
      if ans = 'i'
       then new_line; Display_Info;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(file,p);
      when '2' => Variable_Scaling(file,p,bas,res_scvc);
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(file,p,bas,res_scvc);
      when others    => null;
    end case;
    basis := bas; scvc := res_scvc;
  end Main;

  procedure Main ( file : in file_type;
                   p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                   basis : out natural32;
                   scvc : out DoblDobl_Complex_Vectors.Link_to_Vector;
                   verbose : in integer32 := 0 ) is

    ans : character;
    res_scvc : DoblDobl_Complex_Vectors.Link_to_Vector;
    bas : natural32 := 0;

  begin
    if verbose > 0
     then put_line("-> in scaling_methods.Main 2 ...");
    end if;
    loop
      new_line;
      put_line("MENU for Scaling Polynomial Systems :");
      put_line("  0 : No Scaling       : leave the menu                    ");
      put_line("  1 : Equation Scaling : divide by average coefficient     ");
      put_line("  2 : Variable Scaling : change of variables, z = (10^c)*x ");
      put("Type 0, 1, or 2 to select scaling, or i for info : ");
      Ask_Alternative(ans,"012i");
      if ans = 'i'
       then new_line; Display_Info;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(file,p);
      when '2' => Variable_Scaling(file,p,bas,res_scvc);
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(file,p,bas,res_scvc);
      when others    => null;
    end case;
    basis := bas; scvc := res_scvc;
  end Main;

  procedure Main ( file : in file_type;
                   p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                   basis : out natural32;
                   scvc : out QuadDobl_Complex_Vectors.Link_to_Vector;
                   verbose : in integer32 := 0 ) is

    ans : character;
    res_scvc : QuadDobl_Complex_Vectors.Link_to_Vector;
    bas : natural32 := 0;

  begin
    if verbose > 0
     then put_line("-> in scaling_methods.Main 3 ...");
    end if;
    loop
      new_line;
      put_line("MENU for Scaling Polynomial Systems :");
      put_line("  0 : No Scaling       : leave the menu                    ");
      put_line("  1 : Equation Scaling : divide by average coefficient     ");
      put_line("  2 : Variable Scaling : change of variables, z = (10^c)*x ");
      put("Type 0, 1, or 2 to select scaling, or i for info : ");
      Ask_Alternative(ans,"012i");
      if ans = 'i'
       then new_line; Display_Info;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(file,p);
      when '2' => Variable_Scaling(file,p,bas,res_scvc);
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(file,p,bas,res_scvc);
      when others    => null;
    end case;
    basis := bas; scvc := res_scvc;
  end Main;

  procedure Main ( file : in file_type;
                   p : in out Multprec_Complex_Poly_Systems.Poly_Sys;
                   basis : out natural32;
                   scvc : out Multprec_Complex_Vectors.Link_to_Vector;
                   verbose : in integer32 := 0 ) is

    ans : character;
    res_scvc : Multprec_Complex_Vectors.Link_to_Vector;
    bas : natural32 := 0;

  begin
    if verbose > 0
     then put_line("-> in scaling_methods.Main 4 ...");
    end if;
    loop
      new_line;
      put_line("MENU for Scaling Polynomial Systems :");
      put_line("  0 : No Scaling       : leave the menu                    ");
      put_line("  1 : Equation Scaling : divide by average coefficient     ");
      put_line("  2 : Variable Scaling : change of variables, z = (10^c)*x ");
      put("Type 0, 1, or 2 to select scaling, or i for info : ");
      Ask_Alternative(ans,"012i");
      if ans = 'i'
       then new_line; Display_Info;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(file,p);
      when '2' => Variable_Scaling(file,p,bas,res_scvc);
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(file,p,bas,res_scvc);
      when others    => null;
    end case;
    basis := bas; scvc := res_scvc;
  end Main;

end Scaling_Methods;
