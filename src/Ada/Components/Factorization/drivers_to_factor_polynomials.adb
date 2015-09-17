with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Monodromy_Partitions;               use Monodromy_Partitions;
with Interpolate_Multivariate_Factor;    use Interpolate_Multivariate_Factor;
with Multivariate_Factorization;         use Multivariate_Factorization;

package body Drivers_to_Factor_Polynomials is

-- GENERATING THE INPUT :

  procedure Read_Polynomial ( n : out natural32; p : out Poly ) is

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

-- REPORTING THE OUTPUT :

  procedure Write_Factors ( file : in file_type; factors : in Poly_Sys;
                            mu : in Standard_Natural_Vectors.Vector ) is
  begin
    for i in factors'range loop
      new_line(file);
      put(file,"factor "); put(file,i,1);
      put_line(file," :"); put_line(file,factors(i));
      put(file,"with multiplicity = "); put(file,mu(i),1);
      new_line(file);
    end loop;
  end Write_Factors;

  function Maximal_Coefficient_Norm ( p : Poly ) return double_float is

  -- DESCRIPTION :
  --   Returns the maximal absolute value norm over all coefficients of p.

    max : double_float := 0.0;

    procedure Scan_Term ( t : Term; continue : out boolean ) is
    begin
      if AbsVal(t.cf) > max
       then max := AbsVal(t.cf);
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return max;
  end Maximal_Coefficient_Norm; 

  procedure Multiply_Factors ( p : in Poly; f : in Poly_Sys;
                               mu : in Standard_Natural_Vectors.Vector;
                               maxres : out double_float ) is

  -- DESCRIPTION :
  --   Multiplies the factors as many times as their multiplicities
  --   and compares with the given normalized polynomial p.

    mf : constant Poly := Multiply(f,mu);
    res : constant Poly := p - mf;

  begin
   -- put_line("The original polynomial :"); put_line(p);
   -- put_line("The multiplied factors :"); put_line(mf);
   -- put_line("The residual polynomial :"); put_line(res);
    maxres := Maximal_Coefficient_Norm(res);
  end Multiply_Factors; 

  procedure Multiply_Factors ( file : in file_type;
                               p : in Poly; f : in Poly_Sys;
                               mu : in Standard_Natural_Vectors.Vector;
                               maxres : out double_float ) is

  -- DESCRIPTION :
  --   Multiplies the factors as many times as their multiplicities
  --   and compares with the given polynomial p.

    mf : constant Poly := Multiply(f,mu);
    res : constant Poly := p - mf;

  begin
    put_line(file,"The original polynomial :"); put_line(file,p);
    put_line(file,"The multiplied factors :"); put_line(file,mf);
    put_line(file,"The residual polynomial :"); put_line(file,res);
    maxres := Maximal_Coefficient_Norm(res);
  end Multiply_Factors;

  procedure Write_Timing_Summary
               ( file : in file_type;
                 mongrp,lintrc,itrpol,mulval,total : in duration ) is

  -- DESCRIPTION :
  --   Writes the elapsed user times for monodromy groupings,
  --   linear trace certification, interpolation at the factors,
  --   and validation by multiplying the factors.

  begin
    put(file,"User time for finding the factorization   : ");
    print_hms(file,mongrp); new_line(file);
    put(file,"User time for linear traces certification : ");
    print_hms(file,lintrc); new_line(file);
    put(file,"User time for interpolation at factors    : ");
    print_hms(file,itrpol); new_line(file);
    put(file,"User time for multiplication validation   : ");
    print_hms(file,mulval); new_line(file);
    put(file,"Total elapsed user time for all stages    : ");
    print_hms(file,total); new_line(file);               
  end Write_Timing_Summary;

-- TARGET ROUTINES :

  procedure Factor ( monodromy : in boolean; n : in natural32; p : in Poly;
                     b,v : out Standard_Complex_Vectors.Vector;
                     wp : out Standard_Complex_Vectors.Link_to_Vector;
                     mw,mf : out Standard_Natural_Vectors.Link_to_Vector;
                     deco : out Standard_Natural_VecVecs.Link_to_VecVec;
                     fail : out boolean; factors : out Link_to_Poly_Sys ) is

    d : constant integer32 := Degree(p);
    tol : constant double_float := 1.0E-8;
    rad,dst : Standard_Floating_Vectors.Vector(1..d);
    rdp : Link_to_Poly_Sys;
    maxdif : double_float;

  begin
    if monodromy then
      Factor(p,n,natural32(d),deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
      if not fail then
        Certify(p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,maxdif);
        fail := (maxdif > tol);
      end if;
    else
      Trace_Factor(p,n,natural32(d),deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
    end if;
    if not fail
     then Interpolate(p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,factors);
    end if;
  end Factor;

  procedure Driver_to_Factor 
               ( file : in file_type; output,monodromy : in boolean;
                 n : in natural32; p : in Poly;
		 fail : out boolean; factors : out Link_to_Poly_Sys ) is

    d : constant integer32 := Degree(p);
    rad,dst : Standard_Floating_Vectors.Vector(1..d);
    timer,total_timer : Timing_Widget;
    b,v : Standard_Complex_Vectors.Vector(1..integer32(n));
    deco : Standard_Natural_VecVecs.Link_to_VecVec;
    wp : Standard_Complex_Vectors.Link_to_Vector;
    mf,mw : Standard_Natural_Vectors.Link_to_Vector;
    rdp : Link_to_Poly_Sys;
    maxdif,maxres : double_float;
    mongrp_time,lintrc_time,itrpol_time,mulval_time,total : duration;

  begin
    tstart(total_timer);
    tstart(timer);
    if monodromy then
      Factor(file,output,p,n,natural32(d),deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
    else
      Trace_Factor(file,p,n,natural32(d),deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
    end if;
    tstop(timer);
    if fail then
      put_line(file,"Failed to compute witness points.");
    else
      mongrp_time := Elapsed_User_Time(timer);
      put(file,"r = "); put(file,rad,3); new_line(file);
      put(file,"R = "); put(file,dst,3); new_line(file);
      new_line(file);
      print_times(file,timer,"finding the factorization");
      new_line(file);
      if monodromy then
        tstart(timer);
        Certify(file,p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,maxdif);
        tstop(timer);
        lintrc_time := Elapsed_User_Time(timer);
        put(file,"Maximal difference certificate : ");
        put(file,maxdif,3); put_line(file,".");
        new_line(file);
        print_times(file,timer,"linear traces certification");
        new_line(file);
      else
        lintrc_time := 0.0;
      end if;
      flush(file);
      tstart(timer);
      Interpolate(file,p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,factors);
      tstop(timer);
      itrpol_time := Elapsed_User_Time(timer);
      Write_Factors(file,factors.all,mf.all);
      new_line(file);
      print_times(file,timer,"interpolation of factors");
      new_line(file);
      tstart(timer);
      Multiply_Factors(file,p,factors.all,mf.all,maxres);
      tstop(timer);
      mulval_time := Elapsed_User_Time(timer);
      put(file,"Validation by multiplied factors :");
      put(file,maxres,3); put_line(file,".");
      new_line(file);
      print_times(file,timer,"multiplication of factors");
      new_line(file);
      tstop(total_timer);
      total := Elapsed_User_Time(total_timer);
      Write_Timing_Summary
        (file,mongrp_time,lintrc_time,itrpol_time,mulval_time,total);
    end if;
  end Driver_to_Factor;

  procedure Driver_to_Factor 
              ( monodromy : in boolean; n : in natural32; p : in Poly;
                fail : out boolean; factors : out Link_to_Poly_Sys ) is

    d : constant integer32 := Degree(p);
    rad,dst : Standard_Floating_Vectors.Vector(1..d);
    b,v : Standard_Complex_Vectors.Vector(1..integer32(n));
    deco : Standard_Natural_VecVecs.Link_to_VecVec;
    wp : Standard_Complex_Vectors.Link_to_Vector;
    mf,mw : Standard_Natural_Vectors.Link_to_Vector;
    rdp : Link_to_Poly_Sys;
    maxdif,maxres : double_float;
    timer,total_timer : Timing_Widget;
    mongrp_time,lintrc_time,itrpol_time,mulval_time,total : duration;

  begin
    tstart(total_timer);
    tstart(timer);
    if monodromy
     then Factor(p,n,natural32(d),deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
     else Trace_Factor(p,n,natural32(d),deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
    end if;
    tstop(timer);
    if fail then
      put_line("Failed to compute witness points.");
    else
      mongrp_time := Elapsed_User_Time(timer);
      put("r = "); put(rad,3); new_line;
      put("R = "); put(dst,3); new_line;
      put_line("The factorization, with multiplicities : ");
      Write_Factors(Standard_Output,deco.all,mw.all);
      put_line("The witness points : "); put_line(wp.all);
      put("with multiplicities : "); put(mw.all); new_line;
      if monodromy then
        tstart(timer);
        Certify(p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,maxdif);
        tstop(timer);
        lintrc_time := Elapsed_User_Time(timer);
        put("Maximal difference certificate : ");
        put(maxdif,3); put_line(".");
      else
        lintrc_time := 0.0;
      end if;
      flush(Standard_Output);
      tstart(timer);
      Interpolate(p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,factors);
      tstop(timer);
      itrpol_time := Elapsed_User_Time(timer);
      Write_Factors(Standard_Output,factors.all,mf.all);
      tstart(timer);
      Multiply_Factors(p,factors.all,mf.all,maxres);
      tstop(timer);
      mulval_time := Elapsed_User_Time(timer);
      put("Validation by multiplied factors :");
      put(maxres,3); put_line(".");
      tstop(total_timer);
      total := Elapsed_User_Time(total_timer);
      Write_Timing_Summary
        (Standard_Output,mongrp_time,lintrc_time,
         itrpol_time,mulval_time,total);
    end if;
  end Driver_to_Factor;

  procedure Factor ( n : in natural32; p : in Poly;
                     factors : out Link_to_Poly_Sys ) is

    ans : character;
    outfile : file_type;
    output,monodromy,fail : boolean;

  begin
    new_line;
    put("Do you wish intermediate output on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
      new_line;
      put_line("MENU for factorization method : ");
      put_line("  1. using monodromy to group witness points;");
      put_line("  2. combinatorially exploiting linear traces.");
      put("Type 1 or 2 to choose : ");
      Ask_Alternative(ans,"12");
      monodromy := (ans = '1');
      if monodromy then
        new_line;
        put("Do you wish intermediate output"
          & " during continuation ? (y/n) ");
        Ask_Yes_or_No(ans);
        output := (ans = 'y');
      end if;
      new_line;
      put_line("See the output file for results...");
      new_line;
      Driver_to_Factor(outfile,output,monodromy,n,p,fail,factors);
    else
      new_line;
      put_line("MENU for factorization method : ");
      put_line("  1. using monodromy to group witness points;");
      put_line("  2. combinatorially exploiting linear traces.");
      put("Type 1 or 2 to choose : ");
      Ask_Alternative(ans,"12");
      monodromy := (ans = '1');
      Driver_to_Factor(monodromy,n,p,fail,factors);
    end if;
  end Factor;

  procedure Driver_to_Factor_Polynomial is

    n : natural32 := 0;
    p : Poly;
    factors : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Numerical Factorization of Complex Multivariate Polynomials");
    Read_Polynomial(n,p);
    Normalize(p);
    Factor(n,p,factors);
  end Driver_to_Factor_Polynomial;

  procedure Write_Factors
               ( filename : in string; factors : in Poly_Sys ) is

  begin
    for i in factors'range loop
      declare
        suffix : constant string := Convert(i);
        name : constant string := filename & "_f" & suffix;
        file : file_type;
      begin
        Create(file,out_file,name);
        put(file,"1  ");
        put(file,Number_of_Unknowns(factors(i)),1);
        put_line(file,factors(i));
        Close(file);
      end;
    end loop;
  end Write_Factors;

  procedure Driver_to_Factor_Polynomial ( filename : in string ) is

    n : natural32 := 0;
    p : Poly;
    factors : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Numerical Factorization of Complex Multivariate Polynomials");
    Read_Polynomial(n,p);
    Normalize(p);
    Factor(n,p,factors);
    if factors /= null
     then Write_Factors(filename,factors.all);
    end if;
  end Driver_to_Factor_Polynomial;

end Drivers_to_Factor_Polynomials;
