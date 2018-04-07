with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with Standard_Scaling;                   use Standard_Scaling;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Square_and_Embed_Systems;           use Square_and_Embed_Systems;
with Standard_Irreducible_Decomp;
with Standard_Irreducible_Decomp_io;     use Standard_Irreducible_Decomp_io;
with Multprec_Irreducible_Decomp;
with Multprec_Irreducible_Decomp_io;     use Multprec_Irreducible_Decomp_io;
with Homotopy_Cascade_Filter;            use Homotopy_Cascade_Filter;
with Drivers_to_Breakup_Components;      use Drivers_to_Breakup_Components;

package body Drivers_to_Witness_Generate is

  procedure Timing_Summary ( file : in file_type;
                             roco,hoco,poco,total : in duration ) is

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b00 : constant string :=
     "  =====================================================================";
    b1 : constant string :=
     "  |         TIMING INFORMATION SUMMARY for Black-Box Solver           |";
    b2 : constant string :=
     "  |   root counts  |  start system  |  continuation  |   total time   |";

  begin
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b00);
    put(file,"  | ");
    print_hms(file,roco); put(file," | ");
    print_hms(file,hoco); put(file," | ");
    print_hms(file,poco); put(file," | ");
    print_hms(file,total); put_line(file," |");
    put_line(file,b0);
  end Timing_Summary;

  procedure Black_Box_Solver
              ( file : in file_type;
                sys : in Standard_Complex_Poly_Systems.Poly_Sys;
                deg : in boolean;
                sols : out Standard_Complex_Solutions.Solution_List;
                rc : out natural32; totaltime : out duration ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    q : Poly_Sys(sys'range);
    roco,hoco,poco,total : duration;
    sols0 : Solution_List;

  begin
    tstart(timer);
    declare
      pp : Poly_Sys(sys'range);
    begin
      Copy(sys,pp);
      Black_Box_Root_Counting(file,0,pp,deg,rc,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        Scale(pp);
        Black_Box_Polynomial_Continuation(file,true,pp,q,sols,poco);
      end if;
      Clear(pp);
    end;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    new_line(file);
    print_times(file,timer,"Solving the polynomial system");
    new_line(file);
    Timing_Summary(file,roco,hoco,poco,total);
    totaltime := total;
  end Black_Box_Solver;

  procedure Write_Generate_Summary
              ( file : in file_type;
                npaths,nsols0,nsols1,ndiv : in Standard_Natural_Vectors.Vector;
                timings : in Array_of_Duration ) is

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b00 : constant string :=
     "  =====================================================================";
    b1 : constant string :=
     "  |       TIMING INFORMATION SUMMARY for Cascade of Homotopies        |";
    b2 : constant string :=
     "  | level | #paths | #isols | #comps | #infty |     user cpu time     |";

    total : duration := timings(timings'first);

  begin
    for i in timings'first+1..timings'last loop
      total := total + timings(i);
    end loop;
    new_line(file);
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b00);
    for i in reverse npaths'range loop
      put(file,"  | ");
      put(file,i,4); put(file,"  | ");
      put(file,npaths(i),5); put(file,"  | ");
      put(file,nsols1(i),5); put(file,"  | ");
      put(file,nsols0(i),5); put(file,"  | ");
      put(file,ndiv(i),5); put(file,"  |    ");
      print_hms(file,timings(integer(i))); put_line(file,"     |");
    end loop;
    put_line(file,b0);
    put(file,"  | total | ");
    put(file,Standard_Natural_Vectors.Sum(npaths),5);  put(file,"  | ");
    put(file,Standard_Natural_Vectors.Sum(nsols1),5);  put(file,"  | ");
    put(file,Standard_Natural_Vectors.Sum(nsols0),5);  put(file,"  | ");
    put(file,Standard_Natural_Vectors.Sum(ndiv),5); put(file,"  |    ");
    print_hms(file,total); put_line(file,"     |");
    put_line(file,b0);
  end Write_Generate_Summary;

  procedure Write ( file : in file_type; fp : in List; i : in integer32 ) is

    tmp : List := fp;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      put(file,lpt(i),4); put(file," |");
      tmp := Tail_Of(tmp);
    end loop;
  end Write;

  procedure Write_Banner
              ( file : in file_type; m : in natural32; sep : character ) is
  begin
    for i in 1..m loop
      for j in 1..6 loop
        put(file,sep);
      end loop;
    end loop;
  end Write_Banner;

  procedure Write_Classify_Summary
               ( file : in file_type; fp : in List;
                 timings : in Array_of_Duration ) is

    n : integer32;
    m : constant natural32 := Length_Of(fp);
    total : duration := timings(timings'first);

  begin
    for i in timings'first+1..timings'last loop
      total := total + timings(i);
    end loop;
    new_line(file);
    put_line(file,"  -------------------------------------------------------");
    put_line(file,"  |  TIMING INFORMATION SUMMARY for Classifying Points  |");
    put_line(file,"  -------------------------------------------------------");
    if m /= 0 then
      n := Head_Of(fp)'last;
      put(file,"  | dimension "); put(file," |");
      Write(file,fp,n); 
      put_line(file,"  user cpu time  |");
      put(file,"  =============="); Write_Banner(file,m,'=');
      put_line(file,"==================");
      if n > 0 then
        put(file,"  | level ");
        put(file,n-1,1);  put(file,"    |");
        Write(file,fp,n-1);
        put(file," "); print_hms(file,timings(integer(n)-1));
        put_line(file,"  |");
        for i in reverse 0..n-2 loop
          put(file,"  |       ");
          put(file,i,1);  put(file,"    |");
          Write(file,fp,i);
          put(file," "); print_hms(file,timings(integer(i)));
          put_line(file,"  |");
        end loop;
      end if;
      put(file,"  --------------"); Write_Banner(file,m,'-');
      put_line(file,"------------------");
    end if;
    put(file,"  | total time : ");
    Write_Banner(file,m,' ');
    print_hms(file,total); put_line(file,"  |");
    put(file,"  --------------"); Write_Banner(file,m,'-');
    put_line(file,"------------------");
  end Write_Classify_Summary;

  procedure Driver_for_Cascade_Filter
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                k : in integer32 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    nbequ : constant natural32 := natural32(p'length);
    nbunk : constant natural32 := Number_of_Unknowns(p(p'first));

  -- ROUTINES :

    function Max_Dim ( p : Poly_Sys ) return natural32 is

    -- DESCRIPTION :
    --   Returns either number of equations or the number of unknowns,
    --   depending on what is maximal.

    begin
      if nbequ >= nbunk
       then return nbequ;
       else return nbunk;
      end if;
    end Max_Dim;

    function Eliminate_Slices
                ( sqp : Standard_Complex_Poly_Systems.Poly_Sys;
                  m : natural32 )
                return Standard_Complex_Poly_Systems.Poly_Sys is

    -- DESCRIPTION :
    --   Given is a polynomial system sqp that is the result of making
    --   p square with nbunk - nbequ additional slices.  The system on
    --   return will have min(nbunk-nbequ-2,m) slices and variables less.
    --   The last two slices remain.

      use Standard_Complex_Polynomials;

      res : Poly_Sys(sqp'range);
      cnt : natural32 := nbunk - nbequ;
      reslast : integer32 := sqp'last;

    begin
      Copy(sqp,res);
      if cnt > 2 then
        for i in 1..m loop
          declare
            elim : Poly_Sys(res'first..reslast-1)
                 := Eliminate_Slice
                      (res,natural32(reslast-2),natural32(reslast));
          begin
            for i in elim'range loop
              Copy(elim(i),res(i));
            end loop;
            Clear(res(reslast));
            Clear(elim);
          end;
          cnt := cnt - 1;
          reslast := reslast - 1;
          exit when cnt <= 2;
        end loop;
      end if;
      return res(res'first..reslast);
    end Eliminate_Slices;

    procedure Get_Multprec_System 
              ( stsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                mpsys : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                size : in natural32 ) is

    -- DESCRIPTION :
    --   Depending on the answer of the user, the mpsys on return is read in
    --   from file or is the conversion of the given stsys.

      res : Multprec_Complex_Poly_Systems.Poly_Sys(stsys'range);
      ans : character;
      infile : file_type;
      m : natural32 := 0;
  
    begin
      put("Do you wish to read in the multi-precision coefficients? (y/n) "); 
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put_line("Reading the polynomial system"
                  & " with multi-precision coefficients.");
        Read_Name_and_Open_File(infile);
        get(infile,m,res);
        Set_Size(res,size);
      else
        res := Convert(stsys);
        Set_Size(res,size);
      end if;
      put_line(file,"The system in multi-precision format : ");
      put_line(file,res);
      mpsys := new Multprec_Complex_Poly_Systems.Poly_Sys'(res);
    end Get_Multprec_System;

    procedure Cascade_Filter
                ( stsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                  mpsys : in Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                  embp : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                  n,size,itp : natural32; skewproj,deg : in boolean ) is

    -- DESCRIPTION :
    --   Calls the routines to perform the filtering in the homotopy cascade.
  
    -- ON ENTRY :
    --   stsys    target system with standard complex coefficients;
    --   mpsys    system with multi-precision coefficients, null if size = 0;
    --   embp     sequence of embedded polynomial systems;
    --   n        original dimension;
    --   size     size of the numbers;
    --   itp      interpolator type;
    --   deg      if true, only degree based root counting, otherwise full;

      use Standard_Complex_Solutions;

      k : constant integer32 := embp'last;
      sols,sols0,sols1 : Solution_List;
      stansoco : Standard_Irreducible_Decomp.Solution_Components(0..k);
      multsoco : Multprec_Irreducible_Decomp.Solution_Components(0..k);
      tol : constant double_float := 10.0E-10;
      tol_sing : constant double_float := 10.0E-8;
      tol_eval : constant double_float := 10.0E-5;
      npa,ns0,ns1,div : Standard_Natural_Vectors.Vector(0..k) := (0..k => 0);
        -- npa = #paths,
        -- ns0 = #solutions with zz == 0
        -- ns0 = #solutions with zz /= 0
        -- dvi = #solutions diverged to infinity
      fp,fp_last : List;
      gentimes,clatimes : Array_of_Duration(0..integer(k))
                        := (0..integer(k) => 0.0);
      firstzero : boolean := false;  -- flag for first time zero added var

    begin
      if size = 0 
       then Standard_Initialize(stansoco,integer32(n),stsys);
       else Multprec_Initialize(multsoco,integer32(n),mpsys.all);
      end if;
      put(file,"EMBEDDED SYSTEM at the top dimension k = ");
      put(file,k,1); put_line(file," :");
      put_line(file,embp(k).all);
      Black_Box_Solver(file,embp(k).all,deg,sols,npa(k),gentimes(integer(k)));
      Filter_and_Split_Solutions(file,sols,integer32(n),k,tol,sols0,sols1);
      ns0(k) := Length_Of(sols0);
      ns1(k) := Length_Of(sols1);
      div(k) := npa(k) - ns0(k) - ns1(k);
      firstzero := (ns0(k) /= 0);
      if firstzero then
        if size = 0 then
          Standard_Update_Hypersurfaces
            (file,stansoco,n,natural32(k),natural32(k),itp,skewproj,
             embp(k).all,sols0,tol_sing,clatimes(integer(k)),fp,fp_last);
        else
          Multprec_Update_Hypersurfaces
            (file,multsoco,n,natural32(k),natural32(k),size,itp,skewproj,
             embp(k).all,mpsys.all,sols0,tol_sing,clatimes(integer(k)),
             fp,fp_last);
        end if;
      end if;
      new_line(file);
      if not Is_Null(sols1) then
        if k = 0 then
          if size = 0
           then Copy(sols1,stansoco(0).pts);
           else Copy(sols1,multsoco(0).pts);
          end if;
        else
          Copy(sols1,sols);
          if size = 0 then
            Standard_Cascade_Loop
              (file,n,natural32(k),itp,skewproj,embp,sols,stansoco,tol,tol_sing,
               tol_eval,npa,ns0,ns1,div,fp,fp_last,gentimes,clatimes);
          else
            Multprec_Cascade_Loop
              (file,n,natural32(k),size,itp,skewproj,embp,mpsys.all,sols,
               multsoco,tol,tol_sing,tol_eval,npa,ns0,ns1,div,
               fp,fp_last,gentimes,clatimes);
          end if;
        end if;
      end if;
      new_line(file);
      put_line(file,"THE SOLUTION COMPONENTS : ");
      if size = 0
       then put(file,stansoco);
       else put(file,multsoco);
      end if;
      Write_Generate_Summary(file,npa,ns0,ns1,div,gentimes);
      Write_Classify_Summary(file,fp,clatimes);
    end Cascade_Filter;

    procedure Main_Driver is

      n : natural32 := natural32(p'last);
      sqp : Poly_Sys(1..integer32(Max_Dim(p))) := Square(p);
      embp : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..k);
      stsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
      mpsys : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
      ans : character;
      itp,size : natural32;
      deg,skewproj : boolean;
      timer : Timing_Widget;

    begin
      tstart(timer);
      if nbequ < nbunk then
        --embp := Slice_and_Embed(sqp,k,1);
        Copy(sqp(sqp'last-1),sqp(sqp'last));
        if nbequ < nbunk - 1 then
          declare
            esqp : constant Poly_Sys := Eliminate_Slices(sqp,natural32(k));
          begin
            stsys := new Standard_Complex_Poly_Systems.Poly_Sys'(esqp);
          end;
        else
          stsys := new Standard_Complex_Poly_Systems.Poly_Sys'(sqp);
        end if;
       else stsys := new Standard_Complex_Poly_Systems.Poly_Sys'(sqp);
      end if;
      put_line(file,"The original system after squaring : ");
      put_line(file,stsys.all);
      embp := Slice_and_Embed(stsys.all,natural32(k));
      n := natural32(stsys'last);     -- WARNING : this could be bad....
      Breakup_Menu(itp,size,skewproj);
      if size > 0
       then Get_Multprec_System(stsys.all,mpsys,size);
      end if;
      put("Full root counting or only based on degrees? (f/d) ");
      Ask_Alternative(ans,"fd");
      deg := (ans = 'd');
      new_line;
      put_line("See the output file for results...");
      new_line;
      tstart(timer);
     -- put_line(file,"The original system (after squaring) : ");
     -- if nbequ > nbunk
     --  then Add_New_Symbols(k+nbequ-nbunk);
     --  else
      Add_Embed_Symbols(natural32(k));
     -- end if;
     -- put(file,stsys'last,stsys.all);
      Cascade_Filter(stsys.all,mpsys,embp,n,size,itp,skewproj,deg);
      tstop(timer);
      new_line(file);
      print_times(file,timer,"Numerical Irreducible Decomposition");
    end Main_Driver;

  begin
    Main_Driver;
  end Driver_for_Cascade_Filter;

end Drivers_to_Witness_Generate;
