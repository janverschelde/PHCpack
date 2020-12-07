with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_to_Real_Poly;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_to_Real_Poly;
with DoblDobl_Complex_Poly_Strings;
with DoblDobl_Complex_Laur_Strings;
with QuadDobl_Complex_to_Real_Poly;
with QuadDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Laur_Strings;
with String_System_Readers;
with Symbol_Table;
with Standard_Complex_Laur_Strings;
with Standard_Homotopy;
with Standard_Laurent_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Laurent_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Laurent_Homotopy;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Projective_Transformations;         use Projective_Transformations;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with DoblDobl_Root_Refiners;             use DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;             use QuadDobl_Root_Refiners;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Standard_Parameter_Systems;         use Standard_Parameter_Systems;
with DoblDobl_Parameter_Systems;         use DoblDobl_Parameter_Systems;
with QuadDobl_Parameter_Systems;         use QuadDobl_Parameter_Systems;
with Parameter_Homotopy_Continuation;    use Parameter_Homotopy_Continuation;
with Multitasking_Continuation;
with Write_Seed_Number;
with Greeting_Banners;

package body Main_Homotopy_Continuation is

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                target : in Complex_Number;
                sols,refsols : in out Standard_Complex_Solutions.Solution_List;
                solsfile : in boolean; nbequ : in natural32 := 0 ) is

    epsxa,epsfa,tolsing : constant double_float := 1.0E-8;
    nb : natural32 := 0;
    deflate : boolean := false;

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

  begin
    if Head_Of(sols).n > p'last
     then Affine_Transformation(sols);
    end if;
    if target = Create(1.0) then
      if solsfile then
        if nbequ = 0 then
          Reporting_Root_Refiner
            (outft,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
        else
          Reporting_Root_Sharpener
            (outft,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
        end if;
      else
        if nbequ = 0 then
          Reporting_Root_Refiner
            (outft,p,sols,epsxa,epsfa,tolsing,nb,5,deflate,false);
        else
          Reporting_Root_Sharpener
            (outft,p,sols,epsxa,epsfa,tolsing,nb,5,deflate,false);
        end if;
      end if;
    else
      declare
        pt : Poly_Sys(p'range);
      begin
        pt := Standard_Homotopy.Eval(target);
        if solsfile then
          if nbequ = 0 then
            Reporting_Root_Refiner
              (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
          else
            Reporting_Root_Sharpener
              (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
          end if;
        else
          if nbequ = 0 then
            Reporting_Root_Refiner
              (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,deflate,false);
          else
            Reporting_Root_Sharpener
              (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,deflate,false);
          end if;
        end if;
        Clear(pt);
      end;
    end if;
  end Refine_Solutions;

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                target : in Complex_Number;
                sols,refsols : in out DoblDobl_Complex_Solutions.Solution_List;
                solsfile : in boolean ) is

    epsxa,epsfa,tolsing : constant double_float := 1.0E-12;
    nb : natural32 := 0;
    deflate : boolean := false;
    ddtre : constant double_double := create(REAL_PART(target));
    ddtim : constant double_double := create(IMAG_PART(target));
    ddtarget : constant DoblDobl_Complex_Numbers.Complex_Number
             := DoblDobl_Complex_Numbers.Create(ddtre,ddtim);

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

  begin
    if Head_Of(sols).n > p'last
     then Affine_Transformation(sols);
    end if;
    if target = Create(1.0) then
      if solsfile then
        Reporting_Root_Refiner
          (outft,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
      else
        Reporting_Root_Refiner
          (outft,p,sols,epsxa,epsfa,tolsing,nb,5,deflate,false);
      end if;
    else
      declare
        pt : Poly_Sys(p'range);
      begin
        pt := DoblDobl_Homotopy.Eval(ddtarget);
        if solsfile then
          Reporting_Root_Refiner
            (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
        else
          Reporting_Root_Refiner
            (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,deflate,false);
        end if;
        Clear(pt);
      end;
    end if;
  end Refine_Solutions;

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                target : in Complex_Number;
                sols,refsols : in out QuadDobl_Complex_Solutions.Solution_List;
                solsfile : in boolean ) is

    epsxa,epsfa,tolsing : constant double_float := 1.0E-16;
    nb : natural32 := 0;
    deflate : boolean := false;
    qdtre : constant quad_double := create(REAL_PART(target));
    qdtim : constant quad_double := create(IMAG_PART(target));
    qdtarget : constant QuadDobl_Complex_Numbers.Complex_Number
             := QuadDobl_Complex_Numbers.Create(qdtre,qdtim);

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

  begin
    if Head_Of(sols).n > p'last
     then Affine_Transformation(sols);
    end if;
    if target = Create(1.0) then
      if solsfile then
        Reporting_Root_Refiner
          (outft,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
      else
        Reporting_Root_Refiner
          (outft,p,sols,epsxa,epsfa,tolsing,nb,5,deflate,false);
      end if;
    else
      declare
        pt : Poly_Sys(p'range);
      begin
        pt := QuadDobl_Homotopy.Eval(qdtarget);
        if solsfile then
          Reporting_Root_Refiner
            (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,deflate,false);
        else
          Reporting_Root_Refiner
            (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,deflate,false);
        end if;
        Clear(pt);
      end;
    end if;
  end Refine_Solutions;

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                target : in Complex_Number;
                sols,refsols : in out Standard_Complex_Solutions.Solution_List;
               -- solsfile : in boolean;
                nbequ : in natural32 := 0 ) is

    epsxa,epsfa,tolsing : constant double_float := 1.0E-8;
    nb : natural32 := 0;
   -- deflate : boolean := false;

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

  begin
   -- if Head_Of(sols).n > p'last
   --  then Affine_Transformation(sols);
   -- end if;
    if target = Create(1.0) then
      --if solsfile
      -- then
      if nbequ = 0 then
        Reporting_Root_Refiner
          (outft,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
      else
        Reporting_Root_Sharpener
          (outft,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
      end if;
      -- else Reporting_Root_Refiner
      --        (outft,p,sols,epsxa,epsfa,tolsing,nb,5,false);
      -- end if;
    else
      declare
        pt : Laur_Sys(p'range);
      begin
        pt := Standard_Laurent_Homotopy.Eval(target);
       -- if solsfile 
       --  then
        if nbequ = 0 then
          Reporting_Root_Refiner
            (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
        else
          Reporting_Root_Sharpener
            (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
        end if;
       --  else Reporting_Root_Refiner
       --         (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,false);
       -- end if;
        Clear(pt);
      end;
    end if;
  end Refine_Solutions;

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                target : in Complex_Number;
                sols,refsols : in out DoblDobl_Complex_Solutions.Solution_List;
                solsfile : in boolean ) is

    epsxa,epsfa,tolsing : constant double_float := 1.0E-12;
    nb : natural32 := 0;
    ddtre : constant double_double := create(REAL_PART(target));
    ddtim : constant double_double := create(IMAG_PART(target));
    ddtarget : constant DoblDobl_Complex_Numbers.Complex_Number
             := DoblDobl_Complex_Numbers.Create(ddtre,ddtim);

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

  begin
    if Head_Of(sols).n > p'last
     then Affine_Transformation(sols);
    end if;
    if target = Create(1.0) then
      if solsfile then
        Reporting_Root_Refiner
          (outft,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
      else
        Reporting_Root_Refiner
          (outft,p,sols,epsxa,epsfa,tolsing,nb,5,false);
      end if;
    else
      declare
        pt : Laur_Sys(p'range);
      begin
        pt := DoblDobl_Laurent_Homotopy.Eval(ddtarget);
        if solsfile then
          Reporting_Root_Refiner
            (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
        else
          Reporting_Root_Refiner
            (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,false);
        end if;
        Clear(pt);
      end;
    end if;
  end Refine_Solutions;

  procedure Refine_Solutions
              ( outft : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                target : in Complex_Number;
                sols,refsols : in out QuadDobl_Complex_Solutions.Solution_List;
                solsfile : in boolean ) is

    epsxa,epsfa,tolsing : constant double_float := 1.0E-12;
    nb : natural32 := 0;
    qdtre : constant quad_double := create(REAL_PART(target));
    qdtim : constant quad_double := create(IMAG_PART(target));
    qdtarget : constant QuadDobl_Complex_Numbers.Complex_Number
             := QuadDobl_Complex_Numbers.Create(qdtre,qdtim);

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

  begin
    if Head_Of(sols).n > p'last
     then Affine_Transformation(sols);
    end if;
    if target = Create(1.0) then
      if solsfile then
        Reporting_Root_Refiner
          (outft,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
      else
        Reporting_Root_Refiner
          (outft,p,sols,epsxa,epsfa,tolsing,nb,5,false);
      end if;
    else
      declare
        pt : Laur_Sys(p'range);
      begin
        pt := QuadDobl_Laurent_Homotopy.Eval(qdtarget);
        if solsfile then
          Reporting_Root_Refiner
            (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
        else
          Reporting_Root_Refiner
            (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,false);
        end if;
        Clear(pt);
      end;
    end if;
  end Refine_Solutions;

  procedure Write_Solutions
             ( solsft : in file_type;
               stsols : in Standard_Complex_Solutions.Solution_List;
               ddsols : in DoblDobl_Complex_Solutions.Solution_List;
               qdsols : in QuadDobl_Complex_Solutions.Solution_List;
               mpsols : in Multprec_Complex_Solutions.Solution_List ) is

    len : natural32;

  begin
    len := Standard_Complex_Solutions.Length_Of(stsols);
    if len > 0 then
      put(solsft,len,
          natural32(Standard_Complex_Solutions.Head_Of(stsols).n),stsols);
    else
      len := DoblDobl_Complex_Solutions.Length_Of(ddsols);
      if len > 0 then
        put(solsft,len,
            natural32(DoblDobl_Complex_Solutions.Head_Of(ddsols).n),ddsols);
      else
        len := QuadDobl_Complex_Solutions.Length_Of(qdsols);
        if len > 0 then
          put(solsft,len,
              natural32(QuadDobl_Complex_Solutions.Head_Of(qdsols).n),
              qdsols);
        else
          len := Multprec_Complex_Solutions.Length_Of(mpsols);
          if len > 0 then
            put(solsft,len,
                natural32(Multprec_Complex_Solutions.Head_Of(mpsols).n),
                mpsols);
          end if;
        end if;
      end if;
    end if;
  end Write_Solutions;

  procedure Write_Conclusion
              ( file : in file_type;
                start_moment,ended_moment : in Ada.Calendar.Time ) is
  begin
    new_line(file);
    put(file,"PHC ran from ");
    Time_Stamps.Write_Time_Stamp(file,start_moment);
    put(file," till ");
    Time_Stamps.Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    Time_Stamps.Write_Elapsed_Time(file,start_moment,ended_moment);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end Write_Conclusion;

  procedure Secant_Polynomial_Homotopy
              ( outfilename : in string;
                nbequ,nbvar,prclvl : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    solsft,outft : file_type;
    sols,refsols : Solution_List;
    ddsols,ddrefsols : DoblDobl_Complex_Solutions.Solution_List;
    qdsols,qdrefsols : QuadDobl_Complex_Solutions.Solution_List;
    mpsols : Multprec_Complex_Solutions.Solution_List;
    solsfile : boolean;
    ans : character;
    target : Complex_Number;
    nv_p : natural32;
    dd_p : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qd_p : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 then
      put("-> in main_homotopy_continuation.");
      put_line("Secant_Polynomial_Homotopy ...");
    end if;
    Create_Output_File(outft,outfilename);
    if nbequ = nbvar
     then put(outft,nbequ,p);
     else put(outft,nbequ,nbvar,p);
    end if;
    new_line(outft);
    new_line;
    put("Do you want the solutions on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Reading the name of the file to write the solutions on.");
      Read_Name_and_Create_File(solsft);
      solsfile := true;
    else
      solsfile := false;
    end if;
    Driver_for_Polynomial_Continuation
      (outft,p,prclvl,ls,sols,ddsols,qdsols,mpsols,target);
    if Length_Of(sols) > 0 then
      if nbequ = nbvar
       then Refine_Solutions(outft,p,target,sols,refsols,solsfile);
       else Refine_Solutions(outft,p,target,sols,refsols,solsfile,nbequ);
      end if;
    elsif DoblDobl_Complex_Solutions.Length_Of(ddsols) > 0 then
      nv_p := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      dd_p := DoblDobl_Complex_Poly_Strings.Parse(nv_p,ls.all);
      Refine_Solutions(outft,dd_p,target,ddsols,ddrefsols,solsfile);
    elsif QuadDobl_Complex_Solutions.Length_Of(qdsols) > 0 then
      nv_p := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      qd_p := QuadDobl_Complex_Poly_Strings.Parse(nv_p,ls.all);
      Refine_Solutions(outft,qd_p,target,qdsols,qdrefsols,solsfile);
    end if;
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outft,start_moment,ended_moment);
    Close(outft);
    if solsfile then
      Write_Solutions(solsft,refsols,ddrefsols,qdrefsols,mpsols);
      Close(solsft);
    end if;
  end Secant_Polynomial_Homotopy;

  procedure Secant_Laurent_Homotopy
              ( outfilename : in string;
                nbequ,nbvar,prclvl : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    solsft,outft : file_type;
    sols,refsols : Solution_List;
    ddsols,ddrefsols : DoblDobl_Complex_Solutions.Solution_List;
    qdsols,qdrefsols : QuadDobl_Complex_Solutions.Solution_List;
    mpsols : Multprec_Complex_Solutions.Solution_List;
    solsfile : boolean;
    ans : character;
    target : Complex_Number;
    nv_p : natural32;
    dd_p : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qd_p : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 then
      put("-> in main_homotopy_continuation.");
      put_line("Secant_Laurent_Homotopy ...");
    end if;
    Create_Output_File(outft,outfilename);
    if nbequ = nbvar
     then put(outft,nbequ,p);
     else put(outft,nbequ,nbvar,p);
    end if;
    new_line;
    put("Do you want the solutions on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Reading the name of the file to write the solutions on.");
      Read_Name_and_Create_File(solsft);
      solsfile := true;
    else
      solsfile := false;
    end if;
    Driver_for_Laurent_Continuation
      (outft,p,prclvl,ls,sols,ddsols,qdsols,mpsols,target);
    if Length_Of(sols) > 0 then
      if nbequ = nbvar
       then Refine_Solutions(outft,p,target,sols,refsols); -- ,solsfile);
       else Refine_Solutions(outft,p,target,sols,refsols,  --  solsfile,
                             nbequ);
      end if;
    elsif DoblDobl_Complex_Solutions.Length_Of(ddsols) > 0 then
      nv_p := Standard_Complex_Laurentials.Number_of_Unknowns(p(p'first));
      dd_p := DoblDobl_Complex_Laur_Strings.Parse(nv_p,ls.all);
      Refine_Solutions(outft,dd_p,target,ddsols,ddrefsols,solsfile);
    elsif QuadDobl_Complex_Solutions.Length_Of(qdsols) > 0 then
      nv_p := Standard_Complex_Laurentials.Number_of_Unknowns(p(p'first));
      qd_p := QuadDobl_Complex_Laur_Strings.Parse(nv_p,ls.all);
      Refine_Solutions(outft,qd_p,target,qdsols,qdrefsols,solsfile);
    end if;
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outft,start_moment,ended_moment);
    Close(outft);
    if solsfile then
      Write_Solutions(solsft,refsols,ddrefsols,qdrefsols,mpsols);
      Close(solsft);
    end if;
  end Secant_Laurent_Homotopy;

  procedure Parameter_Homotopy
              ( infile : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 ) is

    outfile : file_type;
    sols : Standard_Complex_Solutions.Solution_List;
    nb_equ,nb_unk,nb_par : integer32;

  begin
    if vrb > 0 
     then put_line("-> in mainpoco.Parameter_Homotopy 1 ...");
    end if;
    Read_Solution_Parameters(infile,outfile,p,sols,nb_equ,nb_unk,nb_par);
    Coefficient_Parameter_Homotopy_Continuation
      (outfile,p,sols,nb_equ,nb_unk,nb_par);
  end Parameter_Homotopy;

  procedure Parameter_Homotopy
              ( infile : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 ) is

    outfile : file_type;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    nb_equ,nb_unk,nb_par : integer32;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 
     then put_line("-> in mainpoco.Parameter_Homotopy 2 ...");
    end if;
    Read_Solution_Parameters(infile,outfile,p,sols,nb_equ,nb_unk,nb_par);
    Coefficient_Parameter_Homotopy_Continuation
      (outfile,p,sols,nb_equ,nb_unk,nb_par);
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outfile,start_moment,ended_moment);
    close(outfile);
  end Parameter_Homotopy;

  procedure Parameter_Homotopy
              ( infile : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 ) is

    outfile : file_type;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    nb_equ,nb_unk,nb_par : integer32;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 
     then put_line("-> in mainpoco.Parameter_Homotopy 3 ...");
    end if;
    Read_Solution_Parameters(infile,outfile,p,sols,nb_equ,nb_unk,nb_par);
    Coefficient_Parameter_Homotopy_Continuation
      (outfile,p,sols,nb_equ,nb_unk,nb_par);
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outfile,start_moment,ended_moment);
    close(outfile);
  end Parameter_Homotopy;

  procedure Sweep_Homotopy
              ( infile : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 ) is

    outfile : file_type;
    sols : Standard_Complex_Solutions.Solution_List;
    nb_equ,nb_unk,nb_par : integer32;
    isreal : boolean := Standard_Complex_to_Real_Poly.Is_Real(p);

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 
     then put_line("-> in mainpoco.Sweep_Homotopy 1 ...");
    end if;
    Read_Solution_Parameters(infile,outfile,p,sols,nb_equ,nb_unk,nb_par);
    Sweep(outfile,isreal,p,sols,nb_equ,nb_unk,nb_par);
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outfile,start_moment,ended_moment);
    close(outfile);
  end Sweep_Homotopy;

  procedure Sweep_Homotopy
              ( infile : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 ) is

    outfile : file_type;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    nb_equ,nb_unk,nb_par : integer32;
    isreal : boolean := DoblDobl_Complex_to_Real_Poly.Is_Real(p);

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 
     then put_line("-> in mainpoco.Sweep_Homotopy 2 ...");
    end if;
    Read_Solution_Parameters(infile,outfile,p,sols,nb_equ,nb_unk,nb_par);
    Sweep(outfile,isreal,p,sols,nb_equ,nb_unk,nb_par);
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outfile,start_moment,ended_moment);
    close(outfile);
  end Sweep_Homotopy;

  procedure Sweep_Homotopy
              ( infile : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 ) is

    outfile : file_type;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    nb_equ,nb_unk,nb_par : integer32;
    isreal : boolean := QuadDobl_Complex_to_Real_Poly.Is_Real(p);

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 
     then put_line("-> in mainpoco.Sweep_Homotopy 3 ...");
    end if;
    Read_Solution_Parameters(infile,outfile,p,sols,nb_equ,nb_unk,nb_par);
    Sweep(outfile,isreal,p,sols,nb_equ,nb_unk,nb_par);
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outfile,start_moment,ended_moment);
    close(outfile);
  end Sweep_Homotopy;

  procedure Multitasking_Secant_Homotopy
               ( outfilename : in string;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 ls : in Link_to_Array_of_Strings;
                 nt,nbequ,nbvar,prclvl : in natural32;
                 vrb : in integer32 := 0 ) is

    outft : file_type;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 then
      put("-> in main_homotopy_continuation.");
      put_line("Multitasking_Secant_Homotopy 1 ...");
    end if;
    Create_Output_File(outft,outfilename);
    Multitasking_Continuation.Driver_to_Path_Tracker
      (outft,p,prclvl,ls,integer32(nt),integer32(nbequ),integer32(nbvar));
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outft,start_moment,ended_moment);
    close(outft);
  end Multitasking_Secant_Homotopy;

  procedure Multitasking_Secant_Homotopy
               ( outfilename : in string;
                 p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 ls : in Link_to_Array_of_Strings;
                 nt,nbequ,nbvar,prclvl : in natural32;
                 vrb : in integer32 := 0 ) is

    outft : file_type;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0 then
      put("-> in main_homotopy_continuation.");
      put_line("Multitasking_Secant_Homotopy 2 ...");
    end if;
    Create_Output_File(outft,outfilename);
    Multitasking_Continuation.Driver_to_Path_Tracker
      (outft,p,prclvl,ls,integer32(nt),integer32(nbequ),integer32(nbvar));
    ended_moment := Ada.Calendar.Clock;
    Write_Conclusion(outft,start_moment,ended_moment);
    close(outft);
  end Multitasking_Secant_Homotopy;

  procedure Parameter_or_Sweep_Homotopy
              ( inft : in out file_type; prclvl : in natural32;
                lp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 ) is

    pos : constant character := Parameter_Homotopy_Continuation.Show_Menu;
    prc : character;

  begin
    if vrb > 0 
     then put_line("-> in mainpoco.Parameter_or_Sweep_Homotopy ...");
    end if;
    if prclvl = 1 then
      prc := Prompt_for_Precision;
    elsif prclvl = 2 then
      prc := '1';
    elsif prclvl = 4 then
      prc := '2';
    else
      prc := Prompt_for_Precision;
    end if;
    case prc is 
      when '0' =>
        if pos = '1'
         then Parameter_Homotopy(inft,lp.all,vrb-1);
         else Sweep_Homotopy(inft,lp.all,vrb-1);
        end if;
      when '1' =>
        declare
          nvr : constant natural32 
              := Standard_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
          ddp : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(lp'range)
              := DoblDobl_Complex_Poly_Strings.Parse(nvr,ls.all);
        begin
          if pos = '1'
           then Parameter_Homotopy(inft,ddp,vrb-1);
           else Sweep_Homotopy(inft,ddp,vrb-1);
          end if;
        end;
      when '2' =>
        declare
          nvr : constant natural32 
              := Standard_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
          qdp : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(lp'range)
              := QuadDobl_Complex_Poly_Strings.Parse(nvr,ls.all);
        begin
          if pos = '1'
           then Parameter_Homotopy(inft,qdp,vrb-1);
           else Sweep_Homotopy(inft,qdp,vrb-1);
          end if;
        end;
      when others => null;
    end case;
  end Parameter_or_Sweep_Homotopy;

  procedure Polynomial_Tracker
              ( outfilename : in string;
                inft : in out file_type; nt,prclvl : in natural32;
                lp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 ) is

    nva,neq : natural32;

  begin
    if vrb > 0 
     then put_line("-> in main_homotopy_continuation.Polynomial_Tracker ...");
    end if;
    neq := natural32(lp'last);
    nva := Standard_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    if nt = 0 then
      if nva <= neq then
        close(inft);
        Secant_Polynomial_Homotopy
          (outfilename,neq,nva,prclvl,lp.all,ls,vrb-1);
      else
        new_line;
        put("Found "); put(neq,1);
        put(" equations in "); put(nva,1); put_line(" unknowns...");
        new_line;
        Parameter_or_Sweep_Homotopy(inft,prclvl,lp,ls);
      end if;
    else
      Multitasking_Secant_Homotopy
        (outfilename,lp.all,ls,nt,neq,nva,prclvl,vrb-1);
    end if;
  end Polynomial_Tracker;

  procedure Standard_Laurent_Tracker
              ( outfilename : in string;
                inft : in out file_type; nt,prclvl : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 ) is

    nva,neq : natural32;

  begin
    if vrb > 0 then
      put("-> in main_homotopy_continuation.");
      put_line("Standard_Laurent_Tracker ...");
    end if;
    neq := natural32(p'last);
    nva := Number_of_Unknowns(p(p'first));
    if nt = 0 then
      if nva <= neq then
        close(inft);
        Secant_Laurent_Homotopy(outfilename,neq,nva,prclvl,p,ls,vrb-1);
      else
        new_line;
        put("Found "); put(neq,1);
        put(" equations in "); put(nva,1); put_line(" unknowns...");
        put_line("Laurent parameter homotopies are not yet supported ...");
       -- Parameter_or_Sweep_Homotopy(inft,prclvl,p,ls);
      end if;
    else
      Multitasking_Secant_Homotopy(outfilename,p,ls,nt,neq,nva,prclvl,vrb-1);
    end if;
  end Standard_Laurent_Tracker;

  procedure Main_Dispatch
              ( outfilename : in string;
                inft : in out file_type; nt,prclvl : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                ls : in Link_to_Array_of_Strings;
                vrb : in integer32 := 0 ) is
  begin
    if vrb > 0
     then put_line("-> in mainpoco.Main_Dispatch ...");
    end if;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(p) then
      Standard_Laurent_Tracker(outfilename,inft,nt,prclvl,p,ls,vrb-1);
    else
      declare
        use Standard_Complex_Poly_Systems;
        use Standard_Laur_Poly_Convertors;
        lp : Link_to_Poly_Sys;
      begin
        lp := new Poly_Sys'(Positive_Laurent_Polynomial_System(p));
        Polynomial_Tracker(outfilename,inft,nt,prclvl,lp,ls,vrb-1);
      end;
    end if;
  end Main_Dispatch;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   prclvl : in natural32; verbose : in integer32 := 0 ) is

    inft : file_type;
    ls : Link_to_Array_of_Strings;
    n,m : natural32;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in main_homotopy_continuation.Main ...");
    end if;
    String_System_Readers.Read_System(inft,infilename,n,m,ls);
    if ls = null then
      new_line;
      put_line("Reading the target polynomial system...");
      Read_Name_and_Open_File(inft);
      String_Splitters.get(inft,natural(n),natural(m),ls);
    end if;
   -- put("number of equations, n = "); put(n,1); new_line;
   -- put("number of variables, m = "); put(m,1); new_line;
   -- put_line("the strings : ");
   -- for i in ls'range loop
   --   put_line(ls(i).all);
   -- end loop;
   -- put_line("parsing the strings into a Laurent system ...");
    Symbol_Table.Init(m);
    declare
      q : constant Standard_Complex_Laur_Systems.Laur_Sys(1..integer32(n))
         := Standard_Complex_Laur_Strings.Parse(m,ls.all);
    begin
     -- put_line("the Laurent polynomial system : "); put(q);
      Main_Dispatch(outfilename,inft,nt,prclvl,q,ls,verbose-1);
    end;
  end Main;

end Main_Homotopy_Continuation;
