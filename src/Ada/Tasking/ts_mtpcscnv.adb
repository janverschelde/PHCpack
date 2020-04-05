with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Homotopy;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy;
with QuadDobl_Homotopy_Convolutions_io;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Corrector_Convolutions;             use Corrector_Convolutions;
with Predictor_Corrector_Loops;          use Predictor_Corrector_Loops;
with Standard_Solutions_Queue;
with DoblDobl_Solutions_Queue;
with QuadDobl_Solutions_Queue;
with Multitasking;

procedure ts_mtpcscnv is

-- DESCRIPTION :
--   Development of multitasked tracking with predictor-corrector-shift
--   loops on homotopy systems of convolution circuits.

  procedure Allocate ( v : in out Standard_Integer_VecVecs.VecVec;
                       n : in integer32 ) is

  -- DESCRIPTION :
  --   Allocates vectors of range 1..n in v.

  begin
    for k in v'range loop
      v(k) := new Standard_Integer_Vectors.Vector(1..n);
    end loop;
  end Allocate;

  procedure Allocate ( v : in out Standard_Complex_VecVecs.VecVec;
                       n : in integer32 ) is

  -- DESCRIPTION :
  --   Allocates vectors of range 1..n in v.

  begin
    for k in v'range loop
      v(k) := new Standard_Complex_Vectors.Vector(1..n);
    end loop;
  end Allocate;

  procedure Allocate ( v : in out DoblDobl_Complex_VecVecs.VecVec;
                       n : in integer32 ) is

  -- DESCRIPTION :
  --   Allocates vectors of range 1..n in v.

  begin
    for k in v'range loop
      v(k) := new DoblDobl_Complex_Vectors.Vector(1..n);
    end loop;
  end Allocate;

  procedure Allocate ( v : in out QuadDobl_Complex_VecVecs.VecVec;
                       n : in integer32 ) is

  -- DESCRIPTION :
  --   Allocates vectors of range 1..n in v.

  begin
    for k in v'range loop
      v(k) := new QuadDobl_Complex_Vectors.Vector(1..n);
    end loop;
  end Allocate;

  procedure Standard_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Applies multitasking to track all paths in double precision.

  -- ON ENTRY :
  --   file     to write output information to;
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

    use Standard_Complex_Solutions,Standard_Speelpenning_Convolutions;
    use Standard_Predictor_Convolutions;

    maxit : constant integer32 := 4;
    nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail : boolean;
    homsa,abhsa : System_Array(1..nbtasks);
    homlead,abhlead : VecVecVec(1..nbtasks);
    homcff : VecVecVec_Array(1..nbtasks);
    dx : Standard_Complex_VecVecs.VecVec(1..nbtasks);
    wrk : Standard_Complex_VecVecs.VecVec(1..nbtasks);
    ipvt : Standard_Integer_VecVecs.VecVec(1..nbtasks);
    first : constant Link_to_Solution := Head_Of(sols);
    prd : Predictor_Array(1..nbtasks)
        := Create(nbtasks,first.v,hom.neq,hom.deg,
                  integer32(pars.numdeg),integer32(pars.dendeg),SVD);
    psv : Predictor_Vectors_Array(1..nbtasks)
        := Create(nbtasks,hom.dim,hom.neq);
    svh : SVD_Hessians_Array(1..nbtasks) := Create(nbtasks,hom.dim);

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t : double_float;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := 0.0;
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
                       maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
                       wrk(i),t,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.v := psv(i).sol; ls.t := Standard_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;
      t : double_float;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        cnt := Standard_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks  path "
                         & Multitasking.to_string(cnt));
        psv(i).sol := ls.v; t := 0.0;
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
                       maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
                       wrk(i),t,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.v := psv(i).sol; ls.t := Standard_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    for k in 1..nbtasks loop
      Copy(hom,homsa(k)); Copy(abh,abhsa(k));
      Allocate_Leading_Coefficients(hom.crc,homlead(k));
      Allocate_Leading_Coefficients(abh.crc,abhlead(k));
      Store_Leading_Coefficients(hom.crc,homlead(k));
      Store_Leading_Coefficients(abh.crc,abhlead(k));
      Allocate_Coefficients(hom.crc,homcff(k));
      Store_Coefficients(hom.crc,homcff(k));
    end loop;
    Allocate(ipvt,hom.dim); Allocate(dx,hom.dim);
    wrk := Allocate_Coefficients(nbtasks,hom.deg);
    Standard_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
    Clear(homsa); Clear(abhsa);
    Clear(homlead); Clear(abhlead); Clear(homcff);
    Standard_Complex_VecVecs.Clear(dx);
    Standard_Complex_VecVecs.Clear(wrk);
    Standard_Integer_VecVecs.Clear(ipvt);
    Clear(prd); Clear(psv); Clear(svh);
  end Standard_Multitasked_Tracker;

  procedure DoblDobl_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Applies multitasking to track all paths in double double precision.

  -- ON ENTRY :
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

    use DoblDobl_Complex_Solutions,DoblDobl_Speelpenning_Convolutions;
    use DoblDobl_Predictor_Convolutions;

    maxit : constant integer32 := 4;
    nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail : boolean;
    homsa,abhsa : System_Array(1..nbtasks);
    homlead,abhlead : VecVecVec(1..nbtasks);
    homcff : VecVecVec_Array(1..nbtasks);
    dx : DoblDobl_Complex_VecVecs.VecVec(1..nbtasks);
    wrk : DoblDobl_Complex_VecVecs.VecVec(1..nbtasks);
    ipvt : Standard_Integer_VecVecs.VecVec(1..nbtasks);
    first : constant Link_to_Solution := Head_Of(sols);
    prd : Predictor_Array(1..nbtasks)
        := Create(nbtasks,first.v,hom.neq,hom.deg,
                  integer32(pars.numdeg),integer32(pars.dendeg),SVD);
    psv : Predictor_Vectors_Array(1..nbtasks)
        := Create(nbtasks,hom.dim,hom.neq);
    svh : SVD_Hessians_Array(1..nbtasks) := Create(nbtasks,hom.dim);

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t : double_double;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
                       maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
                       wrk(i),t,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.v := psv(i).sol; ls.t := DoblDobl_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;
      t : double_double;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        cnt := DoblDobl_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks  path "
                         & Multitasking.to_string(cnt));
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
                       maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
                       wrk(i),t,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.v := psv(i).sol; ls.t := DoblDobl_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    for k in homsa'range loop
      Copy(hom,homsa(k)); Copy(hom,abhsa(k));
      Allocate_Leading_Coefficients(hom.crc,homlead(k));
      Allocate_Leading_Coefficients(abh.crc,abhlead(k));
      Store_Leading_Coefficients(hom.crc,homlead(k));
      Store_Leading_Coefficients(abh.crc,abhlead(k));
      Allocate_Coefficients(hom.crc,homcff(k));
      Store_Coefficients(hom.crc,homcff(k));
    end loop;
    Allocate(ipvt,hom.dim); Allocate(dx,hom.dim);
    wrk := Allocate_Coefficients(hom.dim,hom.deg);
    DoblDobl_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
    Clear(homsa); Clear(abhsa);
    Clear(homlead); Clear(abhlead); Clear(homcff);
    DoblDobl_Complex_VecVecs.Clear(dx);
    DoblDobl_Complex_VecVecs.Clear(wrk);
    Standard_Integer_VecVecs.Clear(ipvt);
    Clear(prd); Clear(psv); Clear(svh);
  end DoblDobl_Multitasked_Tracker;

  procedure QuadDobl_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Applies multitasking to track all paths in quad double precision.

  -- ON ENTRY :
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

    use QuadDobl_Complex_Solutions,QuadDobl_Speelpenning_Convolutions;
    use QuadDobl_Predictor_Convolutions;

    maxit : constant integer32 := 4;
    nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail : boolean;
    homsa,abhsa : System_Array(1..nbtasks);
    homlead,abhlead : VecVecVec(1..nbtasks);
    homcff : VecVecVec_Array(1..nbtasks);
    dx : QuadDobl_Complex_VecVecs.VecVec(1..nbtasks);
    wrk : QuadDobl_Complex_VecVecs.VecVec(1..nbtasks);
    ipvt : Standard_Integer_VecVecs.VecVec(1..nbtasks);
    first : constant Link_to_Solution := Head_Of(sols);
    prd : Predictor_Array(1..nbtasks)
        := Create(nbtasks,first.v,hom.neq,hom.deg,
                  integer32(pars.numdeg),integer32(pars.dendeg),SVD);
    psv : Predictor_Vectors_Array(1..nbtasks)
        := Create(nbtasks,hom.dim,hom.neq);
    svh : SVD_Hessians_Array(1..nbtasks) := Create(nbtasks,hom.dim);

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t : quad_double;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
                       maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
                       wrk(i),t,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.v := psv(i).sol; ls.t := QuadDobl_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;
      t : quad_double;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        cnt := QuadDobl_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks  path "
                         & Multitasking.to_string(cnt));
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
                       maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
                       wrk(i),t,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.v := psv(i).sol; ls.t := QuadDobl_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    for k in homsa'range loop
      Copy(hom,homsa(k)); Copy(hom,abhsa(k));
      Allocate_Leading_Coefficients(hom.crc,homlead(k));
      Allocate_Leading_Coefficients(abh.crc,abhlead(k));
      Store_Leading_Coefficients(hom.crc,homlead(k));
      Store_Leading_Coefficients(abh.crc,abhlead(k));
      Allocate_Coefficients(hom.crc,homcff(k));
      Store_Coefficients(hom.crc,homcff(k));
    end loop;
    Allocate(ipvt,hom.dim); Allocate(dx,hom.dim);
    wrk := Allocate_Coefficients(hom.dim,hom.deg);
    QuadDobl_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
    Clear(homsa); Clear(abhsa);
    Clear(homlead); Clear(abhlead); Clear(homcff);
    QuadDobl_Complex_VecVecs.Clear(dx);
    QuadDobl_Complex_VecVecs.Clear(wrk);
    Standard_Integer_VecVecs.Clear(ipvt);
    Clear(prd); Clear(psv); Clear(svh);
  end QuadDobl_Multitasked_Tracker;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double precision.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg,nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    ans : character;
    file : file_type;
    verbose : boolean;

  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    Standard_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    pars.gamma := Standard_Homotopy.Accessibility_Constant;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Give the number of tasks : "); get(nbt);
    Standard_Multitasked_Tracker(nbt,cnvhom,abshom,sols,pars,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(sols),
        natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double double precisin.

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg,nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    ddgamma : DoblDobl_Complex_Numbers.Complex_Number;
    ans : character;
    file : file_type;
    verbose : boolean;

    use DoblDobl_Complex_Numbers_cv;
  
  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    DoblDobl_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    ddgamma := DoblDobl_Homotopy.Accessibility_Constant;
    pars.gamma := DoblDobl_Complex_to_Standard(ddgamma);
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Give the number of tasks : "); get(nbt);
    DoblDobl_Multitasked_Tracker(nbt,cnvhom,abshom,sols,pars,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
        natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in quad double precision.

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg,nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    qdgamma : QuadDobl_Complex_Numbers.Complex_Number;
    ans : character;
    file : file_type;
    verbose : boolean;

    use QuadDobl_Complex_Numbers_cv;

  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    QuadDobl_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    qdgamma := QuadDobl_Homotopy.Accessibility_Constant;
    pars.gamma := QuadDobl_Complex_to_Standard(qdgamma);
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Give the number of tasks : "); get(nbt);
    QuadDobl_Multitasked_Tracker(nbt,cnvhom,abshom,sols,pars,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
        natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to the selected precision.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtpcscnv;
