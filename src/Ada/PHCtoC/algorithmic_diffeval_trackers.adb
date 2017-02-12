with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Random_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Random_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;

package body Algorithmic_DiffEval_Trackers is

  function Prompt_for_Verbose return integer32 is

  -- DESCRIPTION :
  --   Asks the user if extra output is wanted and returns 0 or 1.

    ans : character;

  begin
    new_line;
    put("Extra output before and after the computations ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then return 1;
     else return 0;
    end if;
  end Prompt_for_Verbose;

-- WRAPPERS TO THE C INTERFACE FUNCTIONS :

  procedure Standard_ADE_Newton
              ( verbose : in integer32; pars : in Parameters ) is

    return_of_call : integer32;
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function newton ( v,maxstep,npred : integer32;
                      ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                      mxit : integer32; minrnd : double_float;
                      mxitref : integer32; minrndref : double_float )
                    return integer32;
    pragma import(C, newton, "standard_adenewton_with_pars");

  begin
    return_of_call := newton(verbose,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end Standard_ADE_Newton;

  procedure DoblDobl_ADE_Newton
              ( verbose : in integer32; pars : in Parameters ) is

    return_of_call : integer32;
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function newton ( v,maxstep,npred : integer32;
                      ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                      mxit : integer32; minrnd : double_float;
                      mxitref : integer32; minrndref : double_float )
                    return integer32;
    pragma import(C, newton, "dobldobl_adenewton_with_pars");

  begin
    return_of_call := newton(verbose,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end DoblDobl_ADE_Newton;

  procedure QuadDobl_ADE_Newton
              ( verbose : in integer32; pars : in Parameters ) is

    return_of_call : integer32;
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function newton ( v,maxstep,npred : integer32;
                      ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                      mxit : integer32; minrnd : double_float;
                      mxitref : integer32; minrndref : double_float )
                    return integer32;
    pragma import(C, newton, "quaddobl_adenewton_with_pars");

  begin
    return_of_call := newton(verbose,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end QuadDobl_ADE_Newton;

  procedure Standard_ADE_Track_One
              ( verbose : in integer32; 
                gamma : in Standard_Complex_Numbers.Complex_Number;
                pars : in Parameters ) is

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function one_track
               ( v : integer32; r,i : double_float;
                 maxstep,npred : integer32;
                 ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                 mxit : integer32; minrnd : double_float;
                 mxitref : integer32; minrndref : double_float )
               return integer32;
    pragma import(C, one_track, "standard_adeonepath_with_pars");

  begin
    return_of_call := one_track(verbose,regam,imgam,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end Standard_ADE_Track_One;

  procedure DoblDobl_ADE_Track_One
              ( verbose : in integer32; 
                gamma : in Standard_Complex_Numbers.Complex_Number;
                pars : in Parameters ) is

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function one_track
               ( v : integer32; r,i : double_float;
                 maxstep,npred : integer32;
                 ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                 mxit : integer32; minrnd : double_float;
                 mxitref : integer32; minrndref : double_float )
               return integer32;
    pragma import(C, one_track, "dobldobl_adeonepath_with_pars");

  begin
    return_of_call := one_track(verbose,regam,imgam,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end DoblDobl_ADE_Track_One;

  procedure QuadDobl_ADE_Track_One
              ( verbose : in integer32; 
                gamma : in Standard_Complex_Numbers.Complex_Number;
                pars : in Parameters ) is

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function one_track
               ( v : integer32; r,i : double_float;
                 maxstep,npred : integer32;
                 ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                 mxit : integer32; minrnd : double_float;
                 mxitref : integer32; minrndref : double_float )
               return integer32;
    pragma import(C, one_track, "quaddobl_adeonepath_with_pars");

  begin
    return_of_call := one_track(verbose,regam,imgam,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end QuadDobl_ADE_Track_One;

  procedure Standard_ADE_Track_Many
              ( verbose : in integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                pars : in Parameters ) is

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function many_track
               ( v : integer32; r,i : double_float;
                 maxstep,npred : integer32;
                 ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                 mxit : integer32; minrnd : double_float;
                 mxitref : integer32; minrndref : double_float )
               return integer32;
    pragma import(C, many_track, "standard_ademanypaths_with_pars");

  begin
    return_of_call := many_track(verbose,regam,imgam,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end Standard_ADE_Track_Many;

  procedure DoblDobl_ADE_Track_Many
              ( verbose : in integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                pars : in Parameters ) is

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function many_track
               ( v : integer32; r,i : double_float;
                 maxstep,npred : integer32;
                 ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                 mxit : integer32; minrnd : double_float;
                 mxitref : integer32; minrndref : double_float )
               return integer32;
    pragma import(C, many_track, "dobldobl_ademanypaths_with_pars");

  begin
    return_of_call := many_track(verbose,regam,imgam,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end DoblDobl_ADE_Track_Many;

  procedure QuadDobl_ADE_Track_Many
              ( verbose : in integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                pars : in Parameters ) is

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);
   -- flatten the path parameters :
    max_step : constant integer32 := pars.max_step;
    n_predictor : constant integer32 := pars.n_predictor;
    step_increase : constant double_float := pars.step_increase;
    step_decrease : constant double_float := pars.step_decrease;
    max_delta_t : constant double_float := pars.max_delta_t;
    max_delta_t_end : constant double_float := pars.max_delta_t_end;
    min_delta_t : constant double_float := pars.min_delta_t;
    err_max_res : constant double_float := pars.err_max_res;
    err_max_delta_x : constant double_float := pars.err_max_delta_x;
    err_max_first_delta_x : constant double_float
      := pars.err_max_first_delta_x;
    max_it : constant integer32 := pars.max_it;
    err_min_round_off : constant double_float := pars.err_min_round_off;
    max_it_refine : constant integer32 := pars.max_it_refine;
    err_min_round_off_refine : constant double_float
      := pars.err_min_round_off_refine;

    function many_track
               ( v : integer32; r,i : double_float;
                 maxstep,npred : integer32;
                 ic,dc,mxdt,mdtnd,mndt,mxres,mxdltx,mx1st : double_float;
                 mxit : integer32; minrnd : double_float;
                 mxitref : integer32; minrndref : double_float )
               return integer32;
    pragma import(C, many_track, "quaddobl_ademanypaths_with_pars");

  begin
    return_of_call := many_track(verbose,regam,imgam,
      max_step,n_predictor,step_increase,step_decrease,
      max_delta_t,max_delta_t_end,min_delta_t,
      err_max_res,err_max_delta_x,err_max_first_delta_x,
      max_it,err_min_round_off,
      max_it_refine,err_min_round_off_refine);
  end QuadDobl_ADE_Track_Many;

-- INTERACTIVE DRIVERS :

  procedure Standard_Newton is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    outfile : file_type;
    timer : timing_widget;
    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(16);

  begin
    new_line;
    put_line("Reading a system with solutions ...");
    Standard_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,p'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,p.all);
   -- put_line("Initializing the systems container ...");
    Standard_PolySys_Container.Initialize(p.all);
   -- put_line("Initializing the solutions container ...");
    Standard_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    Standard_ADE_Newton(verbose,pars);
    tstop(timer);
    newtsols := Standard_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after Newton's method :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line(outfile);
    print_times(outfile,timer,"Newton's method in double precision");  
  end Standard_Newton;

  procedure DoblDobl_Newton is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    outfile : file_type;
    timer : Timing_Widget;
    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(32);

  begin
    new_line;
    put_line("Reading a system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,p'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,p.all);
   -- put_line("Initializing the systems container ...");
    DoblDobl_PolySys_Container.Initialize(p.all);
   -- put_line("Initializing the solutions container ...");
    DoblDobl_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    DoblDobl_ADE_Newton(verbose,pars);
    tstop(timer);
    newtsols := DoblDobl_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after Newton's method :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line(outfile);
    print_times(outfile,timer,
                "Newton's method in double double precision");  
  end DoblDobl_Newton;

  procedure QuadDobl_Newton is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    outfile : file_type;
    timer : Timing_Widget;
    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(64);

  begin
    new_line;
    put_line("Reading a system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,p'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,p.all);
   -- put_line("Initializing the systems container ...");
    QuadDobl_PolySys_Container.Initialize(p.all);
   -- put_line("Initializing the solutions container ...");
    QuadDobl_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    QuadDobl_ADE_Newton(verbose,pars);
    tstop(timer);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after Newton's method :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line(outfile);
    print_times(outfile,timer,"Newton's method in quad double precision");
  end QuadDobl_Newton;

  procedure Standard_Track_one_Path is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    outfile : file_type;
    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(16);

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,target'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,target.all);
    new_line(outfile);
    put(outfile,"gamma = "); put(outfile,gamma); new_line(outfile);
   -- put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
   -- put_line("Initializing the solutions container ...");
    Standard_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    Standard_ADE_Track_One(verbose,gamma,pars);
    newtsols := Standard_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after path tracking :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    print_times(outfile,timer,"Tracking one path in double precision");
  end Standard_Track_one_Path;

  procedure DoblDobl_Track_one_Path is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    outfile : file_type;
    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(32);

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,target'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,target.all);
    new_line(outfile);
    put(outfile,"gamma = "); put(outfile,gamma); new_line(outfile);
   -- put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
   -- put_line("Initializing the solutions container ...");
    DoblDobl_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    DoblDobl_ADE_Track_One(verbose,gamma,pars);
    tstop(timer);
    newtsols := DoblDobl_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after path tracking :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    print_times(outfile,timer,"Tracking one path in double double precision");
  end DoblDobl_Track_one_Path;

  procedure QuadDobl_Track_one_Path is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    outfile : file_type;
    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(64);

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,target'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,target.all);
    new_line(outfile);
    put(outfile,"gamma = "); put(outfile,gamma); new_line(outfile);
   -- put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
   -- put_line("Initializing the solutions container ...");
    QuadDobl_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    QuadDobl_ADE_Track_One(verbose,gamma,pars);
    tstop(timer);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after path tracking :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    print_times(outfile,timer,"Tracking one path in quad double precision");
  end QuadDobl_Track_one_Path;

  procedure Standard_Track_many_Paths is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    outfile : file_type;
    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(16);

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,target'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,target.all);
    new_line(outfile);
    put(outfile,"gamma = "); put(outfile,gamma); new_line(outfile);
   -- put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
   -- put_line("Initializing the solutions container ...");
    Standard_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    Standard_ADE_Track_Many(verbose,gamma,pars);
    tstop(timer);
    newtsols := Standard_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after path tracking :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
      new_line;
      print_times
        (standard_output,timer,"tracking many paths in double precision");
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line(outfile);
    print_times(outfile,timer,"tracking many paths in double precision");
    close(outfile);
  end Standard_Track_many_Paths;

  procedure DoblDobl_Track_many_Paths is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    outfile : file_type;
    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(32);

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,target'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,target.all);
    new_line(outfile);
    put(outfile,"gamma = "); put(outfile,gamma); new_line(outfile);
   -- put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
   -- put_line("Initializing the solutions container ...");
    DoblDobl_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    DoblDobl_ADE_Track_Many(verbose,gamma,pars);
    tstop(timer);
    newtsols := DoblDobl_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after path tracking :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
      new_line;
      print_times(standard_output,timer,
                  "tracking many paths in double double precision");
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line(outfile);
    print_times(outfile,timer,"tracking many paths in double double precision");
    close(outfile);
  end DoblDobl_Track_many_Paths;

  procedure QuadDobl_Track_many_Paths is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    outfile : file_type;
    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    pars : Path_Parameters.Parameters
         := Path_Parameters.Default_Parameters(64);

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,target'last,1); put(outfile,"  ");
    put(outfile,Head_Of(sols).n,1); new_line(outfile);
    put(outfile,target.all);
    new_line(outfile);
    put(outfile,"gamma = "); put(outfile,gamma); new_line(outfile);
   -- put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
   -- put_line("Initializing the solutions container ...");
    QuadDobl_Solutions_Container.Initialize(sols);
    Path_Parameters.Tune(pars);
    new_line(outfile);
    put_line(outfile,"The values of the path parameters :");
    Path_Parameters.Write(outfile,pars);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    QuadDobl_ADE_Track_Many(verbose,gamma,pars);
    tstop(timer);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    if verbose > 0 then
      put_line("The solutions after path tracking :");
      put(standard_output,Length_Of(newtsols),
          natural32(Head_Of(newtsols).n),newtsols);
      new_line;
      print_times
        (standard_output,timer,"tracking many paths in quad double precision");
    end if;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line(outfile);
    print_times(outfile,timer,"tracking many paths in quad double precision");
    close(outfile);
  end QuadDobl_Track_many_Paths;

end Algorithmic_DiffEval_Trackers;
