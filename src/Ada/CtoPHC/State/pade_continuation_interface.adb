with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_cv;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Homotopy;
with Standard_Coefficient_Homotopy;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with DoblDobl_Homotopy;
with DoblDobl_Coefficient_Homotopy;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with QuadDobl_Homotopy;
with QuadDobl_Coefficient_Homotopy;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with Standard_Pade_Approximants;
with DoblDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Series_Path_Trackers;
with Drivers_to_Series_Trackers;
with Standard_SeriesPade_Tracker;
with DoblDobl_SeriesPade_Tracker;
with QuadDobl_SeriesPade_Tracker;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;

package body Pade_Continuation_Interface is

  function Pade_Continuation_Parameters_Set_Defaults
             ( vrblvl : integer32 := 0 ) return integer32 is

    pars : constant Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Parameters_Set_Defaults ...");
    end if;
    Homotopy_Continuation_Parameters.Construct(pars);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Set_Defaults");
      end if;
      return 735;
  end Pade_Continuation_Parameters_Set_Defaults;

  function Pade_Continuation_Parameters_Get_Value
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    idx : constant natural32 := natural32(v_a(v_a'first));
    fail : integer32 := 0;
    homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

    use Homotopy_Continuation_Parameters;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Parameters_Get_Value ...");
    end if;
    if homconpars = null then
      fail := Pade_Continuation_Parameters_Set_Defaults;
      homconpars := Homotopy_Continuation_Parameters.Retrieve;
    end if;
    case idx is
      when  1 => Assign(homconpars.gamma,c);
      when  2 => Assign(integer32(homconpars.numdeg),b);
      when  3 => Assign(integer32(homconpars.dendeg),b);
      when  4 => Assign(homconpars.maxsize,c);
      when  5 => Assign(homconpars.minsize,c);
      when  6 => Assign(homconpars.pbeta,c);
      when  7 => Assign(homconpars.cbeta,c);
      when  8 => Assign(homconpars.alpha,c);
      when  9 => Assign(homconpars.tolres,c);
      when 10 => Assign(homconpars.epsilon,c);
      when 11 => Assign(integer32(homconpars.corsteps),b);
      when 12 => Assign(integer32(homconpars.maxsteps),b);
      when others => 
        put_line("Index value for the parameter is out of range.");
        fail := 737;
    end case;
    return fail;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Get_Value");
      end if;
      return 737;
  end Pade_Continuation_Parameters_Get_Value;

  function Pade_Continuation_Parameters_Set_Value
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    fail : integer32 := 0;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    idx : constant natural32 := natural32(v_a(v_a'first));
    v_b : C_Integer_Array(0..0);
    v_c : C_Double_Array(0..0);
    v_gamma : C_Double_Array(0..1);
    regamma,imgamma : double_float;

    homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

    use Interfaces.C;
    use Homotopy_Continuation_Parameters;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Parameters_Set_Value ...");
    end if;
    if homconpars = null then
      fail := Pade_Continuation_Parameters_Set_Defaults;
      homconpars := Homotopy_Continuation_Parameters.Retrieve;
    end if;
    case idx is
      when  1 =>
        v_gamma := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2));
        regamma := double_float(v_gamma(v_gamma'first));
        imgamma := double_float(v_gamma(v_gamma'first+1));
        homconpars.gamma := Standard_Complex_Numbers.Create(regamma,imgamma);
      when  2 => v_b := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
                 homconpars.numdeg := natural32(v_b(v_b'first));
      when  3 => v_b := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
                 homconpars.dendeg := natural32(v_b(v_b'first));
      when  4 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.maxsize := double_float(v_c(v_c'first));
      when  5 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.minsize := double_float(v_c(v_c'first));
      when  6 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.pbeta := double_float(v_c(v_c'first));
      when  7 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.cbeta := double_float(v_c(v_c'first));
      when  8 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.alpha := double_float(v_c(v_c'first));
      when  9 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.tolres := double_float(v_c(v_c'first));
      when 10 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.epsilon := double_float(v_c(v_c'first));
      when 11 => v_b := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
                 homconpars.corsteps := natural32(v_b(v_b'first));
      when 12 => v_b := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
                 homconpars.maxsteps := natural32(v_b(v_b'first));
      when others =>
        put_line("Index value for the parameter is out of range.");
        fail := 738;
    end case;
    return fail;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Set_Value");
      end if;
      return 738;
  end Pade_Continuation_Parameters_Set_Value;

  function Pade_Continuation_Parameters_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    pars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
         := Homotopy_Continuation_Parameters.Retrieve;

  begin
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      Homotopy_Continuation_Parameters_io.put
        (PHCpack_Operations.output_file,pars.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      Homotopy_Continuation_Parameters_io.put(pars.all);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Parameters_Write");
      end if;
      return 874;
  end Pade_Continuation_Parameters_Write;

  function Pade_Continuation_Parameters_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Parameters_Clear ...");
    end if;
    Homotopy_Continuation_Parameters.Destruct;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Clear");
      end if;
      return 736;
  end Pade_Continuation_Parameters_Clear;

  function Pade_Continuation_Parameters_Reset_Values
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Parameters_Reset_Values ...");
    end if;
    case prc is
      when 0 => Standard_SeriesPade_Tracker.Init(homconpars.all);
      when 1 => DoblDobl_SeriesPade_Tracker.Init(homconpars.all);
      when 2 => QuadDobl_SeriesPade_Tracker.Init(homconpars.all);
      when others => null;
    end case;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Parameters_Reset_Values");
      end if;
      return 740;
  end Pade_Continuation_Parameters_Reset_Values;

  procedure Standard_Track
              ( name : in string; localfile,verbose : in boolean;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Tracks the solution paths in standard precision,
  --   writing no output if name is empty,
  --   otherwise, creates an output file with the given name.
  --   The file name is the defined output file if localfile is false,
  --   otherwise, if localfile, the the file is a local variable.
  --   The value for mhom is 0, 1, 2 or higher, depending whether
  --   affine, 1-homogeneous, or m-homogeneous coordinates are applied.

    file : file_type;
    start,target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;
    nvr : natural32;

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    PHCpack_Operations.Retrieve_Target_System(target);
    if mhom = 0 then
      Standard_Homotopy.Create(target.all,start.all,tpow,homconpars.gamma);
    else
      Standard_Homotopy.Create(target.all,start.all,1,homconpars.gamma);
      Standard_Coefficient_Homotopy.Create
        (start.all,target.all,1,homconpars.gamma);
    end if;
    if name = "" then
      if not verbose then
        Drivers_to_Series_Trackers.Standard_Track
          (target'last,sols,homconpars.all,mhom,idz); --,verbose);
      else
        Homotopy_Continuation_Parameters_io.put(homconpars.all);
        Drivers_to_Series_Trackers.Standard_Track
          (standard_output,target'last,sols,homconpars.all,mhom,idz,verbose);
      end if;
    elsif localfile then
      Create(file,out_file,name);
      put(file,natural32(target'last),target.all);
      new_line(file);
      put_line(file,"THE START SYSTEM :");
      put(file,natural32(start'last),start.all);
      new_line(file);
      put_line(file,"THE START SOLUTIONS :");
      put(file,Standard_Complex_Solutions.Length_Of(sols),
          natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
      new_line(file);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.Standard_Track
        (file,target'last,sols,homconpars.all,mhom,idz,verbose);
      close(file);
    else
      PHCpack_Operations.Define_Output_File(name);
      put(PHCpack_Operations.output_file,natural32(target'last),target.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,natural32(start'last),start.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SOLUTIONS :");
      put(PHCpack_Operations.output_file,
          Standard_Complex_Solutions.Length_Of(sols),
          natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
      new_line(PHCpack_Operations.output_file);
      Homotopy_Continuation_Parameters_io.put
        (PHCpack_Operations.output_file,homconpars.all);
      if mhom > 1 then
        nvr := natural32(start'last) - mhom;
        new_line(PHCpack_Operations.output_file);
        Series_Path_Trackers.Write_Partition
          (PHCpack_Operations.output_file,nvr,mhom,idz);
        new_line(PHCpack_Operations.output_file);
      end if;
      Drivers_to_Series_Trackers.Standard_Track
        (PHCpack_Operations.output_file,
         target'last,sols,homconpars.all,mhom,idz,verbose);
    end if;
   -- put_line("Clearing the solutions container ...");
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(sols);
  end Standard_Track;

  procedure DoblDobl_Track
              ( name : in string; localfile,verbose : in boolean;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Tracks the solution paths in standard precision,
  --   writing no output if name is empty,
  --   otherwise, creates an output file with the given name.
  --   The file name is the defined output file if localfile is false,
  --   otherwise, if localfile, the the file is a local variable.
  --   The value for mhom is 0, 1, 2 or higher, depending whether
  --   affine, 1-homogeneous, or m-homogeneous coordinates are applied.

    use DoblDobl_Complex_Numbers_cv;

    file : file_type;
    start,target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;
    nvr : natural32;

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

    gamma : constant Standard_Complex_Numbers.Complex_Number
          := homconpars.gamma;
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(gamma);

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    PHCpack_Operations.Retrieve_Target_System(target);
    if mhom = 0 then
      DoblDobl_Homotopy.Create(target.all,start.all,tpow,dd_gamma);
    else
      DoblDobl_Homotopy.Create(target.all,start.all,1,dd_gamma);
      DoblDobl_Coefficient_Homotopy.Create(start.all,target.all,1,dd_gamma);
    end if;
    if name = "" then
      if not verbose then
        Drivers_to_Series_Trackers.DoblDobl_Track
          (target'last,sols,homconpars.all,mhom,idz); -- ,verbose);
      else
        Homotopy_Continuation_Parameters_io.put(homconpars.all);
        Drivers_to_Series_Trackers.DoblDobl_Track
          (standard_output,target'last,sols,homconpars.all,mhom,idz,verbose);
      end if;
    elsif localfile then
      Create(file,out_file,name);
      put(file,natural32(target'last),target.all);
      new_line(file);
      put_line(file,"THE START SYSTEM :");
      put(file,natural32(start'last),start.all);
      new_line(file);
      put_line(file,"THE START SOLUTIONS :");
      put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
          natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
      new_line(file);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.DoblDobl_Track
        (file,target'last,sols,homconpars.all,mhom,idz,verbose);
      close(file);
    else
      PHCpack_Operations.Define_Output_File(name);
      put(PHCpack_Operations.output_file,natural32(target'last),target.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,natural32(start'last),start.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SOLUTIONS :");
      put(PHCpack_Operations.output_file,
          DoblDobl_Complex_Solutions.Length_Of(sols),
          natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
      new_line(PHCpack_Operations.output_file);
      Homotopy_Continuation_Parameters_io.put
        (PHCpack_Operations.output_file,homconpars.all);
      if mhom > 1 then
        nvr := natural32(start'last) - mhom;
        new_line(PHCpack_Operations.output_file);
        Series_Path_Trackers.Write_Partition
          (PHCpack_Operations.output_file,nvr,mhom,idz);
        new_line(PHCpack_Operations.output_file);
      end if;
      Drivers_to_Series_Trackers.DoblDobl_Track
        (PHCpack_Operations.output_file,
         target'last,sols,homconpars.all,mhom,idz,verbose);
    end if;
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(sols);
  end DoblDobl_Track;

  procedure QuadDobl_Track
              ( name : in string; localfile,verbose : in boolean;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Tracks the solution paths in standard precision,
  --   writing no output if name is empty,
  --   otherwise, creates an output file with the given name.
  --   The file name is the defined output file if localfile is false,
  --   otherwise, if localfile, the the file is a local variable.
  --   The value for mhom is 0, 1, 2 or higher, depending whether
  --   affine, 1-homogeneous, or m-homogeneous coordinates are applied.

    use QuadDobl_Complex_Numbers_cv;

    file : file_type;
    start,target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;
    nvr : natural32;

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

    gamma : constant Standard_Complex_Numbers.Complex_Number
          := homconpars.gamma;
    qd_gamma : constant QuadDobl_Complex_Numbers.Complex_Number
             := Standard_to_QuadDobl_Complex(gamma);

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    PHCpack_Operations.Retrieve_Target_System(target);
    if mhom = 0 then
      QuadDobl_Homotopy.Create(target.all,start.all,tpow,qd_gamma);
    else
      QuadDobl_Homotopy.Create(target.all,start.all,1,qd_gamma);
      QuadDobl_Coefficient_Homotopy.Create(start.all,target.all,1,qd_gamma);
    end if;
    if name = "" then
      if not verbose then
        Drivers_to_Series_Trackers.QuadDobl_Track
          (target'last,sols,homconpars.all,mhom,idz); -- ,verbose);
      else
        Homotopy_Continuation_Parameters_io.put(homconpars.all);
        Drivers_to_Series_Trackers.QuadDobl_Track
          (standard_output,target'last,sols,homconpars.all,mhom,idz,verbose);
      end if;
    elsif localfile then
      Create(file,out_file,name);
      put(file,natural32(target'last),target.all);
      new_line(file);
      put_line(file,"THE START SYSTEM :");
      put(file,natural32(start'last),start.all);
      new_line(file);
      put_line(file,"THE START SOLUTIONS :");
      put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
          natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
      new_line(file);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.QuadDobl_Track
        (file,target'last,sols,homconpars.all,mhom,idz,verbose);
      close(file);
    else
      PHCpack_Operations.Define_Output_File(name);
      put(PHCpack_Operations.output_file,natural32(target'last),target.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,natural32(start'last),start.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SOLUTIONS :");
      put(PHCpack_Operations.output_file,
          QuadDobl_Complex_Solutions.Length_Of(sols),
          natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
      new_line(PHCpack_Operations.output_file);
      Homotopy_Continuation_Parameters_io.put
        (PHCpack_Operations.output_file,homconpars.all);
      if mhom > 1 then
        nvr := natural32(start'last) - mhom;
        new_line(PHCpack_Operations.output_file);
        Series_Path_Trackers.Write_Partition
          (PHCpack_Operations.output_file,nvr,mhom,idz);
        new_line(PHCpack_Operations.output_file);
      end if;
      Drivers_to_Series_Trackers.QuadDobl_Track
        (PHCpack_Operations.output_file,
         target'last,sols,homconpars.all,mhom,idz,verbose);
    end if;
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(sols);
  end QuadDobl_Track;

  function Pade_Continuation_Track_Paths
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(6));
    prc : constant natural32 := natural32(v_a(v_a'first));
    nbc : constant natural32 := natural32(v_a(v_a'first+1));
    vrb : constant natural32 := natural32(v_a(v_a'first+2));
    mhm : constant natural32 := natural32(v_a(v_a'first+3));
    lcf : constant natural32 := natural32(v_a(v_a'first+4));
    nvr : constant natural32 := natural32(v_a(v_a'first+5));
    verbose : constant boolean := (vrb > 0);
    localfile : constant boolean := (lcf > 0);
    idz : Standard_Natural_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Track_Paths ...");
    end if;
    if mhm > 1 then
      declare
        v_c : constant C_Double_Array(0..Interfaces.C.size_T(nvr-1))
            := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(nvr));
      begin
        idz := new Standard_Natural_Vectors.Vector'(1..integer32(nvr) => 0);
        for i in idz'range loop
          idz(i) := natural32(v_c(Interfaces.C.size_T(i-1)));
        end loop;
      end;
    end if;
    if nbc = 0 then
      case prc is
        when 0 => Standard_Track("",localfile,verbose,mhm,idz);
        when 1 => DoblDobl_Track("",localfile,verbose,mhm,idz);
        when 2 => QuadDobl_Track("",localfile,verbose,mhm,idz);
        when others => null;
      end case;
    else
      declare
        v_b : constant C_Integer_Array(0..Interfaces.C.size_T(nbc-1))
            := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc));
        name : constant string := C_Integer_Array_to_String(nbc,v_b);
      begin
        case prc is
          when 0 => Standard_Track(name,localfile,verbose,mhm,idz);
          when 1 => DoblDobl_Track(name,localfile,verbose,mhm,idz);
          when 2 => QuadDobl_Track(name,localfile,verbose,mhm,idz);
          when others => null;
        end case;
      end;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Track_Paths");
      end if;
      return 739;
  end Pade_Continuation_Track_Paths;

  procedure Standard_Initialize_Tracker ( homogeneous : in boolean ) is

  -- DESCRIPTION :
  --   Retrieves target and start system and initializes the
  --   Series-Pade tracker in standard double precision.
  --   If homogeneous, then homogeneous coordinates will be used,
  --   otherwise the tracking happens in the original affine coordinates.

    start,target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    Standard_SeriesPade_Tracker.Init(target,start,homogeneous);
  end Standard_Initialize_Tracker;

  procedure DoblDobl_Initialize_Tracker ( homogeneous : in boolean ) is

  -- DESCRIPTION :
  --   Retrieves target and start system and initializes the
  --   Series-Pade tracker in double double precision.
  --   If homogeneous, then homogeneous coordinates will be used,
  --   otherwise the tracking happens in the original affine coordinates.

    start,target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    DoblDobl_SeriesPade_Tracker.Init(target,start,homogeneous);
  end DoblDobl_Initialize_Tracker;

  procedure QuadDobl_Initialize_Tracker ( homogeneous : in boolean ) is

  -- DESCRIPTION :
  --   Retrieves target and start system and initializes the
  --   Series-Pade tracker in quad double precision.
  --   If homogeneous, then homogeneous coordinates will be used,
  --   otherwise the tracking happens in the original affine coordinates.

    start,target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    QuadDobl_SeriesPade_Tracker.Init(target,start,homogeneous);
  end QuadDobl_Initialize_Tracker;

  procedure Standard_Initialize_Tracker ( idx : in integer32 ) is

  -- DESCRIPTION :
  --   Retrieves the target system and initializes the
  --   Series-Pade tracker in standard double precision,
  --   with the continuation parameter defined in idx.

    target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(target);
    Standard_SeriesPade_Tracker.Init(target,idx);
  end Standard_Initialize_Tracker;

  procedure DoblDobl_Initialize_Tracker ( idx : in integer32 ) is

  -- DESCRIPTION :
  --   Retrieves the target system and initializes the
  --   Series-Pade tracker in double double precision,
  --   with the continuation parameter defined in idx.

    target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(target);
    DoblDobl_SeriesPade_Tracker.Init(target,idx);
  end DoblDobl_Initialize_Tracker;

  procedure QuadDobl_Initialize_Tracker ( idx : in integer32 ) is

  -- DESCRIPTION :
  --   Retrieves the target system and initializes the
  --   Series-Pade tracker in quad double precision,
  --   with the continuation parameter defined in idx.

    target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(target);
    QuadDobl_SeriesPade_Tracker.Init(target,idx);
  end QuadDobl_Initialize_Tracker;

  function Pade_Continuation_Artificial_Homotopy
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    hmg : constant natural32 := natural32(v_b(v_b'first+1));
    verbose : constant boolean := (vrb = 1);
    homogeneous : constant boolean := (hmg = 1);

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Artificial_Homotopy ...");
    end if;
    case prc is
      when 0 =>
        if verbose then
          put("Initializing homotopy in Series-Pade tracker ");
          put("in double precision ...");
        end if;
        Standard_SeriesPade_Tracker.Init(homconpars.all);
        Standard_Initialize_Tracker(homogeneous);
      when 1 =>
        if verbose then
          put("Initializing homotopy in Series-Pade tracker ");
          put_line("in double double precision ...");
        end if;
        DoblDobl_SeriesPade_Tracker.Init(homconpars.all);
        DoblDobl_Initialize_Tracker(homogeneous);
      when 2 =>
        if verbose then
          put("Initializing homotopy in Series-Pade tracker ");
          put_line("in quad double precision ...");
        end if;
        QuadDobl_SeriesPade_Tracker.Init(homconpars.all);
        QuadDobl_Initialize_Tracker(homogeneous);
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    if verbose
     then put_line(" done!");
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Artificial_Homotopy");
      end if;
      return 860;
  end Pade_Continuation_Artificial_Homotopy;

  function Pade_Continuation_Natural_Homotopy
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    idx : constant integer32 := integer32(v_b(v_b'first));
    vrb : constant natural32 := natural32(v_b(v_b'first+1));
    verbose : constant boolean := (vrb = 1);

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Natural_Homotopy ...");
    end if;
    case prc is
      when 0 =>
        if verbose then
          put("Initializing parameter homotopy in Series-Pade tracker ");
          put("in double precision ...");
        end if;
        Standard_SeriesPade_Tracker.Init(homconpars.all);
        Standard_Initialize_Tracker(idx);
      when 1 =>
        if verbose then
          put("Initializing parameter homotopy in Series-Pade tracker ");
          put("in double double precision ...");
        end if;
        DoblDobl_SeriesPade_Tracker.Init(homconpars.all);
        DoblDobl_Initialize_Tracker(idx);
      when 2 =>
        if verbose then
          put("Initializing parameter homotopy in Series-Pade tracker ");
          put("in quad double precision ...");
        end if;
        QuadDobl_SeriesPade_Tracker.Init(homconpars.all);
        QuadDobl_Initialize_Tracker(idx);
      when others => null;
    end case;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Natural_Homotopy");
      end if;
      return 878;
  end Pade_Continuation_Natural_Homotopy;

  procedure Standard_Initialize_Solution
              ( idx : in natural32; verbose : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Retrieves the solution at index idx from the solutions container
  --   in standard double precision and writes output if verbose.

    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    if verbose
     then put("initializing solution "); put(idx,1); put_line(" :");
    end if;
    Standard_Solutions_Container.Retrieve(idx,ls,fail);
    Standard_SeriesPade_Tracker.Init(ls);
    if verbose
     then put(ls.all); new_line;
    end if;
  end Standard_Initialize_Solution;

  procedure DoblDobl_Initialize_Solution
              ( idx : in natural32; verbose : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Retrieves the solution at index idx from the solutions container
  --   in double double precision and writes output if verbose.

    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    if verbose
     then put("initializing solution "); put(idx,1); put_line(" :");
    end if;
    DoblDobl_Solutions_Container.Retrieve(idx,ls,fail);
    DoblDobl_SeriesPade_Tracker.Init(ls);
    if verbose
     then put(ls.all); new_line;
    end if;
  end DoblDobl_Initialize_Solution;

  procedure QuadDobl_Initialize_Solution
              ( idx : in natural32; verbose : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Retrieves the solution at index idx from the solutions container
  --   in quad double precision and writes output if verbose.
  --   The failure code of the retrieval is returned in fail.

    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    if verbose
     then put("initializing solution "); put(idx,1); put_line(" :");
    end if;
    QuadDobl_Solutions_Container.Retrieve(idx,ls,fail);
    QuadDobl_SeriesPade_Tracker.Init(ls);
    if verbose
     then put(ls.all); new_line;
    end if;
  end QuadDobl_Initialize_Solution;

  function Pade_Continuation_Set_Start_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    prc : constant natural32 := natural32(v_a(v_a'first));
    idx : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Set_Start_Solution ...");
    end if;
    case prc is
      when 0 => Standard_Initialize_Solution(idx,verbose,fail);
      when 1 => DoblDobl_Initialize_Solution(idx,verbose,fail);
      when 2 => QuadDobl_Initialize_Solution(idx,verbose,fail);
      when others => put_line("Wrong value for the precision.");
    end case;
    if fail then
      if verbose
       then put_line("The initialization of the start solution failed.");
      end if;
      assign(1,a);
    else
      if verbose
       then put_line("The initialization of the start solution succeeded.");
      end if;
      assign(0,a);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Set_Start_Solution");
      end if;
      return 861;
  end Pade_Continuation_Set_Start_Solution;

  function Pade_Continuation_Next_Step
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Next_Step ...");
    end if;
    case prc is
      when 0 => Standard_SeriesPade_Tracker.Predict_and_Correct(fail,verbose);
      when 1 => DoblDobl_SeriesPade_Tracker.Predict_and_Correct(fail,verbose);
      when 2 => QuadDobl_SeriesPade_Tracker.Predict_and_Correct(fail,verbose);
      when others => put_line("Wrong value for the precision.");
    end case;
    if fail then
      if verbose
       then put_line("The predict-correct step failed.");
      end if;
      assign(1,a);
    else
      if verbose
       then put_line("The predict-correct step succeeded.");
      end if;
      assign(0,a);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Next_Step");
      end if;
      return 862;
  end Pade_Continuation_Next_Step;

  function Pade_Continuation_Set_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    prc : constant natural32 := natural32(v_a(v_a'first));
    idx : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_ls : Standard_Complex_Solutions.Link_to_Solution;
    dd_ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    qd_ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Set_Solution ...");
    end if;
    if verbose then
      put("Retrieving the current solution, at index "); put(idx,1);
      put_line(" ...");
    end if;
    case prc is
      when 0 =>
        st_ls := Standard_SeriesPade_Tracker.Get_Current_Solution;
        Standard_Solutions_Container.Replace(idx,st_ls,fail);
      when 1 =>
        dd_ls := DoblDobl_SeriesPade_Tracker.Get_Current_Solution;
        DoblDobl_Solutions_Container.Replace(idx,dd_ls,fail);
      when 2 =>
        qd_ls := QuadDobl_SeriesPade_Tracker.Get_Current_Solution;
        QuadDobl_Solutions_Container.Replace(idx,qd_ls,fail);
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    if fail then
      if verbose
       then put_line("Placement of the current solution failed.");
      end if;
      assign(1,a);
    else
      if verbose
       then put_line("Placement of the current solution succeeded.");
      end if;
      assign(0,a);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Set_Solution");
      end if;
      return 863;
  end Pade_Continuation_Set_Solution;

  function Pade_Continuation_Set_Predicted_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    prc : constant natural32 := natural32(v_a(v_a'first));
    idx : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_ls : Standard_Complex_Solutions.Link_to_Solution;
    dd_ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    qd_ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Set_Predicted_Solution ...");
    end if;
    if verbose then
      put("Retrieving the predicted solution, at index "); put(idx,1);
      put_line(" ...");
    end if;
    case prc is
      when 0 =>
        st_ls := Standard_SeriesPade_Tracker.Get_Predicted_Solution;
        Standard_Solutions_Container.Replace(idx,st_ls,fail);
      when 1 =>
        dd_ls := DoblDobl_SeriesPade_Tracker.Get_Predicted_Solution;
        DoblDobl_Solutions_Container.Replace(idx,dd_ls,fail);
      when 2 =>
        qd_ls := QuadDobl_SeriesPade_Tracker.Get_Predicted_Solution;
        QuadDobl_Solutions_Container.Replace(idx,qd_ls,fail);
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    if fail then
      if verbose
       then put_line("Placement of the predicted solution failed.");
      end if;
      assign(1,a);
    else
      if verbose
       then put_line("Placement of the predicted solution succeeded.");
      end if;
      assign(0,a);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Set_Predicted_Solution");
      end if;
      return 919;
  end Pade_Continuation_Set_Predicted_Solution;

  function Pade_Continuation_Pole_Radius
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    frp : double_float;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Pole_Radius ...");
    end if;
    case prc is
      when 0 =>
        frp := Standard_SeriesPade_Tracker.Get_Current_Pole_Radius;
      when 1 =>
        frp := hi_part(DoblDobl_SeriesPade_Tracker.Get_Current_Pole_Radius);
      when 2 =>
        frp := hihi_part(QuadDobl_SeriesPade_Tracker.Get_Current_Pole_Radius);
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    Assign(frp,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Pole_Radius");
      end if;
      return 865;
  end Pade_Continuation_Pole_Radius;

  function Pade_Continuation_Closest_Pole
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    st_cfp : Standard_Complex_Numbers.Complex_Number;
    dd_cfp : DoblDobl_Complex_Numbers.Complex_Number;
    qd_cfp : QuadDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;
    use QuadDobl_Complex_Numbers_cv;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Closest_Pole ...");
    end if;
    case prc is
      when 0 =>
        st_cfp := Standard_SeriesPade_Tracker.Get_Current_Closest_Pole;
      when 1 =>
        dd_cfp := DoblDobl_SeriesPade_Tracker.Get_Current_Closest_Pole;
        st_cfp := DoblDobl_Complex_to_Standard(dd_cfp);
      when 2 =>
        qd_cfp := QuadDobl_SeriesPade_Tracker.Get_Current_Closest_Pole;
        st_cfp := QuadDobl_Complex_to_Standard(qd_cfp);
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    Assign(st_cfp,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Closest_Pole");
      end if;
      return 866;
  end Pade_Continuation_Closest_Pole;

  function Standard_Pole
             ( leadidx,poleidx : integer32; verbose : boolean )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the pole at position cffidx in component leadidx
  --   of the current poles in standard double precision.
  --   If verbose, then extra output is written to screen.

    res : Standard_Complex_Numbers.Complex_Number;
    poles: constant Standard_Complex_VecVecs.Link_to_VecVec
        := Standard_SeriesPade_Tracker.Get_Current_Poles;

  begin
    res := poles(leadidx)(poleidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end Standard_Pole;

  function DoblDobl_Pole
             ( leadidx,poleidx : integer32; verbose : boolean )
             return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the pole at position cffidx in component leadidx
  --   of the current poles in double double precision.
  --   If verbose, then extra output is written to screen.

    res : DoblDobl_Complex_Numbers.Complex_Number;
    poles: constant DoblDobl_Complex_VecVecs.Link_to_VecVec
        := DoblDobl_SeriesPade_Tracker.Get_Current_Poles;

  begin
    res := poles(leadidx)(poleidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end DoblDobl_Pole;

  function QuadDobl_Pole
             ( leadidx,poleidx : integer32; verbose : boolean )
             return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the pole at position cffidx in component leadidx
  --   of the current poles in quad double precision.
  --   If verbose, then extra output is written to screen.

    res : QuadDobl_Complex_Numbers.Complex_Number;
    poles: constant QuadDobl_Complex_VecVecs.Link_to_VecVec
        := QuadDobl_SeriesPade_Tracker.Get_Current_Poles;

  begin
    res := poles(leadidx)(poleidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end QuadDobl_Pole;

  function Pade_Continuation_Get_Pole
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    prc : constant natural32 := natural32(v_a(v_a'first));
    leadidx : constant integer32 := integer32(v_a(v_a'first+1));
    poleidx : constant integer32 := integer32(v_a(v_a'first+2));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_pole : Standard_Complex_Numbers.Complex_Number;
    dd_pole : DoblDobl_Complex_Numbers.Complex_Number;
    qd_pole : QuadDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;
    use QuadDobl_Complex_Numbers_cv;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Get_Pole ...");
    end if;
    if verbose then
      put("Precision : ");
      case prc is
        when 0 => put("double");
        when 1 => put("double double");
        when 2 => put("quad double");
        when others => put("invalid!");
      end case;
      put("  lead idx : "); put(leadidx,1);
      put("  pole idx : "); put(poleidx,1); new_line;
    end if;
    case prc is
      when 0 => st_pole := Standard_Pole(leadidx,poleidx,verbose);
      when 1 => dd_pole := DoblDobl_Pole(leadidx,poleidx,verbose);
                st_pole := DoblDobl_Complex_to_Standard(dd_pole);
      when 2 => qd_pole := QuadDobl_Pole(leadidx,poleidx,verbose);
                st_pole := QuadDobl_Complex_to_Standard(qd_pole);
      when others => put_line("Invalid value for the precision.");
    end case;
    Assign(st_pole,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Get_Pole");
      end if;
      return 871;
  end Pade_Continuation_Get_Pole;

  function Pade_Continuation_T_Value
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    tval : double_float;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_T_Value ...");
    end if;
    case prc is
      when 0 => tval := Standard_SeriesPade_Tracker.Get_Current_t_Value;
      when 1 => tval := DoblDobl_SeriesPade_Tracker.Get_Current_t_Value;
      when 2 => tval := QuadDobl_SeriesPade_Tracker.Get_Current_t_Value;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(tval,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_T_Value");
      end if;
      return 867;
  end Pade_Continuation_T_Value;

  function Pade_Continuation_Step_Size
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    step : double_float;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Step_Size ...");
    end if;
    case prc is
      when 0 => step := Standard_SeriesPade_Tracker.Get_Current_Step_Size;
      when 1 => step := DoblDobl_SeriesPade_Tracker.Get_Current_Step_Size;
      when 2 => step := QuadDobl_SeriesPade_Tracker.Get_Current_Step_Size;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(step,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Step_Size");
      end if;
      return 868;
  end Pade_Continuation_Step_Size;

  function Pade_Continuation_Series_Step
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    step : double_float;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Series_Step ...");
    end if;
    case prc is
      when 0 => step := Standard_SeriesPade_Tracker.Get_Current_Series_Step;
      when 1 => step := DoblDobl_SeriesPade_Tracker.Get_Current_Series_Step;
      when 2 => step := QuadDobl_SeriesPade_Tracker.Get_Current_Series_Step;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(step,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Series_Step");
      end if;
      return 885;
  end Pade_Continuation_Series_Step;

  function Pade_Continuation_Pole_Step
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    step : double_float;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Pole_Step ...");
    end if;
    case prc is
      when 0 => step := Standard_SeriesPade_Tracker.Get_Current_Pole_Step;
      when 1 => step := DoblDobl_SeriesPade_Tracker.Get_Current_Pole_Step;
      when 2 => step := QuadDobl_SeriesPade_Tracker.Get_Current_Pole_Step;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(step,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Pole_Step");
      end if;
      return 886;
  end Pade_Continuation_Pole_Step;

  function Pade_Continuation_Estimated_Distance
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    dist : double_float;
    dd_dist : double_double;
    qd_dist : quad_double;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Estimated_Distance ...");
    end if;
    case prc is
      when 0 =>
        dist := Standard_SeriesPade_Tracker.Get_Current_Estimated_Distance;
      when 1 =>
        dd_dist := DoblDobl_SeriesPade_Tracker.Get_Current_Estimated_Distance;
        dist := hi_part(dd_dist);
      when 2 => 
        qd_dist := QuadDobl_SeriesPade_Tracker.Get_Current_Estimated_Distance;
        dist := hihi_part(qd_dist);
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(dist,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Estimated_Distance");
      end if;
      return 887;
  end Pade_Continuation_Estimated_Distance;

  function Pade_Continuation_Hessian_Step
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    step : double_float;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Hessian_Step ...");
    end if;
    case prc is
      when 0 => step := Standard_SeriesPade_Tracker.Get_Current_Hessian_Step;
      when 1 => step := DoblDobl_SeriesPade_Tracker.Get_Current_Hessian_Step;
      when 2 => step := QuadDobl_SeriesPade_Tracker.Get_Current_Hessian_Step;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(step,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Hessian_Step");
      end if;
      return 888;
  end Pade_Continuation_Hessian_Step;

  function Standard_Series_Coefficient
             ( leadidx,cffidx : integer32; verbose : boolean )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current series vectors in standard double precision.
  --   If verbose, then extra output is written to screen.

    res : Standard_Complex_Numbers.Complex_Number;
    srv : constant Standard_Complex_Series_Vectors.Link_to_Vector
        := Standard_SeriesPade_Tracker.Get_Current_Series_Vector;

  begin
    res := srv(leadidx).cff(cffidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end Standard_Series_Coefficient;

  function DoblDobl_Series_Coefficient
             ( leadidx,cffidx : integer32; verbose : boolean )
             return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current series vectors in double double precision.
  --   If verbose, then extra output is written to screen.

    res : DoblDobl_Complex_Numbers.Complex_Number;
    srv : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector
        := DoblDobl_SeriesPade_Tracker.Get_Current_Series_Vector;

  begin
    res := srv(leadidx).cff(cffidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end DoblDobl_Series_Coefficient;

  function QuadDobl_Series_Coefficient
             ( leadidx,cffidx : integer32; verbose : boolean )
             return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current series vectors in quad double precision.
  --   If verbose, then extra output is written to screen.

    res : QuadDobl_Complex_Numbers.Complex_Number;
    srv : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector
        := QuadDobl_SeriesPade_Tracker.Get_Current_Series_Vector;

  begin
    res := srv(leadidx).cff(cffidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end QuadDobl_Series_Coefficient;

  function Pade_Continuation_Series_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    prc : constant natural32 := natural32(v_a(v_a'first));
    leadidx : constant integer32 := integer32(v_a(v_a'first+1));
    cffidx : constant integer32 := integer32(v_a(v_a'first+2));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_cff : Standard_Complex_Numbers.Complex_Number;
    dd_cff : DoblDobl_Complex_Numbers.Complex_Number;
    qd_cff : QuadDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;
    use QuadDobl_Complex_Numbers_cv;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Series_Coefficient ...");
    end if;
    if verbose then
      put("Precision : ");
      case prc is
        when 0 => put("double");
        when 1 => put("double double");
        when 2 => put("quad double");
        when others => put("invalid!");
      end case;
      put("  lead idx : "); put(leadidx,1);
      put("  cff idx : "); put(cffidx,1); new_line;
    end if;
    case prc is
      when 0 => st_cff := Standard_Series_Coefficient(leadidx,cffidx,verbose);
      when 1 => dd_cff := DoblDobl_Series_Coefficient(leadidx,cffidx,verbose);
                st_cff := DoblDobl_Complex_to_Standard(dd_cff);
      when 2 => qd_cff := QuadDobl_Series_Coefficient(leadidx,cffidx,verbose);
                st_cff := QuadDobl_Complex_to_Standard(qd_cff);
      when others => put_line("Invalid precision.");
    end case;
    Assign(st_cff,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Series_Coefficient");
      end if;
      return 869;
  end Pade_Continuation_Series_Coefficient;

  function Standard_Pade_Coefficient
             ( leadidx,cffidx : integer32; num : natural32;
               verbose : boolean )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current Pade vectors in standard double precision.
  --   If num = 0, then the denominator coefficient will be returned,
  --   otherwise, the numerator coefficient will be returned.
  --   If verbose, then extra output is written to screen.

    res : Standard_Complex_Numbers.Complex_Number;
    pv : constant Standard_Pade_Approximants.Link_to_Pade_Vector
       := Standard_SeriesPade_Tracker.Get_Current_Pade_Vector;

    use Standard_Pade_Approximants;

  begin
    if num = 0 then
      declare
        c : constant Standard_Complex_Vectors.Vector
          := Denominator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    else
      declare
        c : constant Standard_Complex_Vectors.Vector
          := Numerator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    end if;
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end Standard_Pade_Coefficient;

  function DoblDobl_Pade_Coefficient
             ( leadidx,cffidx : integer32; num : natural32;
               verbose : boolean )
             return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current Pade vectors in double double precision.
  --   If num = 0, then the denominator coefficient will be returned,
  --   otherwise, the numerator coefficient will be returned.
  --   If verbose, then extra output is written to screen.

    res : DoblDobl_Complex_Numbers.Complex_Number;
    pv : constant DoblDobl_Pade_Approximants.Link_to_Pade_Vector
       := DoblDobl_SeriesPade_Tracker.Get_Current_Pade_Vector;

    use DoblDobl_Pade_Approximants;

  begin
    if num = 0 then
      declare
        c : constant DoblDobl_Complex_Vectors.Vector
          := Denominator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    else
      declare
        c : constant DoblDobl_Complex_Vectors.Vector
          := Numerator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    end if;
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end DoblDobl_Pade_Coefficient;

  function QuadDobl_Pade_Coefficient
             ( leadidx,cffidx : integer32; num : natural32;
               verbose : boolean )
             return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current Pade vectors in quad double precision.
  --   If num = 0, then the denominator coefficient will be returned,
  --   otherwise, the numerator coefficient will be returned.
  --   If verbose, then extra output is written to screen.

    res : QuadDobl_Complex_Numbers.Complex_Number;
    pv : constant QuadDobl_Pade_Approximants.Link_to_Pade_Vector
       := QuadDobl_SeriesPade_Tracker.Get_Current_Pade_Vector;

    use QuadDobl_Pade_Approximants;

  begin
    if num = 0 then
      declare
        c : constant QuadDobl_Complex_Vectors.Vector
          := Denominator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    else
      declare
        c : constant QuadDobl_Complex_Vectors.Vector
          := Numerator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    end if;
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end QuadDobl_Pade_Coefficient;

  function Pade_Continuation_Pade_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    prc : constant natural32 := natural32(v_a(v_a'first));
    num : constant natural32 := natural32(v_a(v_a'first+1));
    leadidx : constant integer32 := integer32(v_a(v_a'first+2));
    cffidx : constant integer32 := integer32(v_a(v_a'first+3));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_cff : Standard_Complex_Numbers.Complex_Number;
    dd_cff : DoblDobl_Complex_Numbers.Complex_Number;
    qd_cff : QuadDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;
    use QuadDobl_Complex_Numbers_cv;

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Pade_Coefficient ...");
    end if;
    if verbose then
      put("Precision : ");
      case prc is
        when 0 => put("double");
        when 1 => put("double double");
        when 2 => put("quad double");
        when others => put("invalid!");
      end case;
      put("  lead idx : "); put(leadidx,1);
      put("  cff idx : "); put(cffidx,1);
      if num = 0
       then put_line("  denominator");
       else put_line("  numernator");
      end if;
    end if;
    case prc is
      when 0 =>
        st_cff := Standard_Pade_Coefficient(leadidx,cffidx,num,verbose);
      when 1 =>
        dd_cff := DoblDobl_Pade_Coefficient(leadidx,cffidx,num,verbose);
        st_cff := DoblDobl_Complex_to_Standard(dd_cff);
      when 2 =>
        qd_cff := QuadDobl_Pade_Coefficient(leadidx,cffidx,num,verbose);
        st_cff := QuadDobl_Complex_to_Standard(qd_cff);
      when others =>
        put_line("Invalid value for the precision.");
    end case;
    Assign(st_cff,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Pade_Coefficient");
      end if;
      return 870;
  end Pade_Continuation_Pade_Coefficient;

  function Pade_Continuation_Clear_Data
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));

  begin
    if vrblvl > 0 then
      put("-> in pade_continuation_interface.");
      put_line("Pade_Continuation_Clear_Data ...");
    end if;
    case prc is
      when 0 => Standard_SeriesPade_Tracker.Clear;
      when 1 => DoblDobl_SeriesPade_Tracker.Clear;
      when 2 => QuadDobl_SeriesPade_Tracker.Clear;
      when others => put_line("Wrong value for the precision.");
    end case;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in pade_continuation_interface.");
        put_line("Pade_Continuation_Clear_Data");
      end if;
      return 864;
  end Pade_Continuation_Clear_Data;

end Pade_Continuation_Interface;
