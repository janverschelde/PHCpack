with Characters_and_Numbers;
with String_Splitters;                  use String_Splitters;
with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Continuation_Parameters;
with Standard_Continuation_Data;
with Standard_Continuation_Data_io;
with DoblDobl_Continuation_Data;
with DoblDobl_Continuation_Data_io;
with QuadDobl_Continuation_Data;
with QuadDobl_Continuation_Data_io;
with Main_Poly_Continuation;            use Main_Poly_Continuation;
with Standard_Plane_Representations;
with Standard_Moving_Planes;
with DoblDobl_Plane_Representations;
with DoblDobl_Moving_Planes;
with QuadDobl_Plane_Representations;
with QuadDobl_Moving_Planes;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Intrinsic_Witness_Sets_io;         use Intrinsic_Witness_Sets_io;
with Standard_Intrinsic_Continuation;
with DoblDobl_Intrinsic_Continuation;
with QuadDobl_Intrinsic_Continuation;

package body Main_Samplers is

  procedure Read_Witness_Set
              ( w : in string;
                p : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    if w = "" then
      new_line;
      put_line("Reading a file name for a witness set ...");
      Read_Name_and_Open_File(file);
    else
      Open_Input_File(file,w);
    end if;
    Standard_Read_Embedding(file,p,sols,dim);
  end Read_Witness_Set;

  procedure Read_Witness_Set
              ( w : in string;
                p : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    if w = "" then
      new_line;
      put_line("Reading a file name for a witness set ...");
      Read_Name_and_Open_File(file);
    else
      Open_Input_File(file,w);
    end if;
    DoblDobl_Read_Embedding(file,p,sols,dim);
  end Read_Witness_Set;

  procedure Read_Witness_Set
              ( w : in string;
                p : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                dim : out natural32 ) is

    file : file_type;

  begin
    if w = "" then
      new_line;
      put_line("Reading a file name for a witness set ...");
      Read_Name_and_Open_File(file);
    else
      Open_Input_File(file,w);
    end if;
    QuadDobl_Read_Embedding(file,p,sols,dim);
  end Read_Witness_Set;

  procedure Local_Intrinsic_Continuation
              ( file : in file_type; name : in string;
                n,k : in integer32; output : in boolean;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                start : in Standard_Complex_Matrices.Matrix;
                esols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Continuation_Data;
    use Standard_Intrinsic_Continuation;

    target : constant Matrix(1..n,0..k)
           := Standard_Moving_Planes.One_Random_Direction(start);
          -- := Standard_Moving_Planes.Random_Plane(n,k);
    q : constant Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : constant Eval_Poly_Sys(1..k) := Create(q);
    jm : constant Jaco_Mat(1..k,1..n) := Create(q);
    jf : constant Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    timer : Timing_Widget;
    deg : constant natural32 := Length_Of(esols);
    degstr : constant string := Characters_and_Numbers.nConvert(deg);
    sa : Solu_Info_Array(1..integer32(deg)) := Shallow_Create(esols);

  begin
    if not output then
      tstart(timer);
      Silent_Local_LU_Continue(eq,jf,start,target,true,sa,pp,cp);
      tstop(timer);
    else
      new_line(file);
      tstart(timer);
      Reporting_Local_LU_Continue(file,eq,jf,start,target,true,sa,pp,cp);
      tstop(timer);
    end if;
    Write_Witness_Set_to_File(name,natural32(n),natural32(k),f,target,sa);
    new_line(file);
    print_times(file,timer,"computing " & degstr & " new witness points");
    new_line(file);
    Standard_Continuation_Data_io.Path_Tracking_Statistics(file,sa);
  end Local_Intrinsic_Continuation;

  procedure Local_Intrinsic_Continuation
              ( file : in file_type; name : in string;
                n,k : in integer32; output : in boolean;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                start : in DoblDobl_Complex_Matrices.Matrix;
                esols : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Continuation_Data;
    use DoblDobl_Intrinsic_Continuation;

    target : constant Matrix(1..n,0..k)
           := DoblDobl_Moving_Planes.One_Random_Direction(start);
    q : constant Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : constant Eval_Poly_Sys(1..k) := Create(q);
    jm : constant Jaco_Mat(1..k,1..n) := Create(q);
    jf : constant Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    timer : Timing_Widget;
    deg : constant natural32 := Length_Of(esols);
    degstr : constant string := Characters_and_Numbers.nConvert(deg);
    sa : Solu_Info_Array(1..integer32(deg)) := Shallow_Create(esols);

  begin
    if not output then
      tstart(timer);
      Silent_Local_LU_Continue(eq,jf,start,target,true,sa,pp,cp);
      tstop(timer);
    else
      new_line(file);
      tstart(timer);
      Reporting_Local_LU_Continue(file,eq,jf,start,target,true,sa,pp,cp);
      tstop(timer);
    end if;
    Write_Witness_Set_to_File(name,natural32(n),natural32(k),f,target,sa);
    new_line(file);
    print_times(file,timer,"computing " & degstr & " new witness points");
    new_line(file);
    DoblDobl_Continuation_Data_io.Path_Tracking_Statistics(file,sa);
  end Local_Intrinsic_Continuation;

  procedure Local_Intrinsic_Continuation
              ( file : in file_type; name : in string;
                n,k : in integer32; output : in boolean;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                start : in QuadDobl_Complex_Matrices.Matrix;
                esols : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Continuation_Data;
    use QuadDobl_Intrinsic_Continuation;

    target : constant Matrix(1..n,0..k)
           := QuadDobl_Moving_Planes.One_Random_Direction(start);
    q : constant Poly_Sys(1..k) := Witness_Sets.Make_Square(f,natural32(k));
    eq : constant Eval_Poly_Sys(1..k) := Create(q);
    jm : constant Jaco_Mat(1..k,1..n) := Create(q);
    jf : constant Eval_Jaco_Mat(1..k,1..n) := Create(jm);
    pp : constant Continuation_Parameters.Pred_Pars
       := Continuation_Parameters.Create_for_Path;
    cp : constant Continuation_Parameters.Corr_Pars
       := Continuation_Parameters.Create_for_Path;
    timer : Timing_Widget;
    deg : constant natural32 := Length_Of(esols);
    degstr : constant string := Characters_and_Numbers.nConvert(deg);
    sa : Solu_Info_Array(1..integer32(deg)) := Shallow_Create(esols);

  begin
    if not output then
      tstart(timer);
      Silent_Local_LU_Continue(eq,jf,start,target,true,sa,pp,cp);
      tstop(timer);
    else
      new_line(file);
      tstart(timer);
      Reporting_Local_LU_Continue(file,eq,jf,start,target,true,sa,pp,cp);
      tstop(timer);
    end if;
    Write_Witness_Set_to_File(name,natural32(n),natural32(k),f,target,sa);
    new_line(file);
    print_times(file,timer,"computing " & degstr & " new witness points");
    new_line(file);
    QuadDobl_Continuation_Data_io.Path_Tracking_Statistics(file,sa);
  end Local_Intrinsic_Continuation;

  procedure byebye is
  begin
    new_line;
    put_line("See the output file(s) for results ...");
    new_line;
  end byebye;

  procedure Setup_Local_Coordinates
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                d : in integer32 ) is

    use Standard_Complex_Matrices;
    use Standard_Complex_VecVecs;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;
    use Standard_Plane_Representations;

    n : constant integer32 := ep'last - d; -- remove slack variables
    k : constant integer32 := n-d;
    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    nbzero : constant integer32 := integer32(Number_of_Zero_Equations(p));
    pp : constant Poly_Sys := p(p'first..p'last-nbzero);
    s : constant VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    output : boolean;
    oc : natural32 := 0;

  begin
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    put("Co-dimension of the component : "); put(k,1); new_line;
    put("Degree of the solution set : "); put(Length_Of(sols),1); new_line;
    new_line;
    Continuation_Parameters.Tune(2);
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    output := (oc > 0);
    byebye;
    Local_Intrinsic_Continuation(file,name,n,k,output,pp,pla,sols);
  end Setup_Local_Coordinates;

  procedure Setup_Local_Coordinates
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                d : in integer32 ) is

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Plane_Representations;

    n : constant integer32 := ep'last - d; -- remove slack variables
    k : constant integer32 := n-d;
    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    nbzero : constant integer32 := integer32(Number_of_Zero_Equations(p));
    pp : constant Poly_Sys := p(p'first..p'last-nbzero);
    s : constant VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    output : boolean;
    oc : natural32 := 0;

  begin
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    put("Co-dimension of the component : "); put(k,1); new_line;
    put("Degree of the solution set : "); put(Length_Of(sols),1); new_line;
    new_line;
    Continuation_Parameters.Tune(2);
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    output := (oc > 0);
    byebye;
    Local_Intrinsic_Continuation(file,name,n,k,output,pp,pla,sols);
  end Setup_Local_Coordinates;

  procedure Setup_Local_Coordinates
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                d : in integer32 ) is

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Plane_Representations;

    n : constant integer32 := ep'last - d; -- remove slack variables
    k : constant integer32 := n-d;
    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    nbzero : constant integer32 := integer32(Number_of_Zero_Equations(p));
    pp : constant Poly_Sys := p(p'first..p'last-nbzero);
    s : constant VecVec(1..d) := Slices(ep,natural32(d));
    eqs : constant Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    pla : constant Matrix(1..n,0..k) := Orthogonalize(gen);
    output : boolean;
    oc : natural32 := 0;

  begin
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    put("Co-dimension of the component : "); put(k,1); new_line;
    put("Degree of the solution set : "); put(Length_Of(sols),1); new_line;
    new_line;
    Continuation_Parameters.Tune(2);
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    output := (oc > 0);
    byebye;
    Local_Intrinsic_Continuation(file,name,n,k,output,pp,pla,sols);
  end Setup_Local_Coordinates;

  procedure Sample_in_Standard_Precision ( witset,logfile : in string ) is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    w : Standard_Complex_Solutions.Solution_List;
    d : integer32;
    file : file_type;
    name : Link_to_String;

  begin
    Read_Witness_Set(witset,p,w,natural32(d));
    if logfile = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file,name);
    else
      Create_Output_File(file,logfile,name);
    end if;
    Setup_Local_Coordinates(file,name.all,p.all,w,d);
  end Sample_in_Standard_Precision;

  procedure Sample_in_DoblDobl_Precision ( witset,logfile : in string ) is

    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    w : DoblDobl_Complex_Solutions.Solution_List;
    d : integer32;
    file : file_type;
    name : Link_to_String;

  begin
    Read_Witness_Set(witset,p,w,natural32(d));
    if logfile = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file,name);
    else
      Create_Output_File(file,logfile,name);
    end if;
    Setup_Local_Coordinates(file,name.all,p.all,w,d);
  end Sample_in_DoblDobl_Precision;

  procedure Sample_in_QuadDobl_Precision ( witset,logfile : in string ) is

    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    w : QuadDobl_Complex_Solutions.Solution_List;
    d : integer32;
    file : file_type;
    name : Link_to_String;

  begin
    Read_Witness_Set(witset,p,w,natural32(d));
    if logfile = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file,name);
    else
      Create_Output_File(file,logfile,name);
    end if;
    Setup_Local_Coordinates(file,name.all,p.all,w,d);
  end Sample_in_QuadDobl_Precision;

  procedure Main ( witset,logfile : in string ) is

    ans : character;

  begin
    new_line;
    put_line("MENU to set the precision for sampling :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the level of precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Sample_in_Standard_Precision(witset,logfile);
      when '1' => Sample_in_DoblDobl_Precision(witset,logfile);
      when '2' => Sample_in_QuadDobl_Precision(witset,logfile);
      when others => null;
    end case;
  end Main;

end Main_Samplers;
