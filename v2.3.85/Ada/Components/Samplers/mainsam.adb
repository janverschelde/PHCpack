with text_io;                           use text_io;
with Characters_and_Numbers;
with String_Splitters;                  use String_Splitters;
with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Continuation_Parameters;
with Standard_Continuation_Data;        use Standard_Continuation_Data;
with Standard_Continuation_Data_io;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Moving_Planes;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Intrinsic_Witness_Sets_io;         use Intrinsic_Witness_Sets_io;
with Standard_Intrinsic_Continuation;   use Standard_Intrinsic_Continuation;

procedure mainsam ( witset,logfile : in string ) is

  procedure Read_Witness_Set
              ( w : in string; p : out Link_to_Poly_Sys;
                sols : out Solution_List; dim : out natural32 ) is

  -- DESCRIPTION :
  --   Reads witness set k from file with name in the string w.

  -- ON ENTRY :
  --   w        name of file, if empty, the user is asked for a name.

  -- ON RETURN :
  --   p        polynomial equations defining the witness set;
  --   sols     points in the witness set;

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

  procedure Local_Intrinsic_Continuation
              ( file : in file_type; name : in string;
                n,k : in integer32; output : in boolean;
                f : in Poly_Sys; start : in Matrix;
                esols : in Solution_List ) is

  -- DESCRIPTION :
  --   Deforms the given plane p into a random target plane,
  --   running intrinsic continuation on all solutions.

  -- ON ENTRY :
  --   file     file to write diagnostics on;
  --   name     name of the output file;
  --   n        ambient dimension, number of original variables;
  --   k        codimension of the solution set;
  --   output   if intermediate output wanted during tracking;
  --   f        original polynomial system;
  --   esols    extrinsic coordinates of the solutions.

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

  procedure Setup_Local_Coordinates
              ( file : in file_type; name : in string;
                ep : in Poly_Sys; sols : in Solution_List;
                d : in integer32 ) is

  -- DESCRIPTION :
  --   Prepares the setup for working with local intrinsic coordinates.

  -- ON ENTRY:
  --   file     file to write diagnostics on;
  --   name     name of the output file;
  --   ep       embedded polynomial system with d slack variables;
  --   sols     generic points on the algebraic set;
  --   d        dimension of the witness set.

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
    Local_Intrinsic_Continuation(file,name,n,k,output,pp,pla,sols);
  end Setup_Local_Coordinates;

  procedure Main is

    p : Link_to_Poly_Sys;
    w : Solution_List;
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
  end Main;

begin
  Main;
end mainsam;
