with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Floating_Mixed_Subdivisions_io;
with Drivers_for_Static_Lifting;
with Drivers_for_MixedVol_Algorithm;
with DEMiCs_Output_Data;
with DEMiCs_Algorithm;                   use DEMiCs_Algorithm;

package body Drivers_for_DEMiCs_Algorithm is

  procedure DEMiCs_Algorithm_Info is

    i : array(1..10) of string(1..65);

  begin
    i(1) :="  The DEMiCs Algorithm calls the code of Tomohiko Mizutani,      ";
    i(2) :="Akiko Takeda, and Masakazu Kojima.  The algorithm is described in";
    i(3) :="Discrete Comput. Geom. 37(3):351-367, 2007.  The software DEMiCs ";
    i(4) :="is published in Software for Algebraic Geometry, Springer, 2008. ";
    i(5) :="DEMiCs stands for Dynamic Enumeration of Mixed Cells and applies ";
    i(6) :="a greedy strategy to run through the tree of face combinations   ";
    i(7) :="which span all mixed cells.  For many different Newton polytopes ";
    i(8) :="DEMiCs is faster than MixedVol, producing cells at a faster pace.";
    i(9) :="Compared to other lift-and-prune strategies, only random lifting ";
   i(10) :="values on the supports are supported.                            ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end DEMiCs_Algorithm_Info;

  procedure BlackBox_DEMiCs_Algorithm
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
              lif : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32 ) is  

    dim : constant integer32 := p'last;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    Extract_Supports(p,mix,sup,false); -- verbose is false
    Call_DEMiCs(mix,sup,false);
    declare
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    begin
      Process_Output(dim,mix,sup,lifsup,mcc,false);
      lif := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(lifsup);
    end;
    mv := natural32(DEMiCs_Output_Data.mixed_volume);
  end BlackBox_DEMiCs_Algorithm;

  procedure BlackBox_DEMiCs_Algorithm
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
              lif : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32 ) is  

    dim : constant integer32 := p'last;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    Extract_Supports(p,mix,sup,false); -- verbose is false
    Call_DEMiCs(mix,sup,false);
    declare
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    begin
      Process_Output(dim,mix,sup,lifsup,mcc,false);
      lif := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(lifsup);
    end;
    mv := natural32(DEMiCs_Output_Data.mixed_volume);
  end BlackBox_DEMiCs_Algorithm;

  procedure Process_DEMiCs_Output
              ( file : in file_type;
                mcc2file,ranstart : in boolean;
                subfile,ranfile : in out file_type;
                dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                timer : in Timing_Widget ) is

  -- DESCRIPTION :
  --   Processes the output of the DEMiCs algorithm.

  -- ON ENTRY :
  --   file     for output;
  --   mcc2file if the mixed subdivision has to be written on file;
  --   ranstart if a random coefficient system is needed;
  --   subfile  file opened for output if mcc2file;
  --   ranfile  file to write the random coefficient system on;
  --   dim      dimension of the points;
  --   mix      type of mixture of the supports;
  --   sup      points in the supports;
  --   timer    contains timings of the DEMiCs algorithm.

    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    mcc : Mixed_Subdivision;
    mv : natural32;
    q : Standard_Complex_Laur_Systems.Laur_Sys(1..dim);
    qsols : Solution_List;
    qtimer : Timing_Widget;

    use Drivers_for_MixedVol_Algorithm;

  begin
    Process_Output(dim,mix,sup,lifsup,mcc,false);
    new_line(file);
    put_line(file,"The lifted supports :");
    Floating_Mixed_Subdivisions_io.put(file,lifsup);
    put_line(file,"The mixed-cell configuration :");
    Floating_Mixed_Subdivisions_io.put(file,natural32(dim),mix.all,mcc,mv);
    if mcc2file then
      Floating_Mixed_Subdivisions_io.put(subfile,natural32(dim),mix.all,mcc);
      close(subfile);
    end if;
    put(file,"The mixed volume : "); put(file,mv,1); new_line(file);
    new_line(file);
    print_times(file,timer,"DEMiCs Algorithm");
    if ranstart then
      tstart(qtimer);
      Random_Coefficient_System(0,dim,mix.all,lifsup,mcc,q,qsols);
      tstop(qtimer);
      new_line(file);
      put_line(file,"RANDOM COEFFICIENT SYSTEM :");
      put_line(file,q);
      put_line(ranfile,q);
      new_line(file);
      put_line(file,"THE SOLUTIONS ");
      new_line(ranfile);
      put_line(ranfile,"THE SOLUTIONS ");
      put(file,Length_Of(qsols),natural32(dim),qsols);
      put(ranfile,Length_Of(qsols),natural32(dim),qsols);
      close(ranfile);
      new_line(file);
      print_times(file,qtimer,"polyhedral homotopies");
    end if;
  end Process_DEMiCs_Output;

  procedure Driver_for_DEMiCs_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

    dim : constant integer32 := p'last;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    timer : Timing_Widget;
    mcc2file : boolean := false;
    ranstart : boolean := false;
    subfile,ranfile : file_type;
    ans : character;

  begin
    new_line;
    Drivers_for_Static_Lifting.Prompt_for_File(mcc2file,subfile);
    new_line;
    put("Monitor the progress of the computations on screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    DEMiCs_Output_Data.monitor := (ans = 'y');
    new_line;
    put("Run polyhedral homotopies ");
    put("to solve a random coefficient system ? (y/n) ");
    Ask_Yes_or_No(ans);
    ranstart := (ans = 'y');
    if ranstart then
      put("Reading the name of the file ");
      put_line("for the random coefficient system ...");
      Read_Name_and_Create_File(ranfile);
    end if;
    new_line;
    put_line("See the output file for results ...");
    new_line;
    Extract_Supports(p,mix,sup,false); -- verbose is false
    tstart(timer);
    Call_DEMiCs(mix,sup,false);
    tstop(timer);
    Process_DEMiCs_Output
      (file,mcc2file,ranstart,subfile,ranfile,dim,mix,sup,timer);
  end Driver_for_DEMiCs_Algorithm;

end Drivers_for_DEMiCs_Algorithm;
