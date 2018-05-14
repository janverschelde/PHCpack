with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
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

  procedure Driver_for_DEMiCs_Algorithm
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

    dim : constant integer32 := p'last;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    mv : natural32;
    mcc : Mixed_Subdivision;
    timer : Timing_Widget;

  begin
    Extract_Supports(p,mix,sup,false); -- verbose is false
    tstart(timer);
    Call_DEMiCs(mix,sup,false);
    tstop(timer);
    declare
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    begin
      Process_Output(dim,mix,sup,lifsup,mcc,false);
      new_line(file);
      put_line(file,"The lifted supports :");
      Floating_Mixed_Subdivisions_io.put(file,lifsup);
    end;
    put_line(file,"The mixed-cell configuration :");
    Floating_Mixed_Subdivisions_io.put(file,natural32(dim),mix.all,mcc,mv);
    put(file,"The mixed volume : "); put(file,mv,1); new_line(file);
    new_line(file);
    print_times(file,timer,"DEMiCs Algorithm");
  end Driver_for_DEMiCs_Algorithm;

end Drivers_for_DEMiCs_Algorithm;
