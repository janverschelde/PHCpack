with Timing_Package;                    use Timing_Package;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with DoblDobl_Random_Numbers;           use DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;           use QuadDobl_Random_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Sample_Points;                     use Sample_Points;
with DoblDobl_Sample_Points;            use DoblDobl_Sample_Points;
with QuadDobl_Sample_Points;            use QuadDobl_Sample_Points;
with Rectangular_Sample_Grids;          use Rectangular_Sample_Grids;
with DoblDobl_Rectangular_Sample_Grids; use DoblDobl_Rectangular_Sample_Grids;
with QuadDobl_Rectangular_Sample_Grids; use QuadDobl_Rectangular_Sample_Grids;

package body Make_Sample_Grids is

  procedure Write_Errors ( file : in file_type;
                           sps : in Standard_Sample_List ) is

    tmp : Standard_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Sample_Point(Head_Of(tmp)).err,3);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Errors;

  procedure Write_Errors ( file : in file_type;
                           sps : in DoblDobl_Sample_List ) is

    tmp : DoblDobl_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Sample_Point(Head_Of(tmp)).err,3);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Errors;

  procedure Write_Errors ( file : in file_type;
                           sps : in QuadDobl_Sample_List ) is

    tmp : QuadDobl_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Sample_Point(Head_Of(tmp)).err,3);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Errors;

  procedure Write_Errors ( file : in file_type;
                           sps : in Multprec_Sample_List ) is

    tmp : Multprec_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Sample_Point(Head_Of(tmp)).err,3);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Errors;

  procedure Standard_Rectangular_Grid_Creator
               ( file : in file_type; sps : in Standard_Sample_List;
                 m : in natural32; grid : out Array_of_Standard_Sample_Lists;
                 eps,dst : out double_float ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    grid := Create1(sps,m);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    tstart(timer);
    eps := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,eps,3); new_line(file);
    dst := Minimal_Distance(grid);
    put(file,"Minimal distance between samples in grid : ");
    put(file,dst,3); new_line(file);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid.");
    new_line(file);
  end Standard_Rectangular_Grid_Creator;

  procedure DoblDobl_Rectangular_Grid_Creator
               ( file : in file_type; sps : in DoblDobl_Sample_List;
                 m : in natural32; grid : out Array_of_DoblDobl_Sample_Lists;
                 eps,dst : out double_double ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    grid := Create1(sps,m);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    tstart(timer);
    eps := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,eps,3); new_line(file);
    dst := Minimal_Distance(grid);
    put(file,"Minimal distance between samples in grid : ");
    put(file,dst,3); new_line(file);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid.");
    new_line(file);
  end DoblDobl_Rectangular_Grid_Creator;

  procedure QuadDobl_Rectangular_Grid_Creator
               ( file : in file_type; sps : in QuadDobl_Sample_List;
                 m : in natural32; grid : out Array_of_QuadDobl_Sample_Lists;
                 eps,dst : out quad_double ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    grid := Create1(sps,m);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    tstart(timer);
    eps := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,eps,3); new_line(file);
    dst := Minimal_Distance(grid);
    put(file,"Minimal distance between samples in grid : ");
    put(file,dst,3); new_line(file);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid.");
    new_line(file);
  end QuadDobl_Rectangular_Grid_Creator;

  procedure Standard_Triangular_Grid_Creator
               ( file : in file_type; sps : in Standard_Sample_List;
                 m : in natural32; grid : out Array_of_Standard_Sample_Lists;
                 eps,dst : out double_float ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    grid := Triangular_Create1(sps,m);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    tstart(timer);
    eps := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,eps,3); new_line(file);
    dst := Minimal_Distance(grid);
    put(file,"Minimal distance between samples in grid : ");
    put(file,dst,3); new_line(file);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid.");
    new_line(file);
  end Standard_Triangular_Grid_Creator;

  procedure DoblDobl_Triangular_Grid_Creator
               ( file : in file_type; sps : in DoblDobl_Sample_List;
                 m : in natural32; grid : out Array_of_DoblDobl_Sample_Lists;
                 eps,dst : out double_double ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    grid := Triangular_Create1(sps,m);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    tstart(timer);
    eps := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,eps,3); new_line(file);
    dst := Minimal_Distance(grid);
    put(file,"Minimal distance between samples in grid : ");
    put(file,dst,3); new_line(file);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid.");
    new_line(file);
  end DoblDobl_Triangular_Grid_Creator;

  procedure QuadDobl_Triangular_Grid_Creator
               ( file : in file_type; sps : in QuadDobl_Sample_List;
                 m : in natural32; grid : out Array_of_QuadDobl_Sample_Lists;
                 eps,dst : out quad_double ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    grid := Triangular_Create1(sps,m);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    tstart(timer);
    eps := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,eps,3); new_line(file);
    dst := Minimal_Distance(grid);
    put(file,"Minimal distance between samples in grid : ");
    put(file,dst,3); new_line(file);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing the quality of the grid.");
    new_line(file);
  end QuadDobl_Triangular_Grid_Creator;

  procedure Multprec_Rectangular_Grid_Creator
               ( file : in file_type; sps : in out Standard_Sample_List;
                 m,size : in natural32;
                 grid : out Array_of_Multprec_Sample_Lists;
                 eps,dst : out double_float ) is

    timer : Timing_Widget;
    maxerr,mindst : Floating_Number;

  begin
    tstart(timer);
    grid := Create1(sps,m,size);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,maxerr,3); new_line(file);
    eps := Round(maxerr);
    Clear(maxerr);
    mindst := Minimal_Distance(grid);
    put(file,"Minimal distance between samples in grid : ");
    put(file,mindst,3); new_line(file);
    dst := Round(mindst);
    Clear(mindst);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing quality of the grid");
    new_line(file);
  end Multprec_Rectangular_Grid_Creator;

  procedure Multprec_Triangular_Grid_Creator
               ( file : in file_type; sps : in out Standard_Sample_List;
                 m,size : in natural32;
                 grid : out Array_of_Multprec_Sample_Lists;
                 eps,dst : out double_float ) is

    timer : Timing_Widget;
    maxerr,mindst : Floating_Number;

  begin
    tstart(timer);
    grid := Triangular_Create1(sps,m,size);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    tstart(timer);
    maxerr := Maximal_Error(grid);
    put(file,"Maximal error of the samples in the grid : ");
    put(file,maxerr,3); new_line(file);
    eps := Round(maxerr);
    Clear(maxerr);
    mindst := Minimal_Distance(grid);
    put(file,"Minimal distance between samples in grid : ");
    put(file,mindst,3); new_line(file);
    dst := Round(mindst);
    Clear(mindst);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Testing quality of the grid");
    new_line(file);
  end Multprec_Triangular_Grid_Creator;

  procedure Standard_Stacked_Grid_Creator
               ( file : in file_type; sps : in Standard_Sample_List;
                 full : in boolean;
                 grid : out
                        Standard_Stacked_Sample_Grids.Stacked_Sample_Grid ) is

    use Standard_Stacked_Sample_Grids;
    timer : Timing_Widget;

  begin
    tstart(timer);
    if full
     then grid := Create_Full(file,sps);
     else grid := Create(file,sps);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    put_line(file,"The errors on samples in the grid : ");
    Write_Errors(file,grid);
    put(file,"The maximal error on samples : ");
    put(file,Maximal_Error(grid),3); new_line(file);
    if full
     then Write_Full_Grid_Values(file,grid);
     else Write_Grid_Values(file,grid);
    end if;
  end Standard_Stacked_Grid_Creator;

  procedure DoblDobl_Stacked_Grid_Creator
               ( file : in file_type; sps : in DoblDobl_Sample_List;
                 full : in boolean;
                 grid : out
                        DoblDobl_Stacked_Sample_Grids.Stacked_Sample_Grid ) is

    use DoblDobl_Stacked_Sample_Grids;
    timer : Timing_Widget;

  begin
    tstart(timer);
    if full
     then grid := Create_Full(file,sps);
     else grid := Create(file,sps);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    put_line(file,"The errors on samples in the grid : ");
    Write_Errors(file,grid);
    put(file,"The maximal error on samples : ");
    put(file,Maximal_Error(grid),3); new_line(file);
    if full
     then Write_Full_Grid_Values(file,grid);
     else Write_Grid_Values(file,grid);
    end if;
  end DoblDobl_Stacked_Grid_Creator;

  procedure QuadDobl_Stacked_Grid_Creator
               ( file : in file_type; sps : in QuadDobl_Sample_List;
                 full : in boolean;
                 grid : out
                        QuadDobl_Stacked_Sample_Grids.Stacked_Sample_Grid ) is

    use QuadDobl_Stacked_Sample_Grids;
    timer : Timing_Widget;

  begin
    tstart(timer);
    if full
     then grid := Create_Full(file,sps);
     else grid := Create(file,sps);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    put_line(file,"The errors on samples in the grid : ");
    Write_Errors(file,grid);
    put(file,"The maximal error on samples : ");
    put(file,Maximal_Error(grid),3); new_line(file);
    if full
     then Write_Full_Grid_Values(file,grid);
     else Write_Grid_Values(file,grid);
    end if;
  end QuadDobl_Stacked_Grid_Creator;

  procedure Multprec_Stacked_Grid_Creator
              ( file : in file_type; sps : in Standard_Sample_List;
                full : in boolean; size : in natural32;
                grid : out 
                       Multprec_Stacked_Sample_Grids.Stacked_Sample_Grid ) is

    use Multprec_Stacked_Sample_Grids;
    timer : Timing_Widget;

  begin
    tstart(timer);
    if full
     then grid := Create_Full(file,sps,size);
     else grid := Create(file,sps,size);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of the grid of samples");
    new_line(file);
    put_line(file,"The errors on samples in the grid : ");
    Write_Errors(file,grid);
    put(file,"The maximal error on samples : ");
    put(file,Maximal_Error(grid),3); new_line(file);
    if full
     then Write_Full_Grid_Values(file,grid);
     else Write_Grid_Values(file,grid);
    end if;
  end Multprec_Stacked_Grid_Creator;

  procedure Standard_Test_Samples
               ( file : in file_type; sps : in Standard_Sample_List;
                 sli : in Standard_Complex_VecVecs.VecVec;
                 testsps : out Standard_Sample_List ) is

    testsps_last : Standard_Sample_List;
    newsli : Standard_Complex_VecVecs.VecVec(sli'range);

  begin
    for i in sli'range loop
      newsli(i) := new Standard_Complex_Vectors.Vector'(sli(i).all);
      newsli(i)(0) := Random1;
    end loop;
    Sample(sps,newsli,testsps,testsps_last);
    put_line(file,"Errors on the test samples : ");
    Write_Errors(file,testsps); new_line(file);
  end Standard_Test_Samples;

  procedure DoblDobl_Test_Samples
               ( file : in file_type; sps : in DoblDobl_Sample_List;
                 sli : in DoblDobl_Complex_VecVecs.VecVec;
                 testsps : out DoblDobl_Sample_List ) is

    testsps_last : DoblDobl_Sample_List;
    newsli : DoblDobl_Complex_VecVecs.VecVec(sli'range);

  begin
    for i in sli'range loop
      newsli(i) := new DoblDobl_Complex_Vectors.Vector'(sli(i).all);
      newsli(i)(0) := Random1;
    end loop;
    Sample(sps,newsli,testsps,testsps_last);
    put_line(file,"Errors on the test samples : ");
    Write_Errors(file,testsps); new_line(file);
  end DoblDobl_Test_Samples;

  procedure QuadDobl_Test_Samples
               ( file : in file_type; sps : in QuadDobl_Sample_List;
                 sli : in QuadDobl_Complex_VecVecs.VecVec;
                 testsps : out QuadDobl_Sample_List ) is

    testsps_last : QuadDobl_Sample_List;
    newsli : QuadDobl_Complex_VecVecs.VecVec(sli'range);

  begin
    for i in sli'range loop
      newsli(i) := new QuadDobl_Complex_Vectors.Vector'(sli(i).all);
      newsli(i)(0) := Random1;
    end loop;
    Sample(sps,newsli,testsps,testsps_last);
    put_line(file,"Errors on the test samples : ");
    Write_Errors(file,testsps); new_line(file);
  end QuadDobl_Test_Samples;

  procedure Multprec_Test_Samples
               ( file : in file_type; sps : in Standard_Sample_List;
                 sli : in Standard_Complex_VecVecs.VecVec;
                 testsps : out Multprec_Sample_List ) is

    testsps_last : Multprec_Sample_List;
    newsli : Standard_Complex_VecVecs.VecVec(sli'range);

  begin
    for i in sli'range loop
      newsli(i) := new Standard_Complex_Vectors.Vector'(sli(i).all);
      newsli(i)(0) := Random1;
    end loop;
    Sample(sps,newsli,testsps,testsps_last);
    put_line(file,"Errors on the test samples : ");
    Write_Errors(file,testsps); new_line(file);
  end Multprec_Test_Samples;

end Make_Sample_Grids;
