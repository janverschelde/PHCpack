with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Sampling_Machine;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Standard_Stacked_Sample_Grids;
with Multprec_Stacked_Sample_Grids;
with Make_Sample_Grids;                 use Make_Sample_Grids;
with Standard_Trace_Interpolators;
with Multprec_Trace_Interpolators;

procedure ts_elim is

-- DESCRIPTION :
--   The numerical eliminator proceeds in three steps :
--     1) read an embedded system with generic points;
--     2) determine the elimination order;
--     3) choose a projection operation;
--     4) launch the sampler/projector/interpolator.

  procedure Standard_Eliminate
                ( file : in file_type;
                  p : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  dim : in natural32 ) is

  -- DESCRIPTION :
  --   Given the embedded system with its solutions, this procedure
  --   initializes the sampler and calls the interpolation routines.
  --   All calculations are done with standard floating-point arithmetic.

    use Standard_Stacked_Sample_Grids,Standard_Trace_Interpolators;

    sli : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    deg : constant integer32 := integer32(Length_Of(sps));
    ip : Standard_Complex_Polynomials.Poly;
 
  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if dim = 1 then
      declare
        grid : Array_of_Standard_Sample_Lists(0..deg);
        eps,dst : double_float;
        t : Trace_Interpolator1;  
      begin
        Standard_Rectangular_Grid_Creator
          (file,sps,natural32(deg),grid,eps,dst);
        put(file,"Maximal error of the samples on the grid : ");
        put(file,eps,3); new_line(file);
        put(file,"Minimal distance between samples in one list in grid :");
        put(file,dst,3); new_line(file);
        t := Create(file,grid);
        put(file,"Maximal residual of evaluation at the grid : ");
        put(file,Maximal_Error(t,grid),3); new_line(file);    
        ip := Expand(t);
      end;
    else
      declare
        grid : Stacked_Sample_Grid(integer32(dim),deg);
        t : Trace_Interpolator;
      begin
        Standard_Stacked_Grid_Creator(file,sps,true,grid); 
        t := Create(file,grid,deg);
        put(file,"Maximal residual of evaluation at the grid : ");
        put(file,Maximal_Error(t,grid),3); new_line(file);    
        ip := Expand(t);
      end;
    end if;
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,ip);
  end Standard_Eliminate;

  procedure Multprec_Eliminate
                ( file : in file_type;
                  ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                  mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  dim,size : in natural32 ) is       

  -- DESCRIPTION :
  --   Given the embedded system with its solutions, this procedure
  --   initializes the sampler and calls the interpolation routines.
  --   Sample refinement and interpolation is done with multi-precision
  --   arithmetic.

    use Multprec_Stacked_Sample_Grids,Multprec_Trace_Interpolators;

    sli : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(ep,dim);
    sps : Standard_Sample_List := Create(sols,sli);
    deg : constant integer32 := integer32(Length_Of(sps));
    ip : Multprec_Complex_Polynomials.Poly;

  begin
    Sampling_Machine.Initialize(ep,mp,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    new_line;
    put_line("See the output file for results...");
    new_line;                                                
    if dim = 1 then
      declare
        grid : Array_of_Multprec_Sample_Lists(0..deg);
        eps,dst : double_float;
        t : Trace_Interpolator1;
       begin
        Multprec_Rectangular_Grid_Creator
          (file,sps,natural32(deg),size,grid,eps,dst);
        put(file,"Maximal error of the samples on the grid : ");
        put(file,eps,3); new_line(file);
        put(file,"Minimal distance between samples in one list in grid :");
        put(file,dst,3); new_line(file);
        t := Create(file,grid);
        ip := Expand(t);
      end;
    else
      declare
        grid : Stacked_Sample_Grid(integer32(dim),deg);
        t : Trace_Interpolator;
      begin
        Multprec_Stacked_Grid_Creator(file,sps,true,size,grid);
        t := Create(file,grid,deg);
        ip := Expand(t);
      end;
    end if;
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,ip);
  end Multprec_Eliminate;

  procedure Main is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim,deci,size : natural32 := 0;

  begin
    new_line;
    put_line
      ("Numerical Elimination by Sampling, Projecting, and Interpolating.");
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    Determine_Order(lp.all,sols);
    new_line;
    put("Give the number of decimal places (<= 16 is standard) : ");
    get(deci);
    new_line;
    if deci > 16 then
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(lp.all,mp,size,dim);
      Multprec_Eliminate(file,lp.all,mp.all,sols,dim,size);
    else
      Standard_Eliminate(file,lp.all,sols,dim);
    end if;
  end Main;

begin
  Main;
end ts_elim;
