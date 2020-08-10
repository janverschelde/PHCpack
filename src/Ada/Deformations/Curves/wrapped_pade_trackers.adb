with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;
with Standard_Homotopy;
with Wrapped_Solution_Vectors;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Standard_Homotopy_Convolutions_io;
with Standard_Speelpenning_Convolutions;
with Standard_Coefficient_Circuits;
with Standard_Coefficient_Convolutions;
with Standard_Convolution_Splitters;
with Standard_Circuit_Makers;
with Predictor_Corrector_Trackers;

with Standard_Natural_Numbers;
 use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;
 use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems_io;
 use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
 use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
 use Standard_Complex_Solutions_io;

package body Wrapped_Pade_Trackers is

  procedure Set_Parameters is

    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;

  begin -- set gamma to one because natural parameter homotopy
    pars.gamma := Standard_Complex_Numbers.Create(1.0);
    Homotopy_Continuation_Parameters_io.Tune(pars);
    Homotopy_Continuation_Parameters.Construct(pars);
  end Set_Parameters;

  procedure Set_Parameters ( file : in file_type ) is

    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;

  begin -- set gamma to one because natural parameter homotopy
    pars.gamma := Standard_Complex_Numbers.Create(1.0);
    Homotopy_Continuation_Parameters_io.Tune(pars);
    Homotopy_Continuation_Parameters.Construct(pars);
    Homotopy_Continuation_Parameters_io.put(file,pars);
  end Set_Parameters;

  procedure Clear is
  begin
    Homotopy_Continuation_Parameters.Destruct;
  end Clear;

-- TRACKING ONE PATH WITH NO OUTPUT :

  procedure Run ( idxpar : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 ) is

    sols : Standard_Complex_Solutions.Solution_List
         := Wrapped_Solution_Vectors.Create(xt);

  begin
    if vrblvl > 0 then
      put_line("-> in wrapped_pade_trackers.Call_Path_Trackers 1 ...");
    end if;
    Run(idxpar,h,sols,vrblvl-1);
    xt(xt'first..xt'last-1) := Standard_Complex_Solutions.Head_Of(sols).v;
    xt(xt'last) := Standard_Complex_Solutions.Head_Of(sols).t;
    sol := Standard_Complex_Solutions.Head_Of(sols);
  end Run;

-- TRACKING ONE PATH WITH OUTPUT TO FILE :

  procedure Run ( file : in file_type; idxpar : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  verbose : in boolean := false;
                  vrblvl : in integer32 := 0 ) is

    sols : Standard_Complex_Solutions.Solution_List
         := Wrapped_Solution_Vectors.Create(xt);

  begin
    if vrblvl > 0 then
      put_line("-> in wrapped_pade_trackers.Call_Path_Trackers 2 ...");
    end if;
    Run(file,idxpar,h,sols,verbose,vrblvl-1);
    xt(xt'first..xt'last-1) := Standard_Complex_Solutions.Head_Of(sols).v;
    xt(xt'last) := Standard_Complex_Solutions.Head_Of(sols).t;
    sol := Standard_Complex_Solutions.Head_Of(sols);
  end Run;

-- TRACKING MANY PATHS WITHOUT OUTPUT TO FILE :

  procedure Run ( idxpar : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Polynomials;

    nbequ : constant integer32 := h'last;
    nbvar : constant integer32
          := integer32(Number_Of_Unknowns(h(h'first))) - 1;
    pars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
         := Homotopy_Continuation_Parameters.Retrieve;
    deg : constant integer32
        := integer32(pars.numdeg) + integer32(pars.dendeg) + 2;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    cffhom : Standard_Coefficient_Convolutions.Link_to_System;
    cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
    idz : constant Standard_Natural_Vectors.Link_to_Vector := null;

    use Standard_Homotopy_Convolutions_io;

  begin
    if vrblvl > 0 then
      put_line("-> in wrapped_pade_trackers.Call_Path_Trackers 3 ...");
    end if;
    Standard_Homotopy.Create(h,idxpar);
    if nbequ > nbvar then -- overdetermined homotopy
      put_line("will skip for now ...");
    else
      cnvhom := Make_Homotopy(nbequ,idxpar,deg);
      cffhom := Standard_Convolution_Splitters.Split(cnvhom);
      cfs := Standard_Circuit_Makers.Make_Coefficient_System(cffhom);
      abh := Standard_Coefficient_Circuits.Copy(cfs);
      Standard_Coefficient_Circuits.AbsVal(abh);
      Predictor_Corrector_Trackers.Track_All_Paths
        (cffhom,cfs,abh,sols,pars.all,0,idz,vrblvl-1);
    end if;
    Standard_Speelpenning_Convolutions.Clear(cnvhom);
    Standard_Coefficient_Convolutions.Clear(cffhom);
    Standard_Coefficient_Circuits.Clear(cfs);
    Standard_Coefficient_Circuits.Clear(abh);
    Standard_Homotopy.Clear;
  end Run;

-- TRACKING MANY PATHS WITH OUTPUT TO FILE :

  procedure Run ( file : in file_type; idxpar : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  verbose : in boolean := false;
                  vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Polynomials;

    nbequ : constant integer32 := h'last;
    nbvar : constant integer32
          := integer32(Number_Of_Unknowns(h(h'first))) - 1;
    pars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
         := Homotopy_Continuation_Parameters.Retrieve;
    deg : constant integer32
        := integer32(pars.numdeg) + integer32(pars.dendeg) + 2;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    cffhom : Standard_Coefficient_Convolutions.Link_to_System;
    cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
    idz : constant Standard_Natural_Vectors.Link_to_Vector := null;

    use Standard_Homotopy_Convolutions_io;

  begin
    if vrblvl > 0 then
      put_line("-> in wrapped_pade_trackers.Call_Path_Trackers 4 ...");
    end if;
    Standard_Homotopy.Create(h,idxpar);
    put_line(file,"homotopy for wrapped Pade trackers :");
    put(file,"nbequ : "); put(file,nbequ,1);
    put(file,"  idxpar : "); put(file,idxpar,1); new_line(file);
    put(file,h);
    new_line(file);
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    if nbequ > nbvar then -- overdetermined homotopy
      put_line("will skip for now ...");
    else
      cnvhom := Make_Homotopy(nbequ,idxpar,deg);
      cffhom := Standard_Convolution_Splitters.Split(cnvhom);
      cfs := Standard_Circuit_Makers.Make_Coefficient_System(cffhom);
      abh := Standard_Coefficient_Circuits.Copy(cfs);
      Standard_Coefficient_Circuits.AbsVal(abh);
      Predictor_Corrector_Trackers.Track_All_Paths
        (file,cffhom,cfs,abh,sols,pars.all,0,idz,true,vrblvl-1);
       -- verbose,vrblvl-1);
    end if;
    Standard_Speelpenning_Convolutions.Clear(cnvhom);
    Standard_Coefficient_Convolutions.Clear(cffhom);
    Standard_Coefficient_Circuits.Clear(cfs);
    Standard_Coefficient_Circuits.Clear(abh);
    Standard_Homotopy.Clear;
  end Run;

end Wrapped_Pade_Trackers;
