with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
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

  procedure Run ( n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    sols : Solution_List := Wrapped_Solution_Vectors.Create(xt);
    nbequ : constant integer32 := h'last;
    pars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
         := Homotopy_Continuation_Parameters.Retrieve;
    deg : constant integer32
        := integer32(pars.numdeg) + integer32(pars.dendeg) + 2;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    cffhom : Standard_Coefficient_Convolutions.Link_to_System;
    cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
    idz : constant Standard_Natural_Vectors.Link_to_Vector := null;

  begin
    if vrblvl > 0 then
      put_line("-> in wrapped_pade_trackers.Call_Path_Trackers 1 ...");
    end if;
    Standard_Homotopy.Create(h,n+1);
    if nbequ > n then -- overdetermined homotopy
      put_line("will skip for now ...");
    else
      cnvhom := Standard_Homotopy_Convolutions_io.Make_Homotopy(nbequ,n,deg);
      cffhom := Standard_Convolution_Splitters.Split(cnvhom);
      cfs := Standard_Circuit_Makers.Make_Coefficient_System(cffhom);
      abh := Standard_Coefficient_Circuits.Copy(cfs);
      Standard_Coefficient_Circuits.AbsVal(abh);
      Predictor_Corrector_Trackers.Track_All_Paths
        (cffhom,cfs,abh,sols,pars.all,0,idz,vrblvl-1);
    end if;
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    Standard_Speelpenning_Convolutions.Clear(cnvhom);
    Standard_Coefficient_Convolutions.Clear(cffhom);
    Standard_Coefficient_Circuits.Clear(cfs);
    Standard_Coefficient_Circuits.Clear(abh);
    Standard_Homotopy.Clear;
  end Run;

-- TRACKING ONE PATH WITH OUTPUT TO FILE :

  procedure Run ( file : in file_type; n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  xt : in out Standard_Complex_Vectors.Vector;
                  sol : out Standard_Complex_Solutions.Link_to_Solution;
                  verbose : in boolean := false;
                  vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    sols : Solution_List := Wrapped_Solution_Vectors.Create(xt);
    nbequ : constant integer32 := h'last;
    pars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
         := Homotopy_Continuation_Parameters.Retrieve;
    deg : constant integer32
        := integer32(pars.numdeg) + integer32(pars.dendeg) + 2;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    cffhom : Standard_Coefficient_Convolutions.Link_to_System;
    cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
    idz : constant Standard_Natural_Vectors.Link_to_Vector := null;

  begin
    if vrblvl > 0 then
      put_line("-> in wrapped_pade_trackers.Call_Path_Trackers 2 ...");
    end if;
    Standard_Homotopy.Create(h,n+1);
    if nbequ > n then -- overdetermined homotopy
      put_line("will skip for now ...");
    else
      cnvhom := Standard_Homotopy_Convolutions_io.Make_Homotopy(nbequ,n,deg);
      cffhom := Standard_Convolution_Splitters.Split(cnvhom);
      cfs := Standard_Circuit_Makers.Make_Coefficient_System(cffhom);
      abh := Standard_Coefficient_Circuits.Copy(cfs);
      Standard_Coefficient_Circuits.AbsVal(abh);
      Predictor_Corrector_Trackers.Track_All_Paths
        (file,cffhom,cfs,abh,sols,pars.all,0,idz,verbose,vrblvl-1);
    end if;
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    Standard_Speelpenning_Convolutions.Clear(cnvhom);
    Standard_Coefficient_Convolutions.Clear(cffhom);
    Standard_Coefficient_Circuits.Clear(cfs);
    Standard_Coefficient_Circuits.Clear(abh);
    Standard_Homotopy.Clear;
  end Run;

-- TRACKING MANY PATHS WITHOUT OUTPUT TO FILE :

  procedure Run ( n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  vrblvl : in integer32 := 0 ) is

    nbequ : constant integer32 := h'last;
    pars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
         := Homotopy_Continuation_Parameters.Retrieve;
    deg : constant integer32
        := integer32(pars.numdeg) + integer32(pars.dendeg) + 2;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    cffhom : Standard_Coefficient_Convolutions.Link_to_System;
    cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
    idz : constant Standard_Natural_Vectors.Link_to_Vector := null;

  begin
    if vrblvl > 0 then
      put_line("-> in wrapped_pade_trackers.Call_Path_Trackers 3 ...");
    end if;
    Standard_Homotopy.Create(h,n+1);
    if nbequ > n then -- overdetermined homotopy
      put_line("will skip for now ...");
    else
      cnvhom := Standard_Homotopy_Convolutions_io.Make_Homotopy(nbequ,n,deg);
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

  procedure Run ( file : in file_type; n : in integer32;
                  h : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  verbose : in boolean := false;
                  vrblvl : in integer32 := 0 ) is

    nbequ : constant integer32 := h'last;
    pars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
         := Homotopy_Continuation_Parameters.Retrieve;
    deg : constant integer32
        := integer32(pars.numdeg) + integer32(pars.dendeg) + 2;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    cffhom : Standard_Coefficient_Convolutions.Link_to_System;
    cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
    idz : constant Standard_Natural_Vectors.Link_to_Vector := null;

  begin
    if vrblvl > 0 then
      put_line("-> in wrapped_pade_trackers.Call_Path_Trackers 4 ...");
    end if;
    Standard_Homotopy.Create(h,n+1);
    if nbequ > n then -- overdetermined homotopy
      put_line("will skip for now ...");
    else
      cnvhom := Standard_Homotopy_Convolutions_io.Make_Homotopy(nbequ,n,deg);
      cffhom := Standard_Convolution_Splitters.Split(cnvhom);
      cfs := Standard_Circuit_Makers.Make_Coefficient_System(cffhom);
      abh := Standard_Coefficient_Circuits.Copy(cfs);
      Standard_Coefficient_Circuits.AbsVal(abh);
      Predictor_Corrector_Trackers.Track_All_Paths
        (file,cffhom,cfs,abh,sols,pars.all,0,idz,verbose,vrblvl-1);
    end if;
    Standard_Speelpenning_Convolutions.Clear(cnvhom);
    Standard_Coefficient_Convolutions.Clear(cffhom);
    Standard_Coefficient_Circuits.Clear(cfs);
    Standard_Coefficient_Circuits.Clear(abh);
    Standard_Homotopy.Clear;
  end Run;

end Wrapped_Pade_Trackers;
