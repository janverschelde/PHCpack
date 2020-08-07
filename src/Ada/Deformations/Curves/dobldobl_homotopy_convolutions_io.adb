with text_io;                            use text_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials;
with DoblDobl_Homotopy;
with Series_and_Homotopies;
with Solution_Drops;
with Homotopy_Series_Readers;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Test_Series_Predictors;
with Series_Path_Trackers;               use Series_Path_Trackers;

package body DoblDobl_Homotopy_Convolutions_io is

  function Make_Homotopy ( nq,idx,deg : integer32 ) return Link_to_System is

    res : Link_to_System;
    hom : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
        := DoblDobl_Homotopy.Homotopy_System;
    serhom : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
           := Series_and_Homotopies.Create(hom,idx);

  begin
    Complex_Series_and_Polynomials.Set_Degree(serhom,deg);
    res := Make_Convolution_System(serhom,natural32(deg));
    DoblDobl_Cseries_Poly_Systems.Clear(serhom);
    return res;
  end Make_Homotopy;

  procedure get ( deg : in integer32; h : out Link_to_System;
                  s : out Solution_List; idxpar : out integer32 ) is

    nbeq : integer32;
    sols : Solution_List;

  begin
    Test_Series_Predictors.DoblDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      idxpar := nbeq + 1; -- assuming the system is square ...
      s := sols;
    else
      declare
        dropsols : constant DoblDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        s := dropsols;
      end;
    end if;
    h := Make_Homotopy(nbeq,idxpar,deg);
  end get;

  procedure get ( deg : in integer32; artificial : in boolean;
                  gamma : in Complex_Number; h : out Link_to_System;
                  s : out Solution_List; idxpar : out integer32;
                  mhom : out natural32; z : out Link_to_Partition;
                  idz : out Standard_Natural_Vectors.Link_to_Vector ) is

    nbq,nvr : integer32;
    sols : Solution_List;

  begin
    if artificial then
      DoblDobl_Define_Homotopy(nbq,nvr,gamma,mhom,z,idz,sols);
      idxpar := nbq + 1;
      s := sols;
    else
      mhom := 0;
      Homotopy_Series_Readers.DoblDobl_Parameter_Reader(nbq,nvr,idxpar,sols);
      declare
        dropsols : constant DoblDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        s := dropsols;
      end;
    end if;
    h := Make_Homotopy(nbq,idxpar,deg);
  end get;

end DoblDobl_Homotopy_Convolutions_io;
