with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials;
with Standard_Homotopy;
with Series_and_Homotopies;
with Solution_Drops;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Test_Series_Predictors;

package body Standard_Homotopy_Convolutions_io is

  function Make_Homotopy ( nq,idx,deg : integer32 ) return Link_to_System is

    res : Link_to_System;
    hom : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
        := Standard_Homotopy.Homotopy_System;
    serhom : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
           := Series_and_Homotopies.Create(hom,idx);

  begin
    Complex_Series_and_Polynomials.Set_Degree(serhom,deg);
    res := Make_Convolution_System(serhom,natural32(deg));
    Standard_Cseries_Poly_Systems.Clear(serhom);
    return res;
  end Make_Homotopy;

  procedure get ( deg : in integer32; h : out Link_to_System;
                  s : out Solution_List; idxpar : out integer32 ) is

    nbeq : integer32;
    sols : Solution_List;

  begin
    Test_Series_Predictors.Standard_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      idxpar := nbeq + 1; -- assuming the system is square ...
      s := sols;
    else
      declare
        dropsols : constant Standard_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        s := dropsols;
      end;
    end if;
    h := Make_Homotopy(nbeq,idxpar,deg);
  end get;

end Standard_Homotopy_Convolutions_io;
