with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Root_Refining_Parameters;          use Root_Refining_Parameters;
with Multitasking_Root_Refiners;        use Multitasking_Root_Refiners;

procedure ts_mtverify is

-- DESCRIPTION :
--   Verification of a list of solutions with multitasking.

  procedure Main is

    lq : Link_to_Laur_Sys;
    sols : Solution_List;
    nt : integer32 := 0;
    file : file_type;
    ans : character;
    epsxa,epsfa,tolsing : double_float;
    numbit,maxit : natural32 := 0;
    deflate,wout,reporting : boolean;

  begin
    Standard_System_and_Solutions_io.get(lq,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    put(file,lq'last,1); new_line(file); put(file,lq.all);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    Standard_Menu_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line;
    put("Do you want output to screen during refinement ? (y/n) ");
    Ask_Yes_or_No(ans);
    reporting := (ans = 'y');
    new_line;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      put_line("The system is a genuine Laurent polynomial system.");
      if reporting then
        Reporting_Multitasking_Root_Refiner
          (file,nt,lq.all,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
      else
        Silent_Multitasking_Root_Refiner
          (file,nt,lq.all,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
      end if;
    else
      put_line("The system is a positive Laurent polynomial system.");
      declare
        use Standard_Laur_Poly_Convertors;
        p : Poly_Sys(lq'range) := Positive_Laurent_Polynomial_System(lq.all);
      begin
        if reporting then
          Reporting_Multitasking_Root_Refiner
            (file,nt,p,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
        else
          Silent_Multitasking_Root_Refiner
            (file,nt,p,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
        end if;
        Clear(p);
      end;
    end if;
  end Main;

begin
  Main;
end ts_mtverify;
