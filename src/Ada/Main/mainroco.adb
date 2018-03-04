with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Drivers_for_Root_Counts;            use Drivers_for_Root_Counts;
with Write_Seed_Number;
with Greeting_Banners;
with Bye_Bye_Message;

procedure mainroco ( infilename,outfilename : in string ) is
 
  procedure Read_System
              ( filename : in string;
                lp : out Link_to_Poly_Sys ) is

  -- DESCRIPTION :
  --   If the string on input is not empty,
  --   then the file will be opened for reading.
  --   An exception will be handled if something goes wrong
  --   with the reading of the system on file.
  --   On return is a pointer to a polynomial system,
  --   or to null if something went wrong.

    file : file_type;

  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      Close(file);
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      lp := null; return;
  end Read_System;

  function Is_Square ( p : Poly_Sys ) return boolean is

  -- DESRIPTOIN :
  --   Returns true if the system p has 
  --   as many variables as equations.
  --   Writes an error message and returns false
  --   if the system is not square.

    use Standard_Complex_Polynomials;

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    if nv = nq then
      return true;
    else
      new_line;
      put("The number of equations : "); put(nq,1); new_line;
      put("The number of variables : "); put(nv,1); new_line;
      put_line("The system is not square!");
      new_line;
      return false;
    end if;
  end Is_Square;

  procedure Main is

  -- DESCRIPTION :
  --   Handles the input and output before calling the driver.

    lp : Link_to_Poly_Sys := null;
    outft : file_type;

  begin
    Read_System(infilename,lp);
    if lp = null
     then new_line; get(lp);
    end if;
    if Is_Square(lp.all) then
      Create_Output_File(outft,outfilename);
      put(outft,lp'last,1);
      new_line(outft);
      put(outft,lp.all);
      declare
        q : Poly_Sys(lp'range);
        qsols : Solution_List;
        rc : natural32;
      begin
        Driver_for_Root_Counts(outft,lp.all,q,false,qsols,rc);
      end; 
      new_line(outft);
      put_line(outft,Bye_Bye_Message);
      Write_Seed_Number(outft);
      put_line(outft,Greeting_Banners.Version);
      Close(outft);
    end if;
  end Main;

begin
  Main;
end mainroco;
