with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;     use Standard_Laur_Poly_Convertors;
with Standard_Complex_Laur_Randomizers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Supports_of_Polynomial_Systems;
with Floating_Lifting_Functions;        use Floating_Lifting_Functions;
with Floating_Lifting_Utilities;        use Floating_Lifting_Utilities;
with Floating_Mixed_Subdivisions_io;    use Floating_Mixed_Subdivisions_io;
with Standard_Stable_Homotopies;        use Standard_Stable_Homotopies;
with Stable_Polyhedral_Continuation;    use Stable_Polyhedral_Continuation;

package body Test_Stable_Mixed_Volumes is

  procedure Check_Minimal_Degrees
              ( p : in Laur_Sys; is_Laur : out boolean ) is
  begin
    is_Laur := false;
    for i in p'range loop
      declare
        d : constant Degrees := Minimal_Degrees(p(i));
      begin
        put("Minimal degrees of p("); put(i,1); put(") : ");
        Standard_Integer_Vectors_io.put(d.all);
        if Negative(d) then
          put_line("  Laurent polynomial");
          is_Laur := true;
        else    
          put_line("  common polynomial");
        end if;
      end;
    end loop;
    if is_Laur
     then put_line("The system is a Laurent polynomial system.");
     else put_line("The system is a common polynomial system.");
    end if;
    if Is_Genuine_Laurent(p)
     then put_line("The system is a Laurent polynomial system.");
     else put_line("The system is a common polynomial system.");
    end if;
  end Check_Minimal_Degrees;

  procedure Check_Supports ( p : in Laur_Sys; is_done : out boolean ) is

    s : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);
    origin : constant Standard_Integer_Vectors.Vector(p'range)
           := (p'range => 0);

  begin
    is_done := true;
    for i in s'range loop
      if Lists_of_Integer_Vectors.Is_In(s(i),origin) then
        put("The origin belongs to support ");
      else
        put("The origin does not belong to support ");
        is_done := false;
      end if;
      put(i,1); put_line(".");
    end loop;
  end Check_Supports;

  procedure Check_Zero_Types 
             ( n,r : in integer32; p : in Laur_Sys;
               mix : in Standard_Integer_Vectors.Link_to_Vector;
               mcc : in Mixed_Subdivision ) is

    lif : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r)
        := Lifted_Supports(r,mcc);
    b : constant double_float := Lifting_Bound(lif);
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    ztp : Standard_Integer_Vectors.Vector(1..n);
    cnt : natural32 := 0;
    nbz : integer32;

  begin
    put("The lifting bound : "); put(b); new_line;
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      ztp := Zero_Type(mic.nor.all,b,mic.pts.all);
      nbz := Number_of_Zeroes(ztp);
      cnt := cnt + 1;
      put("Type of Cell "); put(cnt,1); put(" : ");
      Standard_Integer_Vectors_io.put(ztp);
      if nbz < 0 then
        put_line(" spurious cell");
      elsif nbz = 0 then
        put_line(" original cell");
      else
        put(" #zeroes : "); put(nbz,1); new_line;
      end if;
      put_line("-> its inner normal :");
      Standard_Floating_Vectors_io.put_line(mic.nor.all);
      if nbz > 0 then
        put_line("-> inner normal after eliminating zeroes : ");
        put_line(Eliminate_Zeroes(mic.nor.all,ztp,nbz));
        put_line("The cell after elimination :");
        put(natural32(n),mix.all,Substitute_Zeroes(mic,nbz,ztp));
        put_line("The system after elimination :");
        put_line(Substitute_Zeroes(p,ztp,nbz));
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Check_Zero_Types;

  procedure Check_Zero_Types ( p : in Laur_Sys ) is

    file : file_type;
    mcc : Mixed_Subdivision;
    n,r : integer32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;

  begin
    new_line;
    put_line("Reading a file name for mixed cells...");
    Read_Name_and_Open_File(file);
    get(file,natural32(n),natural32(r),mix,mcc);
    put("Read "); put(Length_Of(mcc),1); put_line(" cells from file.");
    Check_Zero_Types(n,r,p,mix,mcc);
  end Check_Zero_Types;
             
  procedure Call_Polyhedral_Continuation ( p : in Laur_Sys ) is

    q : constant Laur_Sys(p'range)
      := Standard_Complex_Laur_Randomizers.Complex_Randomize1(p);
    qsols : Solution_List;
    qfile,infile,outfile : file_type;
    mcc,orgmcc,stbmcc : Mixed_Subdivision;
    n,r : integer32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;

  begin
   -- Random_Coefficient_Systems.Add_Constants(q);
    new_line;
    put_line("Reading a file name to write the random coefficient system...");
    Read_Name_and_Create_File(qfile);
    put_line(qfile,q);
    new_line;
    put_line("Reading a file name for mixed cells...");
    Read_Name_and_Open_File(infile);
    get(infile,natural32(n),natural32(r),mix,mcc);
    put("Read "); put(Length_Of(mcc),1); put_line(" cells from file.");
    declare
      lif : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r)
          := Lifted_Supports(r,mcc);
      b : constant double_float := Lifting_Bound(lif);
      orgcnt,stbcnt : natural32;
      ans : character;
    begin
      Split_Original_Cells(mcc,b,orgmcc,stbmcc,orgcnt,stbcnt);
      put("Number of cells without artificial origin : ");
      put(orgcnt,1); new_line;
      put("#extra stable cells with artificial origin : ");
      put(stbcnt,1); new_line;
      if stbcnt > 0 then
        put("Do you want intermediate output to file ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'n' then
          Silent_Polyhedral_Continuation(q,b,mix,lif,stbmcc,qsols);
        else
          new_line;
          put_line("Reading a file name for intermediate output...");
          Read_Name_and_Create_File(outfile);
          Reporting_Polyhedral_Continuation(outfile,q,b,mix,lif,stbmcc,qsols);
        end if;
      end if;
    end;
    if not Is_Null(qsols) then
      new_line(qfile);
      put_line(qfile,"THE SOLUTIONS :");
      put(qfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    end if;
    close(qfile);
  end Call_Polyhedral_Continuation;

  procedure Main is

    lp : Link_to_Laur_Sys;
    is_Laur,done : boolean;
    b : double_float;
    ans : character;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("MENU for testing stable mixed volumes ...");
    put_line("  1. test minimal degrees of the system;");
    put_line("  2. read a stable mixed subdivision and classify cells;");
    put_line("  3. prepare for polyhedral homotopy continuation.");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    if ans = '1' then
      put_line("Your system : "); put(lp.all);
      Check_Minimal_Degrees(lp.all,is_Laur);
      Check_Supports(lp.all,done);
      if done then
        put_line("The origin belongs to all supports, done.");
      else
        put_line("The origin does not belong to all supports.");
        b := Lifting_Bound(lp.all);
        put("The lifting bound is "); put(b); put_line(".");
      end if;
    elsif ans = '2' then
      Check_Zero_Types(lp.all);
    else
      Call_Polyhedral_Continuation(lp.all);
    end if;
  end Main;

end Test_Stable_Mixed_Volumes;
