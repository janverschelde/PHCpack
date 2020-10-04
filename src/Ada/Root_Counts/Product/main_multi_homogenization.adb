with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Numbers_io;                         use Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Prod_Systems_io;   use Standard_Complex_Prod_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with Degree_Structure;
with Multi_Homogeneous_Start_Systems;    use Multi_Homogeneous_Start_Systems;

package body Main_Multi_Homogenization is

  procedure Multi_Homogenization_Info is

    i : array(1..17) of string(1..65);

  begin
    i( 1):="  A multi-homogeneous Bezout  number  is  based  on  a  tuple  of";
    i( 2):="partitions  of  the set of unknowns.  For every polynomial in the";
    i( 3):="system, a different partition can model its structure.           ";
    i( 4):="  The corresponding start system is a linear-product system:  the";
    i( 5):="i-th  equation  is  the  product  of linear equations with random";
    i( 6):="coefficients in the unknowns of the set of  the  partition.   The";
    i( 7):="number  of  factors  in  the product for the i-th equation of the";
    i( 8):="start system equals the  product  of  the  degrees  of  the  i-th";
    i( 9):="polynomial  in  the  original  system  w.r.t.  every  set  in the";
    i(10):="partition.                                                       ";
    i(11):="  Given a  tuple  of  partitions,  the  multi-homogeneous  Bezout";
    i(12):="number  equals  the  number  of  solutions  of  the corresponding";
    i(13):="linear-product start system.   Before  the  construction  of  the";
    i(14):="start system, a multi-homogeneous Bezout number is first computed";
    i(15):="in a formal way as a generalized permanent of a degree matrix.  A";
    i(16):="heuristic  procedure  is  available  for  generating  a  tuple of";
    i(17):="partitions.                                                      ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Multi_Homogenization_Info;

  procedure Menu_Prompt ( choice : out character; gb : in natural32 ) is

    ans : character;

  begin
    new_line;
    put_line("MENU for Multi-Homogeneous Bezout Numbers :");
    put     ("  0. exit - current Bezout number is "); put(gb,1); new_line;
    put_line("  1. Apply heuristic partitioner");
    put_line("  2. Evaluate your own tuple of partitions.");
    put("Type 0, 1, or 2 to make your choice : "); 
    Ask_Alternative(ans,"012"); choice := ans;
  end Menu_Prompt;

  procedure Menu_Handle
              ( file : in file_type; choice : in character;
                p : in Poly_Sys; gb : in out natural32 ) is

    m : natural32;

  begin
    case choice is
      when '1' => gb := Degree_Structure.Generalized_Bezout_Number(p);
      when '2' =>
        for i in p'range loop
          put("Give the number of sets in partition ");
          put(i,1); put(" : "); Read_Natural(m);
          put("Give "); put(m,1); put(" sets : ");
          declare
            zz : Partition(1..m);
          begin
            Create(zz,p'length); get(zz);
            Degree_Structure.Put(p,natural32(i),m,zz);
            Clear(zz);
          end;
        end loop;
        gb := Degree_Structure.Generalized_Bezout_Number;
      when others => null;
    end case;
    Write_Results(Standard_Output,p,gb);
    Write_Results(file,p,gb);
  end Menu_Handle;

  procedure Define_Partitions
              ( file : in file_type; p : in Poly_Sys;
                gb : in out natural32; b : out natural32 ) is

    timer : timing_widget;
    method : character;

  begin
    new_line(file);
    put_line(file,"MULTI-HOMOGENIZATION :");
    tstart(timer);
    loop
      Menu_Prompt(method,gb);
      exit when method = '0';
      Menu_Handle(file,method,p,gb);
    end loop;
    tstop(timer);
    b := gb;
    new_line(file);
    print_times(file,timer,"Computation of multi-homogeneous Bezout number");
  end Define_Partitions;

  procedure Construct_Start_System
              ( file : in file_type; p : in Poly_Sys;
                q : out Poly_Sys; rq : out Prod_Sys;
                qsols : out Solution_List ) is

    ans : character;
    timer : timing_widget;
    qq : Poly_Sys(p'range);
    qqsols : Solution_List;

  begin
    new_line;
    put("Do you want a multi-homogeneous start system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      tstart(timer);
      GBQ(p,qq,rq,qqsols);
      tstop(timer);
      Save_Results(rq,qqsols);
      new_line(file);
      put_line(file,"MULTI-HOMOGENEOUS START SYSTEM :");
      put_line(file,qq);
      q := qq; qsols := qqsols;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,Length_Of(qqsols),natural32(Head_Of(qqsols).n),qqsols);
      new_line(file);
      print_times(file,timer,
                  "Construction of multi-homogeneous start system");
    end if;
  end Construct_Start_System;

  procedure Write_Results
              ( file : in file_type; p : in Poly_Sys;
                gb : in natural32 ) is

    m : natural32;

  begin
    new_line(file);
    put(file,"  multi-homogeneous Bezout number is ");
    put(file,gb,1); new_line(file);
    put_line(file,"  with partitions :");
    for i in p'range loop
      m := Degree_Structure.Get(natural32(i));
      declare
        z : partition(1..m);
        dg : Standard_Natural_Vectors.Vector(1..integer32(m));
      begin
        Degree_Structure.Get(natural32(i),z,dg);
        put(file,"     partition for equation "); put(file,i,2);
        put(file," : "); put(file,z); new_line(file);
        Clear(z);
      end;
    end loop;
  end Write_Results;

  procedure Save_Results ( q : in Prod_Sys; sols : in Solution_List ) is

    file : file_type;

  begin
    if not Is_Null(sols) then
      new_line;
      put_line("Reading file name to write start system.");
      Read_Name_and_Create_File(file);
      put_line(file,natural32(q'last),q);
      new_line(file);
      put_line(file,"THE SOLUTIONS : ");
      new_line(file);
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Save_Results;

  procedure Main ( file : in file_type; p : in Poly_Sys;
                   b : in out natural32;
                   q : out Poly_Sys; rq : out Prod_Sys;
                   qsols : out Solution_List ) is

    gb : natural32 := b;

  begin
    Define_Partitions(file,p,gb,b);
    if not Degree_Structure.Empty
     then Construct_Start_System(file,p,q,rq,qsols);
    end if;
  end Main; 

  procedure Main is

    file : file_type;
    lp : Link_to_Poly_Sys;
    b : natural32 := 0;

  begin
    get(lp);
    declare
      q : Poly_Sys(lp'range);
      rq : Prod_Sys(q'range);
      qsols : Solution_List;
    begin
      put_line("Reading the output file.");
      Read_Name_and_Create_File(file);
      Main(file,lp.all,b,q,rq,qsols);
    end;
  end Main;

end Main_Multi_Homogenization;
