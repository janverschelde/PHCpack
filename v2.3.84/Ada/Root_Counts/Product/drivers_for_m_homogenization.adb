with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Numbers_io;                         use Numbers_io;
with Timing_Package;                     use Timing_Package;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with m_Homogeneous_Bezout_Numbers;       use m_Homogeneous_Bezout_Numbers;
with m_Homogeneous_Start_Systems;        use m_Homogeneous_Start_Systems;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Prod_Systems_io;   use Standard_Complex_Prod_Systems_io;
with Standard_Complex_Prod_Planes;
with Interpolating_Homotopies_Driver;

package body Drivers_for_m_Homogenization is

  procedure m_Homogenization_Info is

    i : array(1..15) of string(1..65);

  begin
    i( 1):="  An m-homogeneous Bezout number is based on a partition  of  the";
    i( 2):="set of unknowns, with m the number of sets in the partition.     ";
    i( 3):="  The corresponding start system is a linear-product system:  the";
    i( 4):="i-th  equation  is  the  product  of linear equations with random";
    i( 5):="coefficients in the unknowns of the set of  the  partition.   The";
    i( 6):="number  of  factors  in the product of the i-th polynomial equals";
    i( 7):="the product of the degrees of the i-th polynomial in the original";
    i( 8):="system w.r.t. every set in the partition.                        ";
    i( 9):="  Given a partition, the m-homogeneous Bezout number  equals  the";
    i(10):="number  of  solutions  of  the corresponding linear-product start";
    i(11):="system.  Before the construction  of  the  start  system,  an  m-";
    i(12):="homogeneous  Bezout number is first computed in a formal way as a";
    i(13):="generalized permanent of a degree matrix.  Either all  partitions";
    i(14):="of  the  set  of  unknowns  can  be  evaluated,  or the available";
    i(15):="heuristic procedure for generating a partition can be applied.   ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end m_Homogenization_Info;

  procedure Driver_for_m_Homogenization
               ( file : in file_type; p : in Poly_Sys; b : in out natural64;
                 q : out Poly_Sys; qsols : out Solution_List ) is

    n : constant natural32 := natural32(p'length);
    np : constant natural32 := Number_of_Partitions(n);

    procedure Write_Results ( file : in file_type; bz : in natural64;
                              nz : in natural32; z : in partition  ) is
    begin
      new_line(file);
      put(file,"  "); put(file,nz,1);
      put(file,"-homogeneous Bezout number is "); put(file,bz,1);
      new_line(file);
      put(file,"     with partition : "); put(file,z); new_line(file);
    end Write_Results;

    procedure Save_Results ( rq : in Prod_Sys; m : in natural32 ) is

      rqfile : file_type;

    begin
      new_line;
      put_line("Reading file name to write linear-product start system.");
      Read_Name_and_Create_File(rqfile);
      put(rqfile,natural32(rq'last),1); new_line(rqfile);
      put_line(rqfile,rq);
      new_line(rqfile);
      put(rqfile,"TITLE : ");
      put(rqfile,m,1);
      put_line(rqfile,"-homogeneous linear-product start system");
      close(rqfile);
    end Save_Results;

    procedure Save_Results ( qq : in Poly_Sys; qqsols : in Solution_List ) is

      qqfile : file_type;

    begin
      if not Is_Null(qqsols) then
        new_line;
        put_line("Reading file name to write start system.");
        Read_Name_and_Create_File(qqfile);
        put_line(qqfile,qq);
        new_line(qqfile);
        put_line(qqfile,"THE SOLUTIONS : ");
        new_line(qqfile);
        put(qqfile,Length_Of(qqsols),natural32(Head_Of(qqsols).n),qqsols);
        Close(qqfile);
      end if;
    end Save_Results;

    function Universe ( n : natural32 ) return Set is

    -- DESCRIPTION : Returns the set of all unknowns.

      s : Set := Create(n);

    begin
      for i in 1..n loop
        Add(s,i);
      end loop;
      return s;
    end Universe;

    procedure Display_Menu
                ( choice : out character; bz : natural64;
                  nz : in natural32; z : in Partition ) is

      ans : character;

    begin
      new_line;
      put_line("MENU for computing m-Homogeneous Bezout Numbers :");
      put     ("  0. exit - Bezout number : ");
        put(bz,1); put(", partition : "); put(z(1..nz)); new_line;
      put_line("  1. Heuristic generation of a partition");
      put     ("  2. Enumerate all "); put(np,1); put_line(" partitions");
      put_line("  3. Enumerate over #partitions <= given maximum");
      put_line("  4. Enumerate till Bezout number <= given minimum");
      put_line("  5. Combine stop criteria of strategies 3 and 4");
      put_line("  6. Evaluate your own partition");
      put("Type number between 0 and 6 to choose : ");
      Ask_Alternative(ans,"0123456"); choice := ans;
    end Display_Menu;

    procedure Patch ( z : in out Partition; nz : in natural32;
                      bz : in out natural64 ) is

    -- DESCRIPTION :
    --   Because of this bug :
    --    3 
    --     x**2 + x*y + z;
    --     y**2 + x*y + z;
    --     x + y + z;
    --   The heuristic partitioner only returns {z} instead of {x y z}.

    -- NOTE :
    --   Declaring type of bz as simply "out" instead of "in out"
    --   causes the heuristic partitioner to fail.

    begin
      if nz = 1 then
        if Extent_Of(z(1)) < natural32(p'last) then
          Clear(z(1));
          z(1) := Universe(natural32(p'last));
          bz := Total_Degree(p);
        end if;
      end if;
    end Patch;

    procedure Dispatch_Menu
                 ( file : in file_type; choice : in character;
                   bz : in out natural64;
                   nz : in out natural32; z : in out Partition ) is

    -- DESCRIPTION :
    --   Dispatches the choice of the user.

    -- ON ENTRY :
    --   file      must be opened for output;
    --   choice    choice of method;
    --   bz        current Bezout number with partition z;
    --   nz        number of sets in current partition z;
    --   z         current partition is z(1..nz).

      max,min : natural32;

    begin
      new_line(file);
      case choice is
        when '1' => Clear(z); nz := 0;
                    put_line(file,"HEURISTIC PARTITIONER :");
                    PB(p,bz,nz,z);
                    Patch(z,nz,bz);
        when '2' => put(file,"ENUMERATION OF ");
                    put(file,np,1); put_line(file," PARTITIONS :");
                    Bezout_number(p,bz,nz,z);
        when '3' => put("  Give maximum bound : "); Read_Natural(max);
                    put(file,"ENUMERATION OF ");
                    put(file,max,1); put_line(file," PARTITIONS :");
                    Bezout_number(natural64(max),p,bz,nz,z);
        when '4' => put("  Give minimum bound : "); Read_Natural(min);
                    put(file,"ENUMERATION TILL BEZOUT NUMBER <= ");
                    put(file,min,1); new_line(file);
                    Bezout_number(p,natural64(min),bz,nz,z);
        when '5' => put("  Give maximum bound : "); Read_Natural(max);
                    put("  Give minimum bound : "); Read_Natural(min);
                    put(file,"ENUMERATION OF ");
                    put(file,max,1); put_line(file," PARTITIONS");
                    put(file,"         OR TILL BEZOUT NUMBER <= ");
                    put(file,min,1); new_line(file);
                    Bezout_number(natural64(max),p,natural64(min),bz,nz,z);
        when '6' => Clear(z); nz := 0;
                    put_line(file,"PARTITION PROVIDED BY USER :");
                    put("Give the number of sets : ");
                    Read_Positive(integer(nz));
                    put("Give "); put(nz,1); put(" sets : ");
                    Create(z,n); get(z(1..nz));
                    bz := Bezout_Number(p,z);
        when others => null;
      end case;
      Write_Results(file,bz,nz,z(1..nz));
      Write_Results(Standard_Output,bz,nz,z(1..nz));
    end Dispatch_Menu;

    procedure Driver_for_Partition
                 ( file : in file_type; bz : in out natural64;
                   nz : in out natural32; z : in out Partition ) is

      method : character;
      timer : timing_widget;

    begin
      new_line(file);
      put_line(file,"M-HOMOGENIZATION :");
      tstart(timer);
      loop
        Display_Menu(method,bz,nz,z);
        exit when method = '0';
        Dispatch_Menu(file,method,bz,nz,z);
      end loop;
      tstop(timer);
      new_line(file);
      print_times(file,timer,"Computation of m-homogeneous Bezout number");
    end Driver_for_Partition;

    procedure Driver_for_Start_System
                 ( file : in file_type;
                   bz : in natural64; z : in Partition ) is

      ans : character;
      timer : timing_widget;
      qq : Poly_Sys(p'range);
      rq : Prod_Sys(p'range);
      qqsols : Solution_List;
      bb : natural64 := bz;
      save_rq : boolean := false;

    begin
      new_line;
      put_line("MENU for m-Homogeneous Start Systems :");
      put_line("  0. No construction of m-homogeneous start system.");
      put_line("  1. Start system based on interpolation.");
      put_line("  2. Random linear-product start system.");
      put("Type 0, 1, or 2 to choose : "); Ask_Alternative(ans,"012");
      if ans /= '0' then
        tstart(timer);
        if ans = '1' then
          Interpolating_Homotopies_Driver(file,p,z,natural32(bb),qq,qqsols);
        else
          new_line;
          put("Do you want to have all start solutions now ? (y/n) ");
          Ask_Yes_or_No(ans);
          if ans = 'y' then
            m_Homogeneous_Start_System(p,z,qq,qqsols);
          else
            m_Homogeneous_Start_System(p,z);
            rq := Standard_Complex_Prod_Planes.Create;
            save_rq := true;
          end if;
        end if;
        tstop(timer);
        if not save_rq then
          Save_Results(qq,qqsols);
          q := qq; qsols := qqsols;
          new_line(file);
          put(file,z'last,1);
          put_line(file,"-HOMOGENEOUS START SYSTEM : ");
          put_line(file,qq);
          if not Is_Null(qsols) then
            new_line(file);
            put_line(file,"THE SOLUTIONS :");
            new_line(file);
            put(file,Length_Of(qqsols),natural32(Head_Of(qqsols).n),qqsols);
          end if;
        else
          Save_Results(rq,z'last);
          new_line(file);
          put(file,z'last,1);
          put(file,"-HOMOGENOUS LINEAR-PRODUCT SYSTEM : ");
          put(file,natural32(rq'last),1); new_line(file);
          put_line(file,rq);
        end if;
        new_line(file);
        print_times(file,timer,
                    "Construction of m-homogeneous start system");
      end if;
    end Driver_for_Start_System;

    procedure Main_Driver is

      z : Partition(1..n);
      nz : natural32;
      bz : natural64;

    begin
      nz := 1; z(1) := Universe(n); bz := Total_Degree(p);
      Driver_for_Partition(file,bz,nz,z);
      Driver_for_Start_System(file,bz,z(1..nz));
      b := bz;
    end Main_Driver;

  begin
    Main_Driver;
  end Driver_for_m_Homogenization;

end Drivers_for_m_Homogenization;
