with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Linear_Product_System;
with Set_Structure,Set_Structure_io;     use Set_Structure;
with Degree_Sets_Tables;                 use Degree_Sets_Tables;
with Orbits_of_Solutions;                use Orbits_of_Solutions;
with Permutations,Symmetry_Group;        use Permutations,Symmetry_Group;
with Symmetry_Group_io;                  use Symmetry_Group_io;
with Symbolic_Symmetry_Group_io;         use Symbolic_Symmetry_Group_io;
with Drivers_for_Symmetry_Group_io;      use Drivers_for_Symmetry_Group_io;
with Symmetric_Set_Structure;            use Symmetric_Set_Structure;
with Equivariant_Polynomial_Systems;     use Equivariant_Polynomial_Systems;
with Linear_Symmetric_Reduction;         use Linear_Symmetric_Reduction;

package body Driver_for_Symmetric_Set_Structure is

  procedure Symmetric_Set_Structure_Info is

    i : array(1..5) of string(1..65);

  begin
    i(1):="  A symmetric generalized Bezout number is based on  a  symmetric";
    i(2):="supporting  set  structure  and  allows  to  exploit  permutation";
    i(3):="symmetries in the system.  The corresponding linear-product start";
    i(4):="system has the same symmetric structure, so that in the homotopy,";
    i(5):="only the generating solution paths need to be traced.            ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Symmetric_Set_Structure_Info;

  procedure Driver_for_Symmetric_Random_Product_Systems
                  ( file : in file_type; p : in Poly_Sys; q : out Poly_Sys;
                    qsols : out Solution_List; bs : in out natural32;
                    lpos : in out List ) is
  
    tol : constant double_float := 10.0**(-12);
 
    procedure Write_Results ( file : in file_type; bb : in natural32 ) is
    begin
      new_line(file);
      put(file,"  generalized Bezout number is "); put(file,bb,1);
      new_line(file);
      put_line(file,"  based on the set structure :");
      Set_Structure_io.put(file);
    end Write_Results;

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

    procedure Write_Orbits
                ( file : in file_type; qqsols : in Solution_List ) is

      orb : constant Permutation := Orbits(qqsols,tol);

    begin
      new_line;
      put("The number of generating solutions : ");
      put(Length_Of(qqsols),1); new_line;
      new_line(file);
      put(file,"The number of generating solutions : ");
      put(file,Length_Of(qqsols),1); new_line(file);
      put("The orbits : "); Symmetry_Group_io.put(orb); new_line;
      put(file,"The orbits : "); Symmetry_Group_io.put(file,orb);
      new_line(file);
    end Write_Orbits;

    procedure Driver_for_Bezout_Number ( file : in file_type )  is

      timer : timing_widget;
      ns : Standard_Natural_Vectors.Vector(p'range);

    begin
      put_line("Reading the set structure.");
      for i in ns'range loop
        put("  Give the number of sets for polynomial ");
        put(i,1); put(" : ");
        Read_Natural(ns(i));
      end loop;
      Set_Structure.Init(ns);
      put_line("Give the set structure : ");
      Set_Structure_io.get;
     -- Set_Structure.B(bs,lpos);
      tstart(timer);
      bs := natural32(Permanent(Degree_Sets_Tables.Create));
      tstop(timer);
      Write_Results(file,bs);
      Write_Results(Standard_Output,bs);
      new_line(file);
      print_times(file,timer,"computation of generalized permanent");
    end Driver_for_Bezout_Number;

    procedure Construct_Start_System 
                 ( file : in file_type; n : in natural32;
                   allperms : in boolean; v,w : List_of_Permutations;
                   notsymmetric,degenerate : out boolean ) is

      timer : timing_widget;
      notequi,notsym,degen : boolean;

    begin
      tstart(timer);
      if allperms then
        Equivariant_Start_System(n,v,notequi);
        if notequi then
          new_line; new_line(file);
          put_line("The set structure is not equivariant.");
          put_line(file,"The set structure is not equivariant.");
        else
          notsym := false; degen := false;
        end if;
      end if;
      if not allperms or notequi then
        Symmetric_Start_System(n,bs,lpos,v,w,notsym,degen);
        new_line; new_line(file);
        if notsym then
          put_line("The set structure is not symmetric.");
          put_line(file,"The set structure is not symmetric.");
        else
          if degen then
            put_line("The set structure is symmetric but degenerate.");
            put_line(file,"The set structure is symmetric but degenerate.");
          else
            put_line("The set structure is symmetric and not degenerate.");
            put_line(file,"The set structure is symmetric and not degenerate.");
          end if;
        end if;
      end if;
      notsymmetric := notsym;
      degenerate := degen;
      tstop(timer);
      new_line(file);
      print_times(file,timer,"construction of symmetric start system");
    end Construct_Start_System;

    procedure Solve_Start_System 
                 ( file : in file_type;
                   allperms : in boolean; v,w : in List_of_Permutations ) is

      timer : timing_widget;
      nl : natural32;
      qq : Poly_Sys(p'range);
      qqsols : Solution_List;

    begin
     -- Standard_Linear_Product_System_io.put(file,n,2,4,3);
      qq := Standard_Linear_Product_System.Polynomial_System;
      new_line(file);
      put_line(file,"SYMMETRIC LINEAR-PRODUCT SYSTEM : ");
      put_line(file,qq);
     -- put_line(file,"The list of positions : "); put(file,lpos);
     -- if allperms
     --  then Linear_Symmetric_Reduce(lpos,false);
     --  else Linear_Symmetric_Reduce(v,w,lpos);
     -- end if;
      tstart(timer);
      if allperms
       then lpos := Linear_Symmetric_Reduce(false);
       else lpos := Linear_Symmetric_Reduce(v,w);
      end if;
     -- put_line(file,"The reduced list of positions : "); put(file,lpos);
      Standard_Linear_Product_System.Solve(qqsols,nl,lpos);
      tstop(timer);
      Standard_Linear_Product_System.Clear;
      if allperms
       then qqsols := Generating(qqsols,false,tol);
       else Analyze(v,false,tol,qqsols);
      end if;
      Save_Results(qq,qqsols);
      new_line(file);
      put_line(file,"THE GENERATING SOLUTIONS :");
      new_line(file);
      put(file,Length_Of(qqsols),natural32(Head_Of(qqsols).n),qqsols);
      Write_Orbits(file,qqsols);
      new_line(file);
      print_times(file,timer,"solving the linear-product system");
      q := qq; qsols := qqsols;
    end Solve_Start_System;

    procedure Driver_for_Start_System
                 ( file : in file_type; n : in natural32;
                   allperms : in boolean; v,w : List_of_Permutations ) is

      ans : character;
      notsym,degen : boolean;

    begin
      new_line;
      put("Do you want a symmetric linear-product start system ? ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Construct_Start_System(file,n,allperms,v,w,notsym,degen);
        -- new_line(file); Write_Covering(file); new_line(file);
        -- new_line(file); Write_Templates(file,n); new_line(file);
        --  Symmetric_Set_Structure.Clear;
        Set_Structure.Clear;
        if not notsym and not degen
         then Solve_Start_System(file,allperms,v,w);
        end if;
      end if;
    end Driver_for_Start_System;

    procedure Main_Driver is

      totaltimer : timing_widget;
      n : constant natural32 := natural32(p'length);
      allperms,notsym,inva,equi : boolean;
      g,v,w : List_of_Permutations;

    begin
      new_line(file);
      put_line(file,"SYMMETRIC SET STRUCTURE ANALYSIS :");
      new_line(file);
      Read_Permutation_Group(n,g,v,allperms);
      tstart(totaltimer);
      put_line(file,"THE SYMMETRY GROUP :");
      new_line(file);
      put_line(file,"v:"); Symbolic_Symmetry_Group_io.put(file,v);
      new_line(file);
      Act(v,p,w,notsym,inva,equi);
      new_line(file);
      put_line(file,"w:"); Symmetry_Group_io.put(file,w); new_line(file);
      if notsym then
        put_line("The system is not (G,V,W)-symmetric.");
        put_line(file,"The system is not (G,V,W)-symmetric.");
      else
        put_line("The system is (G,V,W)-symmetric.");
        put_line(file,"The system is (G,V,W)-symmetric.");
        if Set_Structure.Empty
         then Driver_for_Bezout_Number(file);
        end if;
        if not Set_Structure.Empty
         then Driver_for_Start_System(file,n,allperms,v,w);
        end if;
      end if;
      tstop(totaltimer);
      new_line(file);
      print_times(file,totaltimer,"symmetric set structure analysis");
    end Main_Driver;

  begin
    Main_Driver;
  end Driver_for_Symmetric_Random_Product_Systems;

end Driver_for_Symmetric_Set_Structure;
