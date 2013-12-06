with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Witness_Sets_io;                   use Witness_Sets_io;
with Intrinsic_Witness_Sets_io;         use Intrinsic_Witness_Sets_io;
with Drivers_to_Intersect_Varieties;    use Drivers_to_Intersect_Varieties;

procedure mainwit ( witset_one,witset_two,logfile : in string ) is

  procedure Read_Witness_Set
              ( w : in string; k : in natural32; p : out Link_to_Poly_Sys;
                sols : out Solution_List; dim : out natural32 ) is

  -- DESCRIPTION :
  --   Reads witness set k from file with name in the string w.

  -- ON ENTRY :
  --   w        name of file, if empty, the user is asked for a name;
  --   k        number of witness set, used for prompting the user.

  -- ON RETURN :
  --   p        polynomial equations defining the witness set;
  --   sols     points in the witness set;

    file : file_type;

  begin
    if w = "" then
      new_line;
      put("Reading a file name for witness set ");
      put(k,1); put_line(".");
      Read_Name_and_Open_File(file);
    else
      Open_Input_File(file,w);
    end if;
    Standard_Read_Embedding(file,p,sols,dim);
  end Read_Witness_Set;

  procedure Intersect_Witness_Sets
              ( file : in file_type; filename : in string;
                p1,p2 : in Poly_Sys; w1,w2 : in Solution_List;
                d1,d2 : in natural32 ) is

  -- DESCRIPTION :
  --   Intersects two witness sets of dimensions d1 and d2, d1 >= d2.

  -- ON ENTRY :
  --   file     log file to be opened for output;
  --   filename is the name of the output file, will be used to create
  --            files with witness sets of the intersection;
  --   p1,p2    polynomials defining the two solution components;
  --   w1,w2    witness points of the two solution components;
  --   d1,d2    dimensions of the two witness sets.

    n : constant natural32 := natural32(p1'last) - d1;
    f : Link_to_Poly_Sys;
    p : Link_to_Matrix;
    s : Solution_List;

  begin
   -- new_line;
   -- put("Intersecting witness sets of dimensions ");
   -- put(d1,1); put(" and "); put(d2,1);
   -- put(" in "); put(n,1); put_line("-space.");
   -- put_line("See the output file for results...");
   -- new_line;
    Intrinsic_Diagonal_Homotopy(file,false,p1,p2,w1,w2,d1,d2,f,p,s);
    Write_Witness_Set(file,filename,n,n-natural32(p'last(2)),f.all,s,p.all);
  end Intersect_Witness_Sets;

  procedure Main is

    p1,p2 : Link_to_Poly_Sys;
    w1,w2 : Solution_List;
    d1,d2 : natural32;
    file : file_type;
    name : Link_to_String;

  begin
    Read_Witness_Set(witset_one,1,p1,w1,d1);
    Read_Witness_Set(witset_two,2,p2,w2,d2);
    if logfile = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file,name);
    else
      Create_Output_File(file,logfile,name);
    end if;
    if d1 >= d2
     then Intersect_Witness_Sets(file,name.all,p1.all,p2.all,w1,w2,d1,d2);
     else Intersect_Witness_Sets(file,name.all,p2.all,p1.all,w2,w1,d2,d1);
    end if;
  end Main;

begin
  Main;
end mainwit;
