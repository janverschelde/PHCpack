with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Vectors;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;

procedure ts_mixlab is

-- DESCRIPTION :
--   A mixed cell can be stored in coordinate form, or using labels
--   with respect to a lifted configurations of points.
--   The coordinate form is self contained and perhaps most natural,
--   while the use of labels is more compact and better for global use.

  procedure Main is

    infile,outfile : file_type;
    n,r : natural32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    crd_sub : Mixed_Subdivision;
    lab_sub : Mixed_Sublabeling;
    ans : character;
 
  begin
    new_line;
    put_line("Coordinate and labeled representations of mixed cells.");
    new_line;
    put_line("MENU to convert representations of mixed-cell configurations :");
    put_line("  1. convert coordinate into labeled representation;");
    put_line("  2. convert labeled into coordinate representation.");
    put("Type 1 or 2 to make your choice : "); Ask_Alternative(ans,"12");
    new_line;
    put_line("Reading a file name for a regular mixed-cell configuration...");
    Read_Name_and_Open_File(infile);
    if ans = '1'
     then get(infile,n,r,mix,crd_sub);
     else get(infile,n,r,mix,lab_sub);
    end if;
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(outfile);
    if ans = '1'
     then lab_sub := Create_Labeled_Subdivision(integer32(r),crd_sub);
          put(outfile,n,mix.all,lab_sub);
     else crd_sub := Create_Coordinate_Subdivision(integer32(r),lab_sub);
          put(outfile,n,mix.all,crd_sub);
    end if;
  end Main;

begin
  Main;
end ts_mixlab;
