with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Polyhedral_Coefficient_Trackers;    use Polyhedral_Coefficient_Trackers;
with Polyhedral_Coefficient_Parameters;
with Jumpstart_Polyhedral_Homotopies;    use Jumpstart_Polyhedral_Homotopies;

procedure ts_exptrack is

-- DESCRIPTION :
--   Incremental development of a path tracker for a homotopy whose
--   continuation parameter occurs as the argument of the exp function.

  procedure Coordinate_Polyhedral_Continuation ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   This procedure is a driver to new path trackers using the
  --   coordinate representations of mixed-cell configurations.

    infile,outfile : file_type;
    n,r,mv : natural32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sub : Mixed_Subdivision;
    q : Poly_Sys(p'range);
    ans : character;
    report,screen : boolean;

  begin
    put_line("Reading a file name for a regular mixed-cell configuration...");
    Read_Name_and_Open_File(infile);
    get(infile,n,r,mix,sub);
    new_line;
    put_line("Computing the volumes of all mixed cells...");
    Mixed_Volume(integer32(n),mix.all,sub,mv);
    put("The mixed volume is "); put(mv,1); put_line(".");
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(outfile);
    new_line;
    Polyhedral_Coefficient_Parameters.Tune;
    new_line;
    put("Do you want reporting correctors during path tracking ? (y/n) ");
    Ask_Yes_or_no(ans);
    report := (ans = 'y');
    put("Do you want to monitor the progress on screen ? (y/n) ");
    ask_Yes_or_No(ans);
    screen := (ans = 'y');
    Polyhedral_Continuation
      (outfile,report,screen,integer32(n),mv,p,mix.all,sub,q);
  end Coordinate_Polyhedral_Continuation;

  procedure Main is

    lp : Link_to_Poly_Sys;
    ans : character;
 
  begin
    new_line;
    put_line("Polyhedral homotopies with t = exp(s), s = -N..0");
    new_line;
    put_line("Reading a polynomial system...");
    get(lp);
    new_line;
    put_line("MENU for different representations of mixed subdivisions :");
    put_line("  1. a coordinate representation of a mixed-cell configuration");
    put_line("     does not permit incremental reading of the mixed cells;");
    put_line("  2. use labeled representation of a mixed-cell configuration");
    put_line("     avoids to keep all mixed cells in memory simultaneously.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    new_line;
    if ans = '1'
     then Coordinate_Polyhedral_Continuation(lp.all);
     else Jumpstart_Polyhedral_Continuation(lp.all);
    end if;
  end Main;

begin
  Main;
end ts_exptrack;
