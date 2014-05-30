with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_VecVecs;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Polyhedral_Coefficient_Parameters;
with Polyhedral_Coefficient_Trackers;    use Polyhedral_Coefficient_Trackers;

package body Jumpstart_Polyhedral_Homotopies is

  function Mixed_Volume ( file : file_type; n,r : integer32;
                          mix : Standard_Integer_Vectors.Link_to_Vector )
                        return natural32 is

    m,mv,sum : natural32 := 0;
    mlb : Mixed_Labels;
    mic : Mixed_Cell;
    fail : boolean;
    ls : Standard_Floating_VecVecs.Array_of_VecVecs(1..r);

  begin
    Read_Lifted_Supports(file,natural32(n),natural32(r),ls,fail);
    if fail then
      put_line("Failed to read the lifted supports.  Please try again...");
      return 0;
    else
      get(file,m);
      put("Processing "); put(m,1); put(" mixed cells...");
      sum := 0;
      for i in 1..m loop
        Read_Next(file,natural32(n),natural32(r),i,mlb,fail);
        exit when fail;
        mic := Create_Coordinates(ls,mlb);
        Mixed_Volume(n,mix.all,mic,mv);
        sum := sum + mv;
        Clear(mlb); Deep_Clear(mic);
      end loop;
      return sum;
    end if;
  end Mixed_Volume;

  procedure Retrieve_Mixed_Volume
              ( file : in out file_type; n,r : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                mv : out natural32 ) is

    ans : character;
    nn,rr : natural32 := 0;
    nmix : Standard_Integer_Vectors.Link_to_Vector;
    fail : boolean;

  begin
    new_line;
    put("Do you know the mixed volume ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      mv := 0;
      put("Give the mixed volume : "); get(mv);
    else
      mv := Mixed_Volume(file,n,r,mix);
      put(" the mixed volume is "); put(mv,1); put_line(".");
      Reset(file);
      Read_Dimensions(file,nn,rr,nmix,fail);
    end if;
  end Retrieve_Mixed_Volume;

  procedure Jumpstart_Polyhedral_Continuation ( p : in Poly_Sys ) is

    infile,outfile : file_type;
    n,r : natural32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    fail,report,screen : boolean;
    q : Poly_Sys(p'range);
    ans : character;
    mv : natural32;

  begin
    put_line("Reading a file name for a regular mixed-cell configuration,");
    put_line("in a labeled representation...");
    Read_Name_and_Open_File(infile);
    Read_Dimensions(infile,n,r,mix,fail);
    if fail then
      put_line("Failed to read the dimensions of mixed-cell configuration.");
    else
      put_line("The dimension of the mixed-cell configuration :");
      Write_Dimensions(standard_output,n,r,mix.all);
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(outfile);
      Retrieve_Mixed_Volume(infile,integer32(n),integer32(r),mix,mv);
      new_line;
      Polyhedral_Coefficient_Parameters.Tune;
      new_line;
      put("Do you want reporting correctors during path tracking ? (y/n) ");
      Ask_Yes_or_No(ans);
      report := (ans = 'y');
      put("Do you want to monitor the progress on screen ? (y/n) ");
      ask_Yes_or_No(ans);
      screen := (ans = 'y');
      Polyhedral_Continuation
        (infile,outfile,report,screen,integer32(n),mv,mix.all,q);
    end if;
  end Jumpstart_Polyhedral_Continuation;

  procedure Read_Cells_and_Track 
              ( p : in Poly_Sys;
                cfile : in out file_type; ofile : in file_type ) is

    n,r : natural32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    fail,report,screen : boolean;
    q : Poly_Sys(p'range);
    ans : character;
    mv : natural32;

  begin
    Read_Dimensions(cfile,n,r,mix,fail);
    if fail then
      put_line("Failed to read the dimensions of mixed-cell configuration.");
    else
      put_line("The dimension of the mixed-cell configuration :");
      Write_Dimensions(standard_output,n,r,mix.all);
      Retrieve_Mixed_Volume(cfile,integer32(n),integer32(r),mix,mv);
      new_line;
      Polyhedral_Coefficient_Parameters.Tune;
      Polyhedral_Coefficient_Parameters.Write(ofile);
      new_line;
      put("Do you want reporting correctors during path tracking ? (y/n) ");
      Ask_Yes_or_No(ans);
      report := (ans = 'y');
      put("Do you want to monitor the progress on screen ? (y/n) ");
      ask_Yes_or_No(ans);
      screen := (ans = 'y');
      Polyhedral_Continuation
        (cfile,ofile,report,screen,integer32(n),mv,mix.all,q);
    end if;
  end Read_Cells_and_Track;

end Jumpstart_Polyhedral_Homotopies;
