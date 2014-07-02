with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;

procedure ts_jpmix is

-- DESCRIPTION :
--   This procedure interactively develops the processing of mixed cells
--   in a subdivision via the jumpstarting principle.

  procedure Process_Coordinate_Subdivision is

  -- DESCRIPTION :
  --   Processing a coordinate subdivision incrementally from file works 
  --   only when the operation is local because the configuration of
  --   lifted points is not stored explicitly at the front of the file,
  --   but must be extracted from scanning through the whole file.

    file : file_type;
    n,r,m,sum : natural32 := 0;
    mv : natural32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mic : Mixed_Cell;
    fail : boolean;

  begin
    put_line("Reading a file name for a coordinate subdivision...");
    Read_Name_and_Open_File(file);
    Read_Dimensions(file,n,r,m,mix,fail);
    new_line;
    if fail then
      put_line("Failed to read the dimensions.  Please check the file.");
    else
      put_line("The dimensions of the mixed-cell configuration : ");
      Write_Dimensions(standard_output,n,r,m,mix.all);
      put("Ready to read "); put(m,1);
      put_line(" mixed cells from file...");
      sum := 0;
      for i in 1..m loop
        Read_Next(file,n,r,mic,fail);
        exit when fail;
        put("Cell "); put(i,1); put(" has mixed volume ");
        Mixed_Volume(integer32(n),mix.all,mic,mv);
        put(mv,1);
        sum := sum + mv;
        put(", sum = "); put(sum,1); new_line;
        Deep_Clear(mic);
      end loop;
      put("Mixed volume : "); put(sum,1); new_line;
    end if;
  end Process_Coordinate_Subdivision;

  procedure Process_Labeled_Subdivision is

  -- DESCRIPTION :
  --   Processing a labeled mixed-cell configuration incrementally from
  --   file is more suitable because the lifted points configuration
  --   occurs right at the beginning, not that it matters for the local
  --   mixed volume computation, but this routine is a template for
  --   memory efficient polyhedral continuation.

    file : file_type;
    n,r,m,sum : natural32 := 0;
    mv : natural32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mlb : Mixed_Labels;
    mic : Mixed_Cell;
    fail : boolean;

  begin
    put_line("Reading a file name for a labeled subdivision...");
    Read_Name_and_Open_File(file);
    Read_Dimensions(file,n,r,mix,fail);
    new_line;
    if fail then
      put_line("Failed to read the dimensions.  Please check the file.");
    else
      put_line("The dimensions of the mixed-cell configuration : ");
      Write_Dimensions(standard_output,n,r,mix.all);
      declare
        ls : Standard_Floating_VecVecs.Array_of_VecVecs(1..integer32(r));
      begin
        Read_Lifted_Supports(file,n,r,ls,fail);
        if fail then
          put_line("Failed to read the lifted supports.  Please try again...");
        else
          put_line("The lifted supports of the mixed-cell configuration :");
          Write_Lifted_Supports(standard_output,ls);
          get(file,m);
          put("Reading and processing "); put(m,1);
          put_line(" mixed cells...");
          sum := 0;
          for i in 1..m loop
            Read_Next(file,n,r,i,mlb,fail);
            exit when fail;
            put("mixed cell "); put(i,1);
           -- put(standard_output,n,mix.all,mlb);
            mic := Create_Coordinates(ls,mlb);
            Mixed_Volume(integer32(n),mix.all,mic,mv);
            put(" has mixed volume "); put(mv,1);
            sum := sum + mv;
            put(", sum = "); put(sum,1); new_line;
            Clear(mlb); Deep_Clear(mic);
          end loop;
          put("Mixed volume : "); put(sum,1); new_line;
        end if;
      end;
    end if;
  end Process_Labeled_Subdivision;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for jumpstarting computation of mixed volume :");
    put_line("  1. coordinate representation of a mixed-cell configuration;");
    put_line("  2. use labeled representation of a mixed-cell configuration.");
    put("Type 1 or 2 to choose representation : ");
    Ask_Alternative(ans,"12");
    new_line;
    if ans = '1'
     then Process_Coordinate_Subdivision;
     else Process_Labeled_Subdivision;
    end if;
  end Main;

begin
  Main;
end ts_jpmix;
