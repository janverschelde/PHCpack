with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Symbol_Table_io;

package body Symbol_Table_Order is

  procedure Main ( infilename,outfilename : in string;
                   vrblvl : in integer32 := 0 ) is

    infile,outfile : file_type;
    lp : Link_to_Laur_Sys;

  begin
    if vrblvl > 0
     then put_line("-> in symbol_table_order.Main ...");
    end if;
    if infilename = "" then
      new_line;
      get(lp);
    else
      Open_Input_File(infile,infilename);
      get(infile,lp);
    end if;
    if outfilename = "" then
      Symbol_Table_io.Write;
    else
      Create_Output_File(outfile,outfilename);
      Symbol_Table_io.Write(outfile);
      close(outfile);
    end if;
  end Main;

end Symbol_Table_Order;
