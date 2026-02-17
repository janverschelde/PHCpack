with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Communications_with_User;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;
with Real_Powered_Homotopy;
with Real_Powered_Homotopy_IO;

package body Main_Real_Powered_Series is

  procedure main ( vrblvl : in integer32 := 0 ) is

    infile,outfile : file_type;

  begin
    new_line;
    put_line("Reading the name of an input file ...");
    Communications_with_User.Read_Name_and_Open_File(infile);
    declare
      lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
      lc : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
      lp : Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
      lpv : Standard_Floating_VecVecs.Link_to_VecVec;
      lpv1 : Standard_Floating_Vectors.Link_to_Vector;
      dim,size : integer32;
    begin
      Real_Powered_Homotopy_IO.get(infile,lq,lc,lp,vrblvl=>vrblvl-1);
      close(infile);
      dim := lq'last;
      lpv := lp(lp'first);
      lpv1 := lpv(lpv'first);
      size := lpv1'last;
      if vrblvl > 0 then
        new_line;
        put("-> read "); put(dim,1); put(" polynomials");
        put(" with coefficients of size "); put(size,1);
        new_line;
        new_line;
        put_line("The system read from file :");
        Real_Powered_Homotopy_IO.put_line(dim,dim,size,lq.all,lc.all,lp.all);
      end if;
      new_line;
      put_line("Reading the name of an output file ...");     
      Communications_with_User.Read_Name_and_Create_File(outfile);
      Real_Powered_Homotopy_IO.put_line
        (outfile,dim,dim,size,lq.all,lc.all,lp.all);
      new_line;
      new_line(outfile);
      if Real_Powered_Homotopy.Is_Linear(lq.all) then
        put_line("The system is linear.");
        put_line(outfile,"The system is linear.");
      else
        put_line("The system is not linear.");
        put_line(outfile,"The system is not linear.");
      end if;
      close(outfile);
    end;
  end main;

end Main_Real_Powered_Series;
