with text_io;                           use text_io;
with File_Management_Interface;
with Standard_Solutions_Interface;
with DoblDobl_Solutions_Interface;
with QuadDobl_Solutions_Interface;
with Multprec_Solutions_Interface;

function use_solcon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use File_Management_Interface;
    use Standard_Solutions_Interface;
    use DoblDobl_Solutions_Interface;
    use QuadDobl_Solutions_Interface;
    use Multprec_Solutions_Interface;

  begin
    case job is
      when 0 => return Standard_Solutions_Read(vrblvl);
      when 1 => return Standard_Solutions_Write(vrblvl);
      when 2 => return Standard_Solutions_Size(b,vrblvl);
      when 3 => return Standard_Solutions_Dimension(b,vrblvl);
      when 4 => return Standard_Solutions_Get(a,b,c,vrblvl);
      when 5 => return Standard_Solutions_Update(a,b,c,vrblvl);
      when 6 => return Standard_Solutions_Add(b,c,vrblvl);
      when 7 => return Standard_Solutions_Clear(vrblvl);
      when 8 => return Standard_Solutions_Drop_by_Index(a,vrblvl);
      when 9 => return Standard_Solutions_Drop_by_Name(a,b,vrblvl);
      when 10 => return File_Management_Prompt_Input_File(vrblvl);
      when 11 => return File_Management_Prompt_Output_File(vrblvl);
      when 12 => return Standard_Solutions_Scan_Banner(vrblvl);
      when 13 => return Standard_Solutions_Read_Dimensions(a,b,vrblvl);
      when 14 => return Standard_Solutions_Write_Dimensions(a,b,vrblvl);
      when 15 => return Standard_Solutions_Read_Next(a,b,c,vrblvl);
      when 16 => return Standard_Solutions_Write_Next(a,b,c,vrblvl);
      when 17 => return Standard_Solutions_Close_Input_File(a,vrblvl);
      when 18 => return Standard_Solutions_Close_Output_File(vrblvl);
      when 19 => return Standard_Solutions_Banner_to_Output(vrblvl);
      when 20 => return Standard_Solutions_Dimensions_to_Output(a,b,vrblvl);
      when 21 => return Standard_Solutions_Next_to_File(a,b,c,vrblvl);
      when 22 => return Standard_Solutions_Total_Degree(a,b,c,vrblvl);
      when 23 => return Standard_Solutions_Next_Product(a,b,c,vrblvl);
      when 24 => return Standard_Solutions_Lex_Product(a,b,c,vrblvl);
      when 25 => return Standard_Solutions_Next_Witness(a,b,c,vrblvl);
      when 30 => return Standard_Solutions_String_Size(a,b,vrblvl);
      when 31 => return Standard_Solutions_Get_String(a,b,vrblvl);
      when 32 => return Standard_Solutions_Intro_String_Size(a,b,vrblvl);
      when 33 => return Standard_Solutions_Intro_String(a,b,vrblvl);
      when 34 => return Standard_Solutions_Vector_String_Size(a,b,vrblvl);
      when 35 => return Standard_Solutions_Vector_String(a,b,vrblvl);
      when 36 => return Standard_Solutions_Diagnostics_String_Size(b,vrblvl);
      when 37 => return Standard_Solutions_Diagnostics_String(a,b,vrblvl);
      when 38 => return Standard_Solutions_Add_String(a,b,vrblvl);
      when 39 => return Standard_Solutions_Replace_String(a,b,vrblvl);
     -- corresponding operations for double double solutions
      when 40 => return DoblDobl_Solutions_Read(vrblvl);
      when 41 => return DoblDobl_Solutions_Write(vrblvl);
      when 42 => return DoblDobl_Solutions_Size(b,vrblvl);
      when 43 => return DoblDobl_Solutions_Dimension(b,vrblvl);
      when 44 => return DoblDobl_Solutions_Get(a,b,c,vrblvl);
      when 45 => return DoblDobl_Solutions_Update(a,b,c,vrblvl);
      when 46 => return DoblDobl_Solutions_Add(b,c,vrblvl);
      when 47 => return DoblDobl_Solutions_Clear(vrblvl);
      when 48 => return DoblDobl_Solutions_Drop_by_Index(a,vrblvl);
      when 49 => return DoblDobl_Solutions_Drop_by_Name(a,b,vrblvl);
      when 70 => return DoblDobl_Solutions_String_Size(a,b,vrblvl);
      when 71 => return DoblDobl_Solutions_Get_String(a,b,vrblvl);
      when 78 => return DoblDobl_Solutions_Add_String(a,b,vrblvl);
     -- corresponding operations for quad double solutions
      when 80 => return QuadDobl_Solutions_Read(vrblvl);
      when 81 => return QuadDobl_Solutions_Write(vrblvl);
      when 82 => return QuadDobl_Solutions_Size(b,vrblvl);
      when 83 => return QuadDobl_Solutions_Dimension(b,vrblvl);
      when 84 => return QuadDobl_Solutions_Get(a,b,c,vrblvl);
      when 85 => return QuadDobl_Solutions_Update(a,b,c,vrblvl);
      when 86 => return QuadDobl_Solutions_Add(b,c,vrblvl);
      when 87 => return QuadDobl_Solutions_Clear(vrblvl);
      when 88 => return QuadDobl_Solutions_Drop_by_Index(a,vrblvl);
      when 89 => return QuadDobl_Solutions_Drop_by_Name(a,b,vrblvl);
      when 110 => return QuadDobl_Solutions_String_Size(a,b,vrblvl);
      when 111 => return QuadDobl_Solutions_Get_String(a,b,vrblvl);
      when 118 => return QuadDobl_Solutions_Add_String(a,b,vrblvl);
     -- corresponding operations for multiprecision solutions
      when 120 => return Multprec_Solutions_Read(vrblvl);
      when 121 => return Multprec_Solutions_Write(vrblvl);
      when 122 => return Multprec_Solutions_Size(b,vrblvl);
      when 123 => return Multprec_Solutions_Dimension(b,vrblvl);
      when 127 => return Multprec_Solutions_Clear(vrblvl);
      when 150 => return Multprec_Solutions_String_Size(a,b,vrblvl);
      when 151 => return Multprec_Solutions_Get_String(a,b,vrblvl);
      when 158 => return Multprec_Solutions_Add_String(a,b,vrblvl);
     -- retrieve next solution
      when 276 => return Standard_Solutions_Retrieve_Next(a,b,c,vrblvl);
      when 277 => return DoblDobl_Solutions_Retrieve_Next(a,b,c,vrblvl);
      when 278 => return QuadDobl_Solutions_Retrieve_Next(a,b,c,vrblvl);
      when 279 => return Multprec_Solutions_Set_Pointer(vrblvl);
     -- move pointer from current to next solution
      when 300 => return Standard_Solutions_Move_Pointer(a,vrblvl);
      when 301 => return DoblDobl_Solutions_Move_Pointer(a,vrblvl);
      when 302 => return QuadDobl_Solutions_Move_Pointer(a,vrblvl);
      when 303 => return Multprec_Solutions_Move_Pointer(a,vrblvl);
     -- return length of current solution string
      when 304 => return Standard_Solutions_Current_Size(a,b,vrblvl);
      when 305 => return DoblDobl_Solutions_Current_Size(a,b,vrblvl);
      when 306 => return QuadDobl_Solutions_Current_Size(a,b,vrblvl);
      when 307 => return Multprec_Solutions_Current_Size(a,b,vrblvl);
     -- return current solution string
      when 308 => return Standard_Solutions_Current_String(a,b,vrblvl);
      when 309 => return DoblDobl_Solutions_Current_String(a,b,vrblvl);
      when 310 => return QuadDobl_Solutions_Current_String(a,b,vrblvl);
      when 311 => return Multprec_Solutions_Current_String(a,b,vrblvl);
     -- reading system and solutions from given file name
      when 544 => return Standard_System_Solutions_Read_from_File(a,b,vrblvl);
      when 545 => return DoblDobl_System_Solutions_Read_from_File(a,b,vrblvl);
      when 546 => return QuadDobl_System_Solutions_Read_from_File(a,b,vrblvl);
      when 547 => return Multprec_System_Solutions_Read_from_File(a,b,vrblvl);
     -- setting value of continuation parameter of the solutions to zero
      when 875 => return Standard_Solutions_Tzero(vrblvl);
      when 876 => return DoblDobl_Solutions_Tzero(vrblvl);
      when 877 => return QuadDobl_Solutions_Tzero(vrblvl);
     -- add one to each solution in 1-homogeneous coordinate transformation
      when 894 => return Standard_Solutions_Make_Homogeneous(vrblvl);
      when 895 => return DoblDobl_Solutions_Make_Homogeneous(vrblvl);
      when 896 => return QuadDobl_Solutions_Make_Homogeneous(vrblvl);
     -- add ones to each solution in m-homogeneous coordinate transformation
      when 910 => return Standard_Solutions_Multi_Homogeneous(a,vrblvl);
      when 911 => return DoblDobl_Solutions_Multi_Homogeneous(a,vrblvl);
      when 912 => return QuadDobl_Solutions_Multi_Homogeneous(a,vrblvl);
     -- affine transformations of the solutions
      when 898 => return Standard_Solutions_1Hom2Affine(vrblvl);
      when 899 => return DoblDobl_Solutions_1Hom2Affine(vrblvl);
      when 900 => return QuadDobl_Solutions_1Hom2Affine(vrblvl);
      when 913 => return Standard_Solutions_mHom2Affine(a,b,vrblvl);
      when 914 => return DoblDobl_Solutions_mHom2Affine(a,b,vrblvl);
      when 915 => return QuadDobl_Solutions_mHom2Affine(a,b,vrblvl);
      when 916 => return Standard_Solutions_Read_from_File(a,b,vrblvl);
      when 917 => return Dobldobl_Solutions_Read_from_File(a,b,vrblvl);
      when 918 => return QuadDobl_Solutions_Read_from_File(a,b,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_solcon;
