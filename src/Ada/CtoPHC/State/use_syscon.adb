with text_io;                           use text_io;
with Symbol_Table_Interface;
with Standard_PolySys_Interface;
with Standard_LaurSys_Interface;
with DoblDobl_PolySys_Interface;
with DoblDobl_LaurSys_Interface;
with QuadDobl_PolySys_Interface;
with QuadDobl_LaurSys_Interface;
with Multprec_PolySys_Interface;
with Multprec_LaurSys_Interface;

function use_syscon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Symbol_Table_Interface;
    use Standard_PolySys_Interface;
    use Standard_LaurSys_Interface;
    use DoblDobl_PolySys_Interface;
    use DoblDobl_LaurSys_Interface;
    use QuadDobl_PolySys_Interface;
    use QuadDobl_LaurSys_Interface;
    use MultPrec_PolySys_Interface;
    use MultPrec_LaurSys_Interface;

  begin
    case job is
      when 0 => return Standard_PolySys_Read(vrblvl);
      when 1 => return Standard_PolySys_Write(vrblvl);
      when 2 => return Standard_PolySys_Get_Dimension(a,vrblvl);
      when 3 => return Standard_PolySys_Set_Dimension(a,vrblvl);
      when 4 => return Standard_PolySys_Size(a,vrblvl);
      when 5 => return Standard_PolySys_Get_Term(a,b,c,vrblvl);
      when 6 => return Standard_PolySys_Add_Term(a,b,c,vrblvl);
      when 7 => return Standard_PolySys_Clear(vrblvl);
      when 8 => return Standard_PolySys_Total_Degree(a,vrblvl);
      when 9 => return Symbol_Table_Clear(vrblvl);
      when 10 => return Standard_PolySys_Make_Function(vrblvl);
      when 11 => return Standard_PolySys_Jacobian_Function(vrblvl);
     -- dropping variables from polynomials
      when 12 => return Standard_PolySys_Drop_by_Index(a,vrblvl);
      when 13 => return DoblDobl_PolySys_Drop_by_Index(a,vrblvl);
      when 14 => return QuadDobl_PolySys_Drop_by_Index(a,vrblvl);
      when 15 => return Standard_PolySys_Drop_by_Name(a,b,vrblvl);
      when 16 => return DoblDobl_PolySys_Drop_by_Name(a,b,vrblvl);
      when 17 => return QuadDobl_PolySys_Drop_by_Name(a,b,vrblvl);
     -- degrees of polynomials :
      when 20 => return Standard_PolySys_Degree(a,b,vrblvl);
     -- dropping variables from Laurent polynomials
      when 22 => return Standard_LaurSys_Drop_by_Index(a,vrblvl);
      when 23 => return DoblDobl_LaurSys_Drop_by_Index(a,vrblvl);
      when 24 => return QuadDobl_LaurSys_Drop_by_Index(a,vrblvl);
      when 25 => return Standard_LaurSys_Drop_by_Name(a,b,vrblvl);
      when 26 => return DoblDobl_LaurSys_Drop_by_Name(a,b,vrblvl);
      when 27 => return QuadDobl_LaurSys_Drop_by_Name(a,b,vrblvl);
     -- jobs for standard double complex Laurent polynomials :
      when 100 => return Standard_LaurSys_Read(vrblvl);
      when 101 => return Standard_LaurSys_Write(vrblvl);
      when 102 => return Standard_LaurSys_Get_Dimension(a,vrblvl);
      when 103 => return Standard_LaurSys_Set_Dimension(a,vrblvl);
      when 104 => return Standard_LaurSys_Size(a,vrblvl);
      when 105 => return Standard_LaurSys_Get_Term(a,b,c,vrblvl);
      when 106 => return Standard_LaurSys_Add_Term(a,b,c,vrblvl);
      when 107 => return Standard_LaurSys_Clear(vrblvl);
     -- jobs for double double complex Laurent polynomials :
      when 110 => return DoblDobl_LaurSys_Read(vrblvl);
      when 111 => return DoblDobl_LaurSys_Write(vrblvl);
      when 112 => return DoblDobl_LaurSys_Get_Dimension(a,vrblvl);
      when 113 => return DoblDobl_LaurSys_Set_Dimension(a,vrblvl);
      when 114 => return DoblDobl_LaurSys_Size(a,vrblvl);
      when 115 => return DoblDobl_LaurSys_Get_Term(a,b,c,vrblvl);
      when 116 => return DoblDobl_LaurSys_Add_Term(a,b,c,vrblvl);
      when 117 => return DoblDobl_LaurSys_Clear(vrblvl);
      when 118 => return DoblDobl_LaurSys_String_Save(a,b,vrblvl);
     -- jobs for quad double complex Laurent polynomials :
      when 120 => return QuadDobl_LaurSys_Read(vrblvl);
      when 121 => return QuadDobl_LaurSys_Write(vrblvl);
      when 122 => return QuadDobl_LaurSys_Get_Dimension(a,vrblvl);
      when 123 => return QuadDobl_LaurSys_Set_Dimension(a,vrblvl);
      when 124 => return QuadDobl_LaurSys_Size(a,vrblvl);
      when 125 => return QuadDobl_LaurSys_Get_Term(a,b,c,vrblvl);
      when 126 => return QuadDobl_LaurSys_Add_Term(a,b,c,vrblvl);
      when 127 => return QuadDobl_LaurSys_Clear(vrblvl);
      when 128 => return QuadDobl_LaurSys_String_Save(a,b,vrblvl);
     -- jobs for multiprecision complex Laurent polynomials :
      when 130 => return Multprec_LaurSys_Read(vrblvl);
      when 131 => return Multprec_LaurSys_Write(vrblvl);
      when 132 => return Multprec_LaurSys_Get_Dimension(a,vrblvl);
      when 133 => return Multprec_LaurSys_Set_Dimension(a,vrblvl);
      when 134 => return Multprec_LaurSys_Size(a,vrblvl);
      when 137 => return Multprec_LaurSys_Clear(vrblvl);
      when 138 => return Multprec_LaurSys_String_Save(a,b,vrblvl);
      when 139 => return Multprec_LaurSys_String_Load(a,b,vrblvl);
     -- jobs for double double complex polynomials
      when 200 => return DoblDobl_PolySys_Read(vrblvl);
      when 201 => return DoblDobl_PolySys_Write(vrblvl);
      when 202 => return DoblDobl_PolySys_Get_Dimension(a,vrblvl);
      when 203 => return DoblDobl_PolySys_Set_Dimension(a,vrblvl);
      when 204 => return DoblDobl_PolySys_Size(a,vrblvl);
      when 205 => return DoblDobl_PolySys_Get_Term(a,b,c,vrblvl);
      when 206 => return DoblDobl_PolySys_Add_Term(a,b,c,vrblvl);
      when 207 => return DoblDobl_PolySys_Clear(vrblvl);
      when 208 => return DoblDobl_PolySys_String_Save(a,b,vrblvl);
      when 209 => return DoblDobl_PolySys_Degree(a,b,vrblvl);
     -- jobs for quad double complex polynomials
      when 210 => return QuadDobl_PolySys_Read(vrblvl);
      when 211 => return QuadDobl_PolySys_Write(vrblvl);
      when 212 => return QuadDobl_PolySys_Get_Dimension(a,vrblvl);
      when 213 => return QuadDobl_PolySys_Set_Dimension(a,vrblvl);
      when 214 => return QuadDobl_PolySys_Size(a,vrblvl);
      when 215 => return QuadDobl_PolySys_Get_Term(a,b,c,vrblvl);
      when 216 => return QuadDobl_PolySys_Add_Term(a,b,c,vrblvl);
      when 217 => return QuadDobl_PolySys_Clear(vrblvl);
      when 218 => return QuadDobl_PolySys_String_Save(a,b,vrblvl);
      when 219 => return QuadDobl_PolySys_Degree(a,b,vrblvl);
     -- jobs for multiprecision complex polynomials
      when 220 => return Multprec_PolySys_Read(vrblvl);
      when 221 => return Multprec_PolySys_Write(vrblvl);
      when 222 => return Multprec_PolySys_Get_Dimension(a,vrblvl);
      when 223 => return Multprec_PolySys_Set_Dimension(a,vrblvl);
      when 224 => return Multprec_PolySys_Size(a,vrblvl);
      when 227 => return Multprec_PolySys_Clear(vrblvl);
      when 228 => return Multprec_PolySys_String_Save(a,b,vrblvl);
      when 229 => return Multprec_PolySys_Degree(a,b,vrblvl);
     -- jobs for interchanging polynomial as strings :
      when 67 => return Standard_PolySys_String_Load(a,b,vrblvl);
      when 68 => return DoblDobl_PolySys_String_Load(a,b,vrblvl);
      when 69 => return QuadDobl_PolySys_String_Load(a,b,vrblvl);
      when 70 => return Multprec_PolySys_String_Load(a,b,vrblvl);
      when 71 => return Standard_PolySys_Random_System(a,b,vrblvl);
      when 72 => return DoblDobl_LaurSys_String_Load(a,b,vrblvl);
      when 73 => return QuadDobl_LaurSys_String_Load(a,b,vrblvl);
      when 74 => return Standard_LaurSys_String_Save(a,b,vrblvl);
      when 76 => return Standard_PolySys_String_Save(a,b,vrblvl);
      when 77 => return Standard_LaurSys_String_Load(a,b,vrblvl);
     -- random systems in double double and quad double precision
      when 78 => return DoblDobl_PolySys_Random_System(a,b,vrblvl);
      when 79 => return QuadDobl_PolySys_Random_System(a,b,vrblvl);
     -- jobs to return the size limit of the string representations
      when 80 => return Standard_PolySys_String_Size(a,b,vrblvl);
      when 81 => return DoblDobl_PolySys_String_Size(a,b,vrblvl);
      when 82 => return QuadDobl_PolySys_String_Size(a,b,vrblvl);
      when 83 => return Multprec_PolySys_String_Size(a,b,vrblvl);
      when 84 => return Standard_LaurSys_String_Size(a,b,vrblvl);
      when 85 => return DoblDobl_LaurSys_String_Size(a,b,vrblvl);
      when 86 => return QuadDobl_LaurSys_String_Size(a,b,vrblvl);
      when 87 => return Multprec_LaurSys_String_Size(a,b,vrblvl);
     -- reading systems into the containers :
      when 540 => return Standard_PolySys_Read_from_File(a,b,vrblvl); 
      when 541 => return DoblDobl_PolySys_Read_from_File(a,b,vrblvl); 
      when 542 => return QuadDobl_PolySys_Read_from_File(a,b,vrblvl); 
      when 543 => return Multprec_PolySys_Read_from_File(a,b,vrblvl); 
     -- projective transformations :
      when 891 => return Standard_PolySys_Make_Homogeneous(a,vrblvl);
      when 892 => return DoblDobl_PolySys_Make_Homogeneous(a,vrblvl);
      when 893 => return QuadDobl_PolySys_Make_Homogeneous(a,vrblvl);
      when 904 => return Standard_PolySys_Multi_Homogeneous(a,b,vrblvl);
      when 905 => return DoblDobl_PolySys_Multi_Homogeneous(a,b,vrblvl);
      when 906 => return QuadDobl_PolySys_Multi_Homogeneous(a,b,vrblvl);
     -- add symbol passed as string to the table
      when 897 => return Symbol_Table_Add(a,b,vrblvl);
     -- affine transformations :
      when 901 => return Standard_PolySys_1Hom2Affine(vrblvl);
      when 902 => return DoblDobl_PolySys_1Hom2Affine(vrblvl);
      when 903 => return QuadDobl_PolySys_1Hom2Affine(vrblvl);
      when 907 => return Standard_PolySys_mHom2Affine(a,vrblvl);
      when 908 => return DoblDobl_PolySys_mHom2Affine(a,vrblvl);
      when 909 => return QuadDobl_PolySys_mHom2Affine(a,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_syscon;
