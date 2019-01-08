with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Floating_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_System_and_Solutions_io;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Prod_Systems_io;   use Standard_Complex_Prod_Systems_io;
with Standard_Complex_Prod_Planes;
with PHCpack_Operations;
with File_Management;
with Standard_Solutions_Container;
with Witness_Sets_io;
with Extrinsic_Diagonal_Homotopies;
with Extrinsic_Diagonal_Homotopies_io;
with Jumpstart_Diagonal_Homotopies;

package body PHCpack_Operations_io is

  procedure Read_Start_System is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the start system...");
    Standard_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Start_System(p.all);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_Start_System;

  procedure Read_Start_Laurent_System is

    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the start system...");
    Standard_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Start_System(p.all);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_Start_Laurent_System;

  procedure Read_DoblDobl_Start_System is

    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the start system...");
    DoblDobl_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Start_System(p.all);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_DoblDobl_Start_System;

  procedure Read_DoblDobl_Start_Laurent_System is

    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the start system...");
    DoblDobl_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Start_System(p.all);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_DoblDobl_Start_Laurent_System;

  procedure Read_QuadDobl_Start_System is

    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the start system...");
    QuadDobl_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Start_System(p.all);
    if not QuadDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_QuadDobl_Start_System;

  procedure Read_QuadDobl_Start_Laurent_System is

    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the start system...");
    QuadDobl_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Start_System(p.all);
    if not QuadDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_QuadDobl_Start_Laurent_System;

  procedure Read_Multprec_Start_System ( decimals : in natural32 ) is

    p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Multprec_Complex_Solutions.Solution_List;
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(decimals);

  begin
    Multprec_Complex_Polynomials_io.Set_Working_Precision(size);
    new_line;
    put_line("Reading the start system...");
    Multprec_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Start_System(p.all);
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      Multprec_Complex_Solutions.Set_Size(sols,size);
      PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_Multprec_Start_System;

  procedure Read_Start_System ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    Standard_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Start_System(p.all);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_Start_System;

  procedure Read_Start_Laurent_System ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    Standard_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Start_System(p.all);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_Start_Laurent_System;

  procedure Read_DoblDobl_Start_System ( filename : in string ) is

    file : file_type;
    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    DoblDobl_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Start_System(p.all);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_DoblDobl_Start_System;

  procedure Read_DoblDobl_Start_Laurent_System ( filename : in string ) is

    file : file_type;
    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    DoblDobl_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Start_System(p.all);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_DoblDobl_Start_Laurent_System;

  procedure Read_QuadDobl_Start_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    QuadDobl_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Start_System(p.all);
    if not QuadDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_QuadDobl_Start_System;

  procedure Read_QuadDobl_Start_Laurent_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    QuadDobl_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Start_System(p.all);
    if not QuadDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_QuadDobl_Start_Laurent_System;

  procedure Read_Multprec_Start_System
              ( filename : in string; decimals : in natural32 ) is

    file : file_type;
    p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Multprec_Complex_Solutions.Solution_List;
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(decimals);

  begin
    Multprec_Complex_Polynomials_io.Set_Working_Precision(size);
    Open(file,in_file,filename);
    Multprec_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Start_System(p.all);
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      Multprec_Complex_Solutions.Set_Size(sols,size);
      PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_Multprec_Start_System;

  procedure Read_Start_System_without_Solutions is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading the name of the file for the start system.");
    File_Management.Silent_Open_Input_File;
    get(File_Management.Link_to_Input.all,p);
    PHCpack_Operations.Store_Start_System(p.all);
  end Read_Start_System_without_Solutions;

  procedure Read_Start_System_without_Solutions ( filename : in string ) is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    File_Management.Silent_Open_Input_File(filename);
    get(File_Management.Link_to_Input.all,p);
    PHCpack_Operations.Store_Start_System(p.all);
  end Read_Start_System_without_Solutions;

  procedure Read_Linear_Product_Start_System ( fail : out boolean ) is

    file : file_type;
    lq : Link_to_Prod_Sys;

  begin
    put_line
      ("Reading the name of a file for a linear-product start system...");
    Read_Name_and_Open_File(file);
    get(file,lq);
    Close(file);
    Standard_Complex_Prod_Planes.Store(lq.all,fail);
    if fail then
      put_line("Storing the system as a linear-product system failed!");
    else
      declare
        p : constant Standard_Complex_Poly_Systems.Poly_Sys(lq'range)
          := Expand(lq.all);
      begin
        PHCpack_Operations.Store_Start_System(p);
      end;
    end if;
  exception
    when others =>
      put_line("Exception raised when reading linear-product start system.");
      fail := true;
  end Read_Linear_Product_Start_System;

  procedure Read_Linear_Product_Start_System
              ( filename : in string; fail : out boolean ) is


    file : file_type;
    lq : Link_to_Prod_Sys;

  begin
    Open(file,in_file,filename);
    get(file,lq);
    Close(file);
    Standard_Complex_Prod_Planes.Store(lq.all,fail);
    if fail then
      put_line("Storing the system as a linear-product system failed!");
    else
      declare
        p : constant Standard_Complex_Poly_Systems.Poly_Sys(lq'range)
          := Expand(lq.all);
      begin
        PHCpack_Operations.Store_Start_System(p);
      end;
    end if;
  exception
    when others =>
      put_line("Exception raised when reading linear-product start system.");
      fail := true;
  end Read_Linear_Product_Start_System;

  procedure Product_of_Symbol_Tables
              ( n1,n2,a,b : in natural32; 
                sbt : in Symbol_Table.Array_of_Symbols ) is

  -- DESCRIPTION :
  --   Adds suffixes '1' and '2' to the symbols in sbt to display
  --   the target and start system to start the cascade to intersect
  --   solution sets of dimensions a and b, with ambient dimensions
  --   respectively equal to n1 and n2.

    use Symbol_Table;
    use Extrinsic_Diagonal_Homotopies_io;

    s : constant Array_of_Symbols := Remove_Embed_Symbols(sbt);
    s1 : constant Array_of_Symbols(s'range) := Add_Suffix(s,'1');
    s2 : constant Array_of_Symbols(s'range) := Add_Suffix(s,'2');
    cd : constant natural32
       := Extrinsic_Diagonal_Homotopies.Cascade_Dimension(n1,n2,a,b);

  begin
    Symbol_Table.Clear;
    Assign_Symbol_Table(s1,s2);
    Witness_Sets_io.Add_Embed_Symbols(cd-Symbol_Table.Number);
  end Product_of_Symbol_Tables;

  procedure Read_Witness_Set_for_Diagonal_Homotopy
              ( k : in natural32; n,dim,deg : out natural32;
                fail : out boolean ) is

    lp1,lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    ls1,lsym,solsym : Symbol_Table.Link_to_Array_of_Symbols;
    n1,dim1 : natural32;

  begin
    n := 0;
    fail := true;
    if k = 1 or k = 2 then
      File_Management.Open_Input_File(k);
      get(File_Management.Link_to_Input(k).all,n);
      put("The ambient dimension : "); put(n,1); new_line;
      lp := new Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
      Symbol_Table.Init(n);
      get(File_Management.Link_to_Input(k).all,lp.all);
      dim := Witness_Sets_io.Count_Embed_Symbols(n,"zz");
      put("  dimension of the witness set : "); put(dim,1); new_line;
      Witness_Sets_io.Swap_Symbols_to_End(n,dim,"zz",lp.all);
      Jumpstart_Diagonal_Homotopies.Read_Degree_of_Witness_Set
        (File_Management.Link_to_Input(k).all,deg,n);
      put("  degree of the solution set : "); put(deg,1); new_line;
      lsym := Extrinsic_Diagonal_Homotopies_io.Get_Link_to_Symbols;
      Standard_Solutions_Container.Store_Symbol_Table(k,lsym.all);
      if k = 1 then
        PHCpack_Operations.Store_Target_System(lp.all);
      else 
        PHCpack_Operations.Retrieve_Target_System(lp1);
        n1 := natural32(lp1'last);
        ls1 := Standard_Solutions_Container.Retrieve_Symbol_Table(1);
        dim1 := Witness_Sets_io.Count_Embed_Symbols(ls1.all,"zz");
        put("  n1 = "); put(n1,1); put("  dim1 = "); put(dim1,1); new_line;
        Jumpstart_Diagonal_Homotopies.Match_Symbols
          (n1,n,dim1,dim,lp.all,ls1,lsym,solsym);
        PHCpack_Operations.Store_Start_System(lp.all);
        put_line("Stored start system ...");
        Standard_Solutions_Container.Store_Symbol_Table(0,solsym.all);
        put_line("Stored symbol table ...");
        Product_of_Symbol_Tables(n1,n,dim1,dim,solsym.all);
        put_line("Did product of symbol tables...");
      end if;
      fail := false;
    end if;
  exception
     when others =>
        put_line("Exception raised in Read_Witness_Set_for_Diagonal_Homotopy");
        raise;
  end Read_Witness_Set_for_Diagonal_Homotopy;

  procedure Reset_Witness_Input_File
              ( k : in natural32; deg,dim : out natural32; 
                fail : out boolean ) is

    found : boolean;

  begin
    File_Management.Reset_Input_File(k);
    Scan_and_Skip(File_Management.Link_to_Input(k).all,"SOLUTIONS",found);
    if not found then
      fail := true;
    else
      Standard_Complex_Solutions_io.Read_First
        (File_Management.Link_to_Input(k).all,deg,dim);
      fail := false;
    end if;
  end Reset_Witness_Input_File;

  procedure Read_Start_Solutions is

    file : file_type;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the start solutions...");
    Read_Name_and_Open_File(file);
    get(file,sols);
    Close(file);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Start_Solutions(sols);
    end if;
  end Read_Start_Solutions;

  procedure Read_Start_Solutions ( filename : in string ) is

    file : file_type;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    get(file,sols);
    Close(file);
    PHCpack_Operations.Store_Start_Solutions(sols);
  end Read_Start_Solutions;

  procedure Read_Target_System is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the target system...");
    Standard_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Target_System(p.all);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_Target_System;

  procedure Read_Target_Laurent_System is

    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the target system...");
    Standard_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Target_System(p.all);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_Target_Laurent_System;

  procedure Read_DoblDobl_Target_System is

    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the target system...");
    DoblDobl_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Target_System(p.all);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_DoblDobl_Target_System;

  procedure Read_DoblDobl_Target_Laurent_System is

    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the target system...");
    DoblDobl_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Target_System(p.all);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_DoblDobl_Target_Laurent_System;

  procedure Read_QuadDobl_Target_System is

    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the target system...");
    QuadDobl_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Target_System(p.all);
    if not QuadDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_QuadDobl_Target_System;

  procedure Read_QuadDobl_Target_Laurent_System is

    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the target system...");
    QuadDobl_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Target_System(p.all);
    if not QuadDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_QuadDobl_Target_Laurent_System;

  procedure Read_Multprec_Target_System ( decimals : in natural32 ) is

    p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Multprec_Complex_Solutions.Solution_List;
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(decimals);

  begin
    new_line;
    put_line("Reading the target system...");
    Multprec_Complex_Polynomials_io.Set_Working_Precision(size);
    Multprec_System_and_Solutions_io.get(p,sols);
    PHCpack_Operations.Store_Target_System(p.all);
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      Multprec_Complex_Solutions.Set_Size(sols,size);
      PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_Multprec_Target_System;

  procedure Read_Target_System ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    Standard_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_Target_System;

  procedure Read_Target_Laurent_System ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    Standard_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_Target_Laurent_System;

  procedure Read_DoblDobl_Target_System ( filename : in string ) is

    file : file_type;
    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    DoblDobl_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_DoblDobl_Target_System;

  procedure Read_DoblDobl_Target_Laurent_System ( filename : in string ) is

    file : file_type;
    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    DoblDobl_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_DoblDobl_Target_Laurent_System;

  procedure Read_QuadDobl_Target_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    QuadDobl_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
    if not QuadDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_QuadDobl_Target_System;

  procedure Read_QuadDobl_Target_Laurent_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    QuadDobl_System_and_Solutions_io.get(file,p,sols);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
    if not QuadDobl_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_QuadDobl_Target_Laurent_System;

  procedure Read_Multprec_Target_System
              ( filename : in string; decimals : in natural32 ) is

    file : file_type;
    p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Multprec_Complex_Solutions.Solution_List;
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(decimals);

  begin
    Open(file,in_file,filename);
    Multprec_Complex_Polynomials_io.Set_Working_Precision(size);
    Multprec_System_and_Solutions_io.get(file,p,sols);
    close(file);
    PHCpack_Operations.Store_Target_System(p.all);
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      Multprec_Complex_Solutions.Set_Size(sols,size);
      PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_Multprec_Target_System;

  procedure Read_Target_System_without_Solutions is

    file : file_type;
    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(file);
    get(file,p);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
  end Read_Target_System_without_Solutions;

  procedure Read_Target_System_without_Solutions ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Open(file,in_file,filename);
    get(file,p);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
  end Read_Target_System_without_Solutions;

  procedure Read_DoblDobl_Target_System_without_Solutions is

    file : file_type;
    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(file);
    get(file,p);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
  end Read_DoblDobl_Target_System_without_Solutions;

  procedure Read_QuadDobl_Target_System_without_Solutions is

    file : file_type;
    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(file);
    get(file,p);
    Close(file);
    PHCpack_Operations.Store_Target_System(p.all);
  end Read_QuadDobl_Target_System_without_Solutions;

  procedure Read_Target_Solutions is

    file : file_type;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the target solutions...");
    Read_Name_and_Open_File(file);
    get(file,sols);
    Close(file);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then PHCpack_Operations.Store_Target_Solutions(sols);
    end if;
  end Read_Target_Solutions;

  procedure Read_Target_Solutions ( filename : in string ) is

    file : file_type;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Open(file,in_file,filename);
    get(file,sols);
    Close(file);
    PHCpack_Operations.Store_Target_Solutions(sols);
  end Read_Target_Solutions;

  procedure Write_Start_System is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(p);
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,natural32(p'last),p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put_line(standard_output,"THE START SYSTEM :");
      put(standard_output,natural32(p'last),p.all);
    end if;
  end Write_Start_System;

  procedure Write_Start_Laurent_System is

    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(p);
    if PHCpack_Operations.Is_File_Defined then
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,natural32(p'last),p.all);
    else
      put_line(standard_output,"THE START SYSTEM :");
      put(standard_output,natural32(p'last),p.all);
    end if;
  end Write_Start_Laurent_System;

  procedure Write_DoblDobl_Start_System is

    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(p);
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put_line(standard_output,"THE START SYSTEM :");
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_DoblDobl_Start_System;

  procedure Write_DoblDobl_Start_Laurent_System is

    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(p);
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put_line(standard_output,"THE START SYSTEM :");
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_DoblDobl_Start_Laurent_System;

  procedure Write_QuadDobl_Start_System is

    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(p);
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put_line(standard_output,"THE START SYSTEM :");
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_QuadDobl_Start_System;

  procedure Write_QuadDobl_Start_Laurent_System is

    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(p);
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put_line(standard_output,"THE START SYSTEM :");
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_QuadDobl_Start_Laurent_System;

  procedure Write_Multprec_Start_System is

    p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(p);
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE START SYSTEM :");
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put_line(standard_output,"THE START SYSTEM :");
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_Multprec_Start_System;

  procedure Write_Start_System ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Start_System(p);
    put_line(file,"THE START SYSTEM :");
    put(file,natural32(p'last),p.all);
    Close(file);
  end Write_Start_System;

  procedure Write_Start_Laurent_System ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Start_System(p);
    put_line(file,"THE START SYSTEM :");
    put(file,natural32(p'last),p.all);
    Close(file);
  end Write_Start_Laurent_System;

  procedure Write_DoblDobl_Start_System ( filename : in string ) is

    file : file_type;
    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Start_System(p);
    put_line(file,"THE START SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_DoblDobl_Start_System;

  procedure Write_DoblDobl_Start_Laurent_System ( filename : in string ) is

    file : file_type;
    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Start_System(p);
    put_line(file,"THE START SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_DoblDobl_Start_Laurent_System;

  procedure Write_QuadDobl_Start_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Start_System(p);
    put_line(file,"THE START SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_QuadDobl_Start_System;

  procedure Write_QuadDobl_Start_Laurent_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Start_System(p);
    put_line(file,"THE START SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_QuadDobl_Start_Laurent_System;

  procedure Write_Multprec_Start_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Start_System(p);
    put_line(file,"THE START SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_Multprec_Start_System;

  procedure Write_Start_Solutions is

    sols : Standard_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if not Standard_Complex_Solutions.Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        new_line(PHCpack_Operations.output_file);
        put_line(PHCpack_Operations.output_file,"THE START SOLUTIONS :");
        put(PHCpack_Operations.output_file,
            Standard_Complex_Solutions.Length_Of(sols),
            natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
        text_io.flush(PHCpack_Operations.output_file);
      else
        put_line(standard_output,"THE START SOLUTIONS :");
        put(standard_output,
            Standard_Complex_Solutions.Length_Of(sols),
            natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    end if;
  end Write_Start_Solutions;

  procedure Write_DoblDobl_Start_Solutions is

    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if not DoblDobl_Complex_Solutions.Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        new_line(PHCpack_Operations.output_file);
        put_line(PHCpack_Operations.output_file,"THE START SOLUTIONS :");
        put(PHCpack_Operations.output_file,
            DoblDobl_Complex_Solutions.Length_Of(sols),
            natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
        text_io.flush(PHCpack_Operations.output_file);
      else
        put_line(standard_output,"THE START SOLUTIONS :");
        put(standard_output,
            DoblDobl_Complex_Solutions.Length_Of(sols),
            natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    end if;
  end Write_DoblDobl_Start_Solutions;

  procedure Write_QuadDobl_Start_Solutions is

    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if not QuadDobl_Complex_Solutions.Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        new_line(PHCpack_Operations.output_file);
        put_line(PHCpack_Operations.output_file,"THE START SOLUTIONS :");
        put(PHCpack_Operations.output_file,
            QuadDobl_Complex_Solutions.Length_Of(sols),
            natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
        text_io.flush(PHCpack_Operations.output_file);
      else
        put_line(standard_output,"THE START SOLUTIONS :");
        put(standard_output,
            QuadDobl_Complex_Solutions.Length_Of(sols),
            natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    end if;
  end Write_QuadDobl_Start_Solutions;

  procedure Write_Multprec_Start_Solutions is

    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        new_line(PHCpack_Operations.output_file);
        put_line(PHCpack_Operations.output_file,"THE START SOLUTIONS :");
        put(PHCpack_Operations.output_file,
            Multprec_Complex_Solutions.Length_Of(sols),
            natural32(Multprec_Complex_Solutions.Head_Of(sols).n),sols);
        text_io.flush(PHCpack_Operations.output_file);
      else
        put_line(standard_output,"THE START SOLUTIONS :");
        put(standard_output,
            Multprec_Complex_Solutions.Length_Of(sols),
            natural32(Multprec_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    end if;
  end Write_Multprec_Start_Solutions;

  procedure Write_Start_Solutions ( filename : in string ) is

    file : file_type;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if not Standard_Complex_Solutions.Is_Null(sols) then
      Create(file,out_file,filename);
      put_line(file,"THE START SOLUTIONS :");
      put(file,Standard_Complex_Solutions.Length_Of(sols),
          natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Write_Start_Solutions;

  procedure Write_DoblDobl_Start_Solutions ( filename : in string ) is

    file : file_type;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if not DoblDobl_Complex_Solutions.Is_Null(sols) then
      Create(file,out_file,filename);
      put_line(file,"THE START SOLUTIONS :");
      put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
          natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Write_DoblDobl_Start_Solutions;

  procedure Write_QuadDobl_Start_Solutions ( filename : in string ) is

    file : file_type;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if not QuadDobl_Complex_Solutions.Is_Null(sols) then
      Create(file,out_file,filename);
      put_line(file,"THE START SOLUTIONS :");
      put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
          natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Write_QuadDobl_Start_Solutions;

  procedure Write_Multprec_Start_Solutions ( filename : in string ) is

    file : file_type;
    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      Create(file,out_file,filename);
      put_line(file,"THE START SOLUTIONS :");
      put(file,Multprec_Complex_Solutions.Length_Of(sols),
          natural32(Multprec_Complex_Solutions.Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Write_Multprec_Start_Solutions;

  procedure Write_Target_System is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(p);
    if PHCpack_Operations.Is_File_Defined
     then put(PHCpack_Operations.output_file,natural32(p'last),p.all);
     else put(standard_output,natural32(p'last),p.all);
    end if;
  end Write_Target_System;

  procedure Write_Target_Laurent_System is

    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(p);
    if PHCpack_Operations.Is_File_Defined
     then put(PHCpack_Operations.output_file,natural32(p'last),p.all);
     else put(standard_output,natural32(p'last),p.all);
    end if;
  end Write_Target_Laurent_System;

  procedure Write_DoblDobl_Target_System is

    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(p);
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_DoblDobl_Target_System;

  procedure Write_DoblDobl_Target_Laurent_System is

    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(p);
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_DoblDobl_Target_Laurent_System;

  procedure Write_QuadDobl_Target_System is

    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(p);
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_QuadDobl_Target_System;

  procedure Write_QuadDobl_Target_Laurent_System is

    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(p);
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_QuadDobl_Target_Laurent_System;

  procedure Write_Multprec_Target_System is

    p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(p);
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,p'last,1);
      new_line(PHCpack_Operations.output_file);
      put(PHCpack_Operations.output_file,p.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      put(standard_output,p'last,1);
      new_line(standard_output);
      put(standard_output,p.all);
    end if;
  end Write_Multprec_Target_System;

  procedure Write_Target_System ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Target_System(p);
    put_line(file,"THE TARGET SYSTEM :");
    put(file,natural32(p'last),p.all);
    Close(file);
  end Write_Target_System;

  procedure Write_Target_Laurent_System ( filename : in string ) is

    file : file_type;
    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Target_System(p);
    put_line(file,"THE TARGET SYSTEM :");
    put(file,natural32(p'last),p.all);
    Close(file);
  end Write_Target_Laurent_System;

  procedure Write_DoblDobl_Target_System ( filename : in string ) is

    file : file_type;
    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Target_System(p);
    put_line(file,"THE TARGET SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_DoblDobl_Target_System;

  procedure Write_DoblDobl_Target_Laurent_System ( filename : in string ) is

    file : file_type;
    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Target_System(p);
    put_line(file,"THE TARGET SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_DoblDobl_Target_Laurent_System;

  procedure Write_QuadDobl_Target_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Target_System(p);
    put_line(file,"THE TARGET SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_QuadDobl_Target_System;

  procedure Write_QuadDobl_Target_Laurent_System ( filename : in string ) is

    file : file_type;
    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Target_System(p);
    put_line(file,"THE TARGET SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_QuadDobl_Target_Laurent_System;

  procedure Write_Multprec_Target_System ( filename : in string ) is

    file : file_type;
    p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Create(file,out_file,filename);
    PHCpack_Operations.Retrieve_Target_System(p);
    put_line(file,"THE TARGET SYSTEM :");
    put(file,p'last,1);
    new_line(file);
    put(file,p.all);
    Close(file);
  end Write_Multprec_Target_System;

  procedure Write_Target_Solutions is

    sols : Standard_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if not Standard_Complex_Solutions.Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put_line(PHCpack_Operations.output_file,"THE TARGET SOLUTIONS :");
        put(PHCpack_Operations.output_file,
            Standard_Complex_Solutions.Length_Of(sols),
            natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
        text_io.flush(PHCpack_Operations.output_file);
      else
        put_line(standard_output,"THE TARGET SOLUTIONS :");
        put(standard_output,
            Standard_Complex_Solutions.Length_Of(sols),
            natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    end if;
  end Write_Target_Solutions;

  procedure Write_DoblDobl_Target_Solutions is

    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if not DoblDobl_Complex_Solutions.Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put_line(PHCpack_Operations.output_file,"THE TARGET SOLUTIONS :");
        put(PHCpack_Operations.output_file,
            DoblDobl_Complex_Solutions.Length_Of(sols),
            natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
        text_io.flush(PHCpack_Operations.output_file);
      else
        put_line(standard_output,"THE TARGET SOLUTIONS :");
        put(standard_output,
            DoblDobl_Complex_Solutions.Length_Of(sols),
            natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    end if;
  end Write_DoblDobl_Target_Solutions;

  procedure Write_QuadDobl_Target_Solutions is

    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if not QuadDobl_Complex_Solutions.Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put_line(PHCpack_Operations.output_file,"THE TARGET SOLUTIONS :");
        put(PHCpack_Operations.output_file,
            QuadDobl_Complex_Solutions.Length_Of(sols),
            natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
        text_io.flush(PHCpack_Operations.output_file);
      else
        put_line(standard_output,"THE TARGET SOLUTIONS :");
        put(standard_output,
            QuadDobl_Complex_Solutions.Length_Of(sols),
            natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    end if;
  end Write_QuadDobl_Target_Solutions;

  procedure Write_Multprec_Target_Solutions is

    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put_line(PHCpack_Operations.output_file,"THE TARGET SOLUTIONS :");
        put(PHCpack_Operations.output_file,
            Multprec_Complex_Solutions.Length_Of(sols),
            natural32(Multprec_Complex_Solutions.Head_Of(sols).n),sols);
        text_io.flush(PHCpack_Operations.output_file);
      else
        put_line(standard_output,"THE TARGET SOLUTIONS :");
        put(standard_output,
            Multprec_Complex_Solutions.Length_Of(sols),
            natural32(Multprec_Complex_Solutions.Head_Of(sols).n),sols);
      end if;
    end if;
  end Write_Multprec_Target_Solutions;

  procedure Write_Target_Solutions ( filename : in string ) is

    file : file_type;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if not Standard_Complex_Solutions.Is_Null(sols) then
      Create(file,out_file,filename);
      put_line(file,"THE TARGET SOLUTIONS :");
      put(file,
          Standard_Complex_Solutions.Length_Of(sols),
          natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Write_Target_Solutions;

  procedure Write_DoblDobl_Target_Solutions ( filename : in string ) is

    file : file_type;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if not DoblDobl_Complex_Solutions.Is_Null(sols) then
      Create(file,out_file,filename);
      put_line(file,"THE TARGET SOLUTIONS :");
      put(file,
          DoblDobl_Complex_Solutions.Length_Of(sols),
          natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Write_DoblDobl_Target_Solutions;

  procedure Write_QuadDobl_Target_Solutions ( filename : in string ) is

    file : file_type;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if not QuadDobl_Complex_Solutions.Is_Null(sols) then
      Create(file,out_file,filename);
      put_line(file,"THE TARGET SOLUTIONS :");
      put(file,
          QuadDobl_Complex_Solutions.Length_Of(sols),
          natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Write_QuadDobl_Target_Solutions;

  procedure Write_Multprec_Target_Solutions ( filename : in string ) is

    file : file_type;
    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if not Multprec_Complex_Solutions.Is_Null(sols) then
      Create(file,out_file,filename);
      put_line(file,"THE TARGET SOLUTIONS :");
      put(file,
          Multprec_Complex_Solutions.Length_Of(sols),
          natural32(Multprec_Complex_Solutions.Head_Of(sols).n),sols);
      Close(file);
    end if;
  end Write_Multprec_Target_Solutions;

end PHCpack_Operations_io;
