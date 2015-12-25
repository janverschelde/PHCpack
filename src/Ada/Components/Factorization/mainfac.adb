with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Numbers_io;                         use Numbers_io;
with Multprec_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Sampling_Machine;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sample_Lists;              use DoblDobl_Sample_Lists;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sample_Lists;              use QuadDobl_Sample_Lists;
with Standard_Stacked_Sample_Grids;
with DoblDobl_Stacked_Sample_Grids;
with QuadDobl_Stacked_Sample_Grids;
with Multprec_Stacked_Sample_Grids;
with Drivers_to_Grid_Creators;           use Drivers_to_Grid_Creators;
with Standard_Trace_Interpolators;
with Multprec_Trace_Interpolators;
with Homotopy_Membership_Tests;          use Homotopy_Membership_Tests;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Combinatorial_Factorization;        use Combinatorial_Factorization;
with Drivers_to_Factor_Components;       use Drivers_to_Factor_Components;
with Drivers_to_Factor_Polynomials;      use Drivers_to_Factor_Polynomials;
with Driver_for_Common_Factor;
with mainfilt;

procedure mainfac ( infilename,outfilename : in string ) is

  procedure Read_Embedding
               ( filename : in string;
                 lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : out Standard_Complex_Solutions.Solution_List;
                 dim : out natural32 ) is

  -- DESCRIPTION :
  --   Attempts to open the file with name in the string filename
  --   in order to read the embedding.

  -- ON ENTRY :
  --   filename  name of file supplied by user from dispatcher.

  -- ON RETURN :
  --   lp        sliced and embedded polynomial system;
  --   sols      list of generic points on the slices;
  --   dim       dimension of the solution component, or in case of
  --             one polynomial, this is the number of variables.

    file : file_type;

  begin
    Open(file,in_file,filename);
    Standard_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others => new_line;
                   put("Could not open file with name ");
                   put(filename); put_line(".");
                   put("Or incorrect formats of system/solutions in ");
                   put(filename); put_line(".");
                   lp := null; return;
  end Read_Embedding;

  procedure Read_Embedding
               ( filename : in string;
                 lp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : out DoblDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 ) is

  -- DESCRIPTION :
  --   Attempts to open the file with name in the string filename
  --   in order to read the embedding.

  -- ON ENTRY :
  --   filename  name of file supplied by user from dispatcher.

  -- ON RETURN :
  --   lp        sliced and embedded polynomial system;
  --   sols      list of generic points on the slices;
  --   dim       dimension of the solution component, or in case of
  --             one polynomial, this is the number of variables.

    file : file_type;

  begin
    Open(file,in_file,filename);
    DoblDobl_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others => new_line;
                   put("Could not open file with name ");
                   put(filename); put_line(".");
                   put("Or incorrect formats of system/solutions in ");
                   put(filename); put_line(".");
                   lp := null; return;
  end Read_Embedding;

  procedure Read_Embedding
               ( filename : in string;
                 lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 sols : out QuadDobl_Complex_Solutions.Solution_List;
                 dim : out natural32 ) is

  -- DESCRIPTION :
  --   Attempts to open the file with name in the string filename
  --   in order to read the embedding.

  -- ON ENTRY :
  --   filename  name of file supplied by user from dispatcher.

  -- ON RETURN :
  --   lp        sliced and embedded polynomial system;
  --   sols      list of generic points on the slices;
  --   dim       dimension of the solution component, or in case of
  --             one polynomial, this is the number of variables.

    file : file_type;

  begin
    Open(file,in_file,filename);
    QuadDobl_Read_Embedding(file,lp,sols,dim);
    Close(file);
  exception
    when others => new_line;
                   put("Could not open file with name ");
                   put(filename); put_line(".");
                   put("Or incorrect formats of system/solutions in ");
                   put(filename); put_line(".");
                   lp := null; return;
  end Read_Embedding;

  procedure Tune_Member_Tolerances ( restol,homtol : in out double_float ) is

  -- DESCRIPTION :
  --   Displays the current values for the tolerances to decide whether
  --   a residual is small enough for the point to lie on a hypersurface,
  --   and for a test point to belong to a witness set.

    ans : character;

  begin
    loop
      new_line;
      put_line("TOLERANCES to decide whether a test point");
      put("  1. satisfies one polynomial : "); put(restol,3); new_line;
      put("  2. belongs to a witness set : "); put(homtol,3); new_line;
      put("Type 1 or 2 to change a tolerance, 0 to exit : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      put("Give new tolerance for point to ");
      if ans = '1' then
        put("satisfy one polynomial : "); Read_Double_Float(restol);
      else
        put("belong to a witness set : "); Read_Double_Float(homtol);
      end if;
    end loop;
  end Tune_Member_Tolerances;

  procedure Standard_Homotopy_Membership_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using standard double precision arithmetic.

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : double_float := 1.0E-10;
    homtol : double_float := 1.0E-6;

  begin
    Standard_Read_Embedding(lp,genpts,dim);
    new_line;
    put_line("The input format of the test points is the solutions format.");
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Tune_Member_Tolerances(restol,homtol);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
  end Standard_Homotopy_Membership_Test;

  procedure DoblDobl_Homotopy_Membership_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using double double precision arithmetic.

    file : file_type;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : DoblDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : double_float := 1.0E-10;
    homtol : double_float := 1.0E-6;

  begin
    DoblDobl_Read_Embedding(lp,genpts,dim);
    new_line;
    put_line("The input format of the test points is the solutions format.");
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Tune_Member_Tolerances(restol,homtol);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
  end DoblDobl_Homotopy_Membership_Test;

  procedure QuadDobl_Homotopy_Membership_Test is

  -- DESCRIPTION :
  --   Prompts the user for a witness set, a list of test points,
  --   and then applies the homotopy membership test
  --   using quad double precision arithmetic.

    file : file_type;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : QuadDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : double_float := 1.0E-10;
    homtol : double_float := 1.0E-6;

  begin
    QuadDobl_Read_Embedding(lp,genpts,dim);
    new_line;
    put_line("The input format of the test points is the solutions format.");
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Tune_Member_Tolerances(restol,homtol);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
  end QuadDobl_Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test is

    ans : character;

  begin
    new_line;
    put_line("Membership test with homotopy :");
    put_line("  Input : embedded polynomial system with generic points, and");
    put_line("          list of test points.");
    put_line("  Output : decision whether test point lies on component.");
    new_line;
    put_line("MENU to choose the precision : ");
    put_line("  0. standard double precision homotopy continuation; or");
    put_line("  1. double double precision homotopy continuation; or");
    put_line("  2. quad double precision homotopy continuation.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Homotopy_Membership_Test;
      when '1' => DoblDobl_Homotopy_Membership_Test;
      when '2' => QuadDobl_Homotopy_Membership_Test;
      when others => null;
    end case;
  end Homotopy_Membership_Test;

  function Create ( file : file_type;
                    p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists is

  -- DESCRIPTION :
  --   Returns a grid of sample points needed for linear traces,
  --   computed in standard double precision.

    res : Array_of_Standard_Sample_Lists(0..2);
    sli : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    eps,dist : double_float;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    Standard_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    Sampling_Machine.Clear;
    return res;
  end Create;

  function Create ( file : file_type;
                    p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists is

  -- DESCRIPTION :
  --   Returns a grid of sample points needed for linear traces,
  --   computed in double double precision.

    res : Array_of_DoblDobl_Sample_Lists(0..2);
    sli : constant DoblDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant DoblDobl_Sample_List := Create(sols,sli);
    eps,dist : double_double;

  begin
    DoblDobl_Sampling_Machine.Initialize(p);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    DoblDobl_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    DoblDobl_Sampling_Machine.Clear;
    return res;
  end Create;

  function Create ( file : file_type;
                    p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists is

  -- DESCRIPTION :
  --   Returns a grid of sample points needed for linear traces,
  --   computed in double double precision.

    res : Array_of_QuadDobl_Sample_Lists(0..2);
    sli : constant QuadDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant QuadDobl_Sample_List := Create(sols,sli);
    eps,dist : quad_double;

  begin
    QuadDobl_Sampling_Machine.Initialize(p);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    QuadDobl_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    QuadDobl_Sampling_Machine.Clear;
    return res;
  end Create;

  function Append_fk ( name : string; k : natural32 ) return string is

  -- DESCRIPTION :
  --   Returns the name_fk, where the k is the string representation
  --   of the given natural number k.

     strk : constant string := convert(integer32(k));
   
  begin
    return name & "_f" & strk;
  end Append_fk;

  function Select_Witness_Set_for_Factor
             ( sols : in Standard_Complex_Solutions.Solution_List;
               f : in Standard_Natural_Vectors.Vector )
             return Standard_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns those solutions whose index occurs in f.

    use Standard_Complex_Solutions;
    res,res_last : Solution_List;
    ls : Link_to_Solution;

  begin
    for i in f'range loop
      ls := Retrieve(sols,f(i));
      if ls /= null
       then Append(res,res_last,ls.all);
      end if;
    end loop;
    return res;
  end Select_Witness_Set_for_Factor;

  function Select_Witness_Set_for_Factor
             ( sols : in DoblDobl_Complex_Solutions.Solution_List;
               f : in Standard_Natural_Vectors.Vector )
             return DoblDobl_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns those solutions whose index occurs in f.

    use DoblDobl_Complex_Solutions;
    res,res_last : Solution_List;
    ls : Link_to_Solution;

  begin
    for i in f'range loop
      ls := Retrieve(sols,f(i));
      if ls /= null
       then Append(res,res_last,ls.all);
      end if;
    end loop;
    return res;
  end Select_Witness_Set_for_Factor;

  function Select_Witness_Set_for_Factor
             ( sols : in QuadDobl_Complex_Solutions.Solution_List;
               f : in Standard_Natural_Vectors.Vector )
             return QuadDobl_Complex_Solutions.Solution_List is

  -- DESCRIPTION :
  --   Returns those solutions whose index occurs in f.

    use QuadDobl_Complex_Solutions;
    res,res_last : Solution_List;
    ls : Link_to_Solution;

  begin
    for i in f'range loop
      ls := Retrieve(sols,f(i));
      if ls /= null
       then Append(res,res_last,ls.all);
      end if;
    end loop;
    return res;
  end Select_Witness_Set_for_Factor;

  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec ) is

  -- DESCRIPTION  :
  --   Writes a witness set for every irreducible factor of the
  --   embedded system p with corresponding solutions in sols.
  --   The factorization is given in f.

    use Standard_Natural_Vectors;

    k : natural32 := 0;

  begin
    for i in f'range loop
      if f(i) /= null then
        k := k+1;
        declare
          file : file_type;
          filename : constant string := Append_fk(name,k);
          ws : Standard_Complex_Solutions.Solution_List;
        begin
          create(file,out_file,filename);
          put_line(file,p);
          new_line(file);
          put_line(file,"THE SOLUTIONS :");
          ws := Select_Witness_Set_for_Factor(sols,f(i).all);
          put(file,Standard_Complex_Solutions.Length_Of(ws),
              natural32(Standard_Complex_Solutions.Head_Of(ws).n),ws);
          close(file);
        end;
      end if;
    end loop;
  end Write_Witness_Sets_for_Factors;

  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec ) is

  -- DESCRIPTION  :
  --   Writes a witness set for every irreducible factor of the
  --   embedded system p with corresponding solutions in sols.
  --   The factorization is given in f.

    use Standard_Natural_Vectors;

    k : natural32 := 0;

  begin
    for i in f'range loop
      if f(i) /= null then
        k := k+1;
        declare
          file : file_type;
          filename : constant string := Append_fk(name,k);
          ws : DoblDobl_Complex_Solutions.Solution_List;
        begin
          create(file,out_file,filename);
          put_line(file,p);
          new_line(file);
          put_line(file,"THE SOLUTIONS :");
          ws := Select_Witness_Set_for_Factor(sols,f(i).all);
          put(file,DoblDobl_Complex_Solutions.Length_Of(ws),
              natural32(DoblDobl_Complex_Solutions.Head_Of(ws).n),ws);
          close(file);
        end;
      end if;
    end loop;
  end Write_Witness_Sets_for_Factors;

  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec ) is

  -- DESCRIPTION  :
  --   Writes a witness set for every irreducible factor of the
  --   embedded system p with corresponding solutions in sols.
  --   The factorization is given in f.

    use Standard_Natural_Vectors;

    k : natural32 := 0;

  begin
    for i in f'range loop
      if f(i) /= null then
        k := k+1;
        declare
          file : file_type;
          filename : constant string := Append_fk(name,k);
          ws : QuadDobl_Complex_Solutions.Solution_List;
        begin
          create(file,out_file,filename);
          put_line(file,p);
          new_line(file);
          put_line(file,"THE SOLUTIONS :");
          ws := Select_Witness_Set_for_Factor(sols,f(i).all);
          put(file,QuadDobl_Complex_Solutions.Length_Of(ws),
              natural32(QuadDobl_Complex_Solutions.Head_Of(ws).n),ws);
          close(file);
        end;
      end if;
    end loop;
  end Write_Witness_Sets_for_Factors;

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in standard double precision.

    grid : Array_of_Standard_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec := Factor(file,n,grid);
 
  begin
    Deep_Clear(grid);
  end Standard_Enumerate_Decomposition;

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in standard double precision.

    grid : Array_of_Standard_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec := Factor(file,n,grid);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
  end Standard_Enumerate_Decomposition;

  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in double double precision.

    grid : Array_of_DoblDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec := Factor(file,n,grid);
 
  begin
    Deep_Clear(grid);
  end DoblDobl_Enumerate_Decomposition;

  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in double double precision.

    grid : Array_of_DoblDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec := Factor(file,n,grid);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
  end DoblDobl_Enumerate_Decomposition;

  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in double double precision.

    grid : Array_of_QuadDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec := Factor(file,n,grid);
 
  begin
    Deep_Clear(grid);
  end QuadDobl_Enumerate_Decomposition;

  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in double double precision.

    grid : Array_of_QuadDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec := Factor(file,n,grid);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
  end QuadDobl_Enumerate_Decomposition;

  procedure Standard_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using monodromy in standard double precision.

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : Standard_Complex_Poly_Systems.Poly_Sys(ep'range);
 
  begin
    Standard_Complex_Poly_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    Standard_Complex_Poly_Systems.Clear(bp);
  end Standard_Monodromy_Decomposition;

  procedure DoblDobl_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using monodromy in standard double precision.

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : DoblDobl_Complex_Poly_Systems.Poly_Sys(ep'range);
 
  begin
    DoblDobl_Complex_Poly_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    DoblDobl_Complex_Poly_Systems.Clear(bp);
  end DoblDobl_Monodromy_Decomposition;

  procedure QuadDobl_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   Factors a component using monodromy in standard double precision.

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : QuadDobl_Complex_Poly_Systems.Poly_Sys(ep'range);
 
  begin
    QuadDobl_Complex_Poly_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    QuadDobl_Complex_Poly_Systems.Clear(bp);
  end QuadDobl_Monodromy_Decomposition;

  procedure Driver_to_Breakup_Components is

  -- DESCRIPTION :
  --   Reads the embedded polynomial system, display the factorization menu
  --   and either does monodromy loops or enumerates factors.

    infile,outfile : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the file with embedding...");
    declare
      name : constant string := Read_String;
    begin
      Open_Input_File(infile,name);
      Standard_Read_Embedding(infile,lp,sols,dim);
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(outfile);
      new_line;
      put_line("MENU to factor pure dimensional solution component :");
      put_line("  1. use monodromy to group witness points;");
      put_line("  2. enumerate factors and validate by linear traces.");
      put("Type 1 or 2 to select factorization method : ");
      Ask_Alternative(ans,"12");
      if ans = '1' then
        Standard_Monodromy_Decomposition(outfile,name,lp.all,sols,dim);
      else
        new_line;
        put_line("See the output file for results...");
        new_line;
        Standard_Enumerate_Decomposition(outfile,name,lp.all,sols,dim);
      end if;
    end;
  end Driver_to_Breakup_Components;

  procedure DoblDobl_Breakup is

  -- DESCRIPTION :
  --   Prompts the user for a filtered witness set a computes
  --   its numerical irreducible decomposition.
  --   Reads the embedded polynomial system, display the factorization menu
  --   and either does monodromy loops or enumerates factors.

    infile,outfile : file_type;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the file with embedding...");
    declare
      name : constant string := Read_String;
    begin
      Open_Input_File(infile,name);
      DoblDobl_Read_Embedding(infile,lp,sols,dim);
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(outfile);
      new_line;
      put_line("MENU to factor pure dimensional solution component :");
      put_line("  1. use monodromy to group witness points;");
      put_line("  2. enumerate factors and validate by linear traces.");
      put("Type 1 or 2 to select factorization method : ");
      Ask_Alternative(ans,"12");
      if ans = '1' then
        DoblDobl_Monodromy_Decomposition(outfile,name,lp.all,sols,dim);
      else
        new_line;
        put_line("See the output file for results...");
        new_line;
        DoblDobl_Enumerate_Decomposition(outfile,name,lp.all,sols,dim);
      end if;
    end;
  end DoblDobl_Breakup;

  procedure QuadDobl_Breakup is

  -- DESCRIPTION :
  --   Prompts the user for a filtered witness set a computes
  --   its numerical irreducible decomposition.
  --   Reads the embedded polynomial system, display the factorization menu
  --   and either does monodromy loops or enumerates factors.

    infile,outfile : file_type;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the file with embedding...");
    declare
      name : constant string := Read_String;
    begin
      Open_Input_File(infile,name);
      QuadDobl_Read_Embedding(infile,lp,sols,dim);
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(outfile);
      new_line;
      put_line("MENU to factor pure dimensional solution component :");
      put_line("  1. use monodromy to group witness points;");
      put_line("  2. enumerate factors and validate by linear traces.");
      put("Type 1 or 2 to select factorization method : ");
      Ask_Alternative(ans,"12");
      if ans = '1' then
        QuadDobl_Monodromy_Decomposition(outfile,name,lp.all,sols,dim);
      else
        new_line;
        put_line("See the output file for results...");
        new_line;
        QuadDobl_Enumerate_Decomposition(outfile,name,lp.all,sols,dim);
      end if;
    end;
  end QuadDobl_Breakup;

  procedure Trace_Form_Interpolation is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("MENU to certify monodromy breakup with by interpolation :");
    put_line
      ("  1. on given decomposition : use bootstrapping Newton to certify;");
    put_line
      ("  2.                        : use full trace form to certify;");
    put_line
      ("  3.                        : use Newton identities on trace form;");
    put_line
      ("  4.                        : use linear trace only to certify.");
    put("Type 1, 2, 3, or 4 to make your choice : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Call_Newton_Interpolate(file,lp.all,sols,dim);
      when '2' => Call_Trace_Interpolate(file,lp.all,sols,dim);
      when '3' => Call_Power_Trace_Interpolate(file,lp.all,sols,dim);
      when '4' => Call_Linear_Trace_Interpolate(file,lp.all,sols,dim);
      when others => null;
    end case;
  end Trace_Form_Interpolation;

  procedure Incremental_Interpolation is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("MENU to decompose with incremental use of interpolation :");
    put_line("  1. massive interpolate with standard arithmetic;");
    put_line("  2. incremental interpolate with standard arithmetic;");
    put_line("  3.  + determine span with standard arithmetic;");
    put_line("  4.  + use central projections;");
    put_line("  5. massive interpolate with multi-precision arithmetic;");
    put_line("  6. incremental interpolate with multi-precision arithmetic;");
    put_line("  7.  + determine span with multi-precision arithmetic;");
    put_line("  8.  + use central projections;");
    Ask_Alternative(ans,"12345678");
    case ans is
      when '1' => Call_Standard_Interpolate(file,lp.all,sols,dim,0);
      when '2' => Call_Standard_Interpolate(file,lp.all,sols,dim,1);
      when '3' => Call_Standard_Interpolate(file,lp.all,sols,dim,2);
      when '4' => Call_Standard_Interpolate(file,lp.all,sols,dim,3);
      when '5' => Call_Multprec_Interpolate(file,lp.all,sols,dim,0);
      when '6' => Call_Multprec_Interpolate(file,lp.all,sols,dim,1);
      when '7' => Call_Multprec_Interpolate(file,lp.all,sols,dim,2);
      when '8' => Call_Multprec_Interpolate(file,lp.all,sols,dim,3);
      when others => null;
    end case;
  end Incremental_Interpolation;

  procedure Standard_Eliminate
                ( file : in file_type;
                  p : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  dim : in integer32 ) is

  -- DESCRIPTION :
  --   Given the embedded system with its solutions, this procedure
  --   initializes the sampler and calls the interpolation routines.
  --   All calculations are done with standard floating-point arithmetic.

    use Standard_Stacked_Sample_Grids,Standard_Trace_Interpolators;

    sli : constant Standard_Complex_VecVecs.VecVec := Slices(p,natural32(dim));
    sps : constant Standard_Sample_List := Create(sols,sli);
    deg : constant integer32 := integer32(Length_Of(sps));
    ip : Standard_Complex_Polynomials.Poly;
 
  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if dim = 1 then
      declare
        grid : Array_of_Standard_Sample_Lists(0..deg);
        eps,dst : double_float;
        t : Trace_Interpolator1;  
      begin
        Standard_Rectangular_Grid_Creator
          (file,sps,natural32(deg),grid,eps,dst);
        put(file,"Maximal error of the samples on the grid : ");
        put(file,eps,3); new_line(file);
        put(file,"Minimal distance between samples in one list in grid :");
        put(file,dst,3); new_line(file);
        t := Create(file,grid);
        put(file,"Maximal residual of evaluation at the grid : ");
        put(file,Maximal_Error(t,grid),3); new_line(file);    
        ip := Expand(t);
      end;
    else
      declare
        grid : Stacked_Sample_Grid(dim,deg);
        t : Trace_Interpolator;
      begin
        Standard_Stacked_Grid_Creator(file,sps,true,grid); 
        t := Create(file,grid,deg);
        put(file,"Maximal residual of evaluation at the grid : ");
        put(file,Maximal_Error(t,grid),3); new_line(file);    
        ip := Expand(t);
      end;
    end if;
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,ip);
  end Standard_Eliminate;

  procedure Multprec_Eliminate
                ( file : in file_type;
                  ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                  mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  dim : in integer32; size : in natural32 ) is       

  -- DESCRIPTION :
  --   Given the embedded system with its solutions, this procedure
  --   initializes the sampler and calls the interpolation routines.
  --   Sample refinement and interpolation is done with multi-precision
  --   arithmetic.

    use Multprec_Stacked_Sample_Grids,Multprec_Trace_Interpolators;

    sli : constant Standard_Complex_VecVecs.VecVec := Slices(ep,natural32(dim));
    sps : Standard_Sample_List := Create(sols,sli);
    deg : constant integer32 := integer32(Length_Of(sps));
    ip : Multprec_Complex_Polynomials.Poly;

  begin
    Sampling_Machine.Initialize(ep,mp,dim,size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    new_line;
    put_line("See the output file for results...");
    new_line;                                                
    if dim = 1 then
      declare
        grid : Array_of_Multprec_Sample_Lists(0..deg);
        eps,dst : double_float;
        t : Trace_Interpolator1;
      begin
        Multprec_Rectangular_Grid_Creator
          (file,sps,natural32(deg),size,grid,eps,dst);
        put(file,"Maximal error of the samples on the grid : ");
        put(file,eps,3); new_line(file);
        put(file,"Minimal distance between samples in one list in grid :");
        put(file,dst,3); new_line(file);
        t := Create(file,grid);
        ip := Expand(t);
      end;
    else
      declare
        grid : Stacked_Sample_Grid(dim,deg);
        t : Trace_Interpolator;
      begin
        Multprec_Stacked_Grid_Creator(file,sps,true,size,grid);
        t := Create(file,grid,deg);
        ip := Expand(t);
      end;
    end if;
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,ip);
  end Multprec_Eliminate;

  procedure Eliminate_Variables is

  -- DESCRIPTION :
  --   Given an embedding with slices parallel to a given subspace,
  --   with interpolation we eliminate variables not belonging to this
  --   subspace from the system.

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim,deci,size : natural32 := 0;

  begin
   -- new_line;
   -- put_line
   --   ("Numerical Elimination by Sampling, Projecting, and Interpolating.");
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    Determine_Order(lp.all,sols);
    new_line;
    put("Give the number of decimal places (<= 16 is standard) : ");
    get(deci);
    new_line;
    if deci > 16 then
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(lp.all,mp,size,dim);
      Multprec_Eliminate(file,lp.all,mp.all,sols,integer32(dim),size);
    else
      Standard_Eliminate(file,lp.all,sols,integer32(dim));
    end if;
  end Eliminate_Variables;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to filter junk, factor components, and eliminate :");
    put_line("  0. filter solution lists subject to criteria;");
    put_line("  1. filter junk with homotopy membership test;");
    put_line("  2. breakup a filtered witness point set into irreducibles;");
    put_line("  3. given partition of breakup, compute trace form of filter;");
    put_line("  4. perform tasks 1, 2, and 3 by incremental interpolation;");
    put_line("  5. eliminate variables by interpolation via special slices;");
    put_line("  6. factor a complex polynomial in several variables;");
    put_line("  7. detect a common factor of two Laurent polynomials;");
    put_line("  8. filtered witness set breakup in double double precision;");
    put_line("  9. filtered witness set breakup in quad double precision.");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, or 9 to select a task : ");
    Ask_Alternative(ans,"0123456789");
    case ans is
      when '0' => mainfilt(infilename,outfilename);
      when '1' => Homotopy_Membership_Test;
      when '2' => Driver_to_Breakup_Components;
      when '3' => Trace_Form_Interpolation;
      when '4' => Incremental_Interpolation;
      when '5' => Eliminate_Variables;
      when '6' =>
        if infilename = ""
         then Driver_to_Factor_Polynomial;
         else Driver_to_Factor_Polynomial(infilename);
        end if;
      when '7' => Driver_for_Common_Factor;
      when '8' => DoblDobl_Breakup;
      when '9' => QuadDobl_Breakup;
      when others => null;
    end case;
  end Main;

begin
  Main;
end mainfac;
