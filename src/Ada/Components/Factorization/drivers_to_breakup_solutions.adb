with Communications_with_User;          use Communications_with_User;
with String_Splitters;                  use String_Splitters;
with Characters_and_Numbers;            use Characters_and_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Witness_Sets,Witness_Sets_io;
with Sampling_Machine; 
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Machine; 
with DoblDobl_Sampling_Laurent_Machine; 
with QuadDobl_Sampling_Machine; 
with QuadDobl_Sampling_Laurent_Machine; 
with Rectangular_Sample_Grids;
with DoblDobl_Rectangular_Sample_Grids;
with QuadDobl_Rectangular_Sample_Grids;
with Make_Sample_Grids;                 use Make_Sample_Grids;
with Combinatorial_Factorization;
with Drivers_to_Factor_Components;
with Greeting_Banners;
with Write_Seed_Number;

package body Drivers_to_Breakup_Solutions is

  function Create ( file : file_type;
                    p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists is

    res : Array_of_Standard_Sample_Lists(0..2);
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(p,dim);
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
                    p : Standard_Complex_Laur_Systems.Laur_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists is

    res : Array_of_Standard_Sample_Lists(0..2);
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    eps,dist : double_float;

  begin
    Sampling_Laurent_Machine.Initialize(p);
    Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    Standard_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    Sampling_Laurent_Machine.Clear;
    return res;
  end Create;

  function Create ( file : file_type;
                    p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(0..2);
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(p,dim);
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
                    p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(0..2);
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(p,dim);
    sps : constant DoblDobl_Sample_List := Create(sols,sli);
    eps,dist : double_double;

  begin
    DoblDobl_Sampling_Laurent_Machine.Initialize(p);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    DoblDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    DoblDobl_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    DoblDobl_Sampling_Laurent_Machine.Clear;
    return res;
  end Create;

  function Create ( file : file_type;
                    p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..2);
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(p,dim);
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

  function Create ( file : file_type;
                    p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..2);
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(p,dim);
    sps : constant QuadDobl_Sample_List := Create(sols,sli);
    eps,dist : quad_double;

  begin
    QuadDobl_Sampling_Laurent_Machine.Initialize(p);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    QuadDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    QuadDobl_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    QuadDobl_Sampling_Laurent_Machine.Clear;
    return res;
  end Create;

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_Standard_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end Standard_Enumerate_Decomposition;

  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_DoblDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end DoblDobl_Enumerate_Decomposition;

  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_QuadDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end QuadDobl_Enumerate_Decomposition;

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_Standard_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end Standard_Enumerate_Decomposition;

  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_DoblDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end DoblDobl_Enumerate_Decomposition;

  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_QuadDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end QuadDobl_Enumerate_Decomposition;

  function Append_fk ( name : string; k : natural32 ) return string is

    strk : constant string := convert(integer32(k));
   
  begin
    return name & "_f" & strk;
  end Append_fk;

  function Select_Witness_Set_for_Factor
             ( sols : in Standard_Complex_Solutions.Solution_List;
               f : in Standard_Natural_Vectors.Vector )
             return Standard_Complex_Solutions.Solution_List is

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

  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec ) is

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
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec ) is

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
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec ) is

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
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_Standard_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end Standard_Enumerate_Decomposition;

  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_DoblDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end DoblDobl_Enumerate_Decomposition;

  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_QuadDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end QuadDobl_Enumerate_Decomposition;

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_Standard_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end Standard_Enumerate_Decomposition;

  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_DoblDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end DoblDobl_Enumerate_Decomposition;

  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    grid : Array_of_QuadDobl_Sample_Lists(0..2) := Create(file,ep,sols,dim);
    n : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    f : constant Standard_Natural_VecVecs.VecVec
      := Combinatorial_Factorization.Factor(file,n,grid);
    lf : Standard_Natural_VecVecs.Link_to_VecVec
       := new Standard_Natural_VecVecs.VecVec'(f);
 
  begin
    Deep_Clear(grid);
    Write_Witness_Sets_for_Factors(name,ep,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(lf);
  end QuadDobl_Enumerate_Decomposition;

  procedure Standard_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : Standard_Complex_Poly_Systems.Poly_Sys(ep'range);

    use Drivers_to_Factor_Components;
 
  begin
    Standard_Complex_Poly_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    Standard_Complex_Poly_Systems.Clear(bp);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end Standard_Monodromy_Decomposition;

  procedure DoblDobl_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : DoblDobl_Complex_Poly_Systems.Poly_Sys(ep'range);

    use Drivers_to_Factor_Components;
 
  begin
    DoblDobl_Complex_Poly_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    DoblDobl_Complex_Poly_Systems.Clear(bp);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end DoblDobl_Monodromy_Decomposition;

  procedure QuadDobl_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : QuadDobl_Complex_Poly_Systems.Poly_Sys(ep'range);

    use Drivers_to_Factor_Components;
 
  begin
    QuadDobl_Complex_Poly_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    QuadDobl_Complex_Poly_Systems.Clear(bp);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end QuadDobl_Monodromy_Decomposition;

  procedure Standard_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : Standard_Complex_Laur_Systems.Laur_Sys(ep'range);

    use Drivers_to_Factor_Components;
 
  begin
    Standard_Complex_Laur_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    Standard_Complex_Laur_Systems.Clear(bp);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end Standard_Monodromy_Decomposition;

  procedure DoblDobl_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : DoblDobl_Complex_Laur_Systems.Laur_Sys(ep'range);

    use Drivers_to_Factor_Components;
 
  begin
    DoblDobl_Complex_Laur_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    DoblDobl_Complex_Laur_Systems.Clear(bp);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end DoblDobl_Monodromy_Decomposition;

  procedure QuadDobl_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;
    bp : QuadDobl_Complex_Laur_Systems.Laur_Sys(ep'range);

    use Drivers_to_Factor_Components;
 
  begin
    QuadDobl_Complex_Laur_Systems.Copy(ep,bp);
    Call_Monodromy_Breakup(file,ep,sols,dim,f);
    Write_Witness_Sets_for_Factors(name,bp,sols,dim,f.all);
    QuadDobl_Complex_Laur_Systems.Clear(bp);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end QuadDobl_Monodromy_Decomposition;

  function Prompt_for_Method return character is

    ans : character;

  begin
    new_line;
    put_line("MENU to factor pure dimensional solution component :");
    put_line("  1. use monodromy to group witness points;");
    put_line("  2. enumerate factors and validate by linear traces.");
    put("Type 1 or 2 to select factorization method : ");
    Ask_Alternative(ans,"12");
    return ans;
  end Prompt_for_Method;

  procedure Standard_Breakup
              ( file : in file_type; name : in string;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    ans : constant character := Prompt_for_Method;

  begin
    if ans = '1' then
      Standard_Monodromy_Decomposition(file,name,p,sols,dim);
    else
      new_line;
      put_line("See the output file for results...");
      new_line;
      Standard_Enumerate_Decomposition(file,name,p,sols,dim);
    end if;
  end Standard_Breakup;

  procedure DoblDobl_Breakup
              ( file : in file_type; name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    ans : constant character := Prompt_for_Method;

  begin
    if ans = '1' then
      DoblDobl_Monodromy_Decomposition(file,name,p,sols,dim);
    else
      new_line;
      put_line("See the output file for results...");
      new_line;
      DoblDobl_Enumerate_Decomposition(file,name,p,sols,dim);
    end if;
  end DoblDobl_Breakup;

  procedure QuadDobl_Breakup
              ( file : in file_type; name : in string;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    ans : constant character := Prompt_for_Method;

  begin
    if ans = '1' then
      QuadDobl_Monodromy_Decomposition(file,name,p,sols,dim);
    else
      new_line;
      put_line("See the output file for results...");
      new_line;
      QuadDobl_Enumerate_Decomposition(file,name,p,sols,dim);
    end if;
  end QuadDobl_Breakup;

  procedure Standard_Breakup
              ( file : in file_type; name : in string;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    ans : constant character := Prompt_for_Method;

  begin
    if ans = '1' then
      Standard_Monodromy_Decomposition(file,name,p,sols,dim);
    else
      new_line;
      put_line("See the output file for results...");
      new_line;
      Standard_Enumerate_Decomposition(file,name,p,sols,dim);
    end if;
  end Standard_Breakup;

  procedure DoblDobl_Breakup
              ( file : in file_type; name : in string;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    ans : constant character := Prompt_for_Method;

  begin
    if ans = '1' then
      DoblDobl_Monodromy_Decomposition(file,name,p,sols,dim);
    else
      new_line;
      put_line("See the output file for results...");
      new_line;
      DoblDobl_Enumerate_Decomposition(file,name,p,sols,dim);
    end if;
  end DoblDobl_Breakup;

  procedure QuadDobl_Breakup
              ( file : in file_type; name : in string;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    ans : constant character := Prompt_for_Method;

  begin
    if ans = '1' then
      QuadDobl_Monodromy_Decomposition(file,name,p,sols,dim);
    else
      new_line;
      put_line("See the output file for results...");
      new_line;
      QuadDobl_Enumerate_Decomposition(file,name,p,sols,dim);
    end if;
  end QuadDobl_Breakup;

  procedure Standard_Breakup ( infilename,outfilename : in string ) is

    use Standard_Laur_Poly_Convertors;

    infile,outfile : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;

  begin
    if infilename /= "" then
      Witness_Sets_io.Standard_Read_Embedding(infilename,lq,sols,dim);
      Create_Output_File(outfile,outfilename);
      if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
        Standard_Breakup(outfile,infilename,lq.all,sols,dim);
      else
        lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                    (Positive_Laurent_Polynomial_System(lq.all));
        Standard_Breakup(outfile,infilename,lp.all,sols,dim);
      end if;
    else
      new_line;
      put_line("Reading the name of the file with embedding...");
      declare
        name : constant string := Read_String;
      begin
        Open_Input_File(infile,name);
        Witness_Sets_io.Standard_Read_Embedding(infile,lq,sols,dim);
        Create_Output_File(outfile,outfilename);
        if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
          Standard_Breakup(outfile,infilename,lq.all,sols,dim);
        else
          lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                      (Positive_Laurent_Polynomial_System(lq.all));
          Standard_Breakup(outfile,name,lp.all,sols,dim);
        end if;
      end;
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Standard_Breakup;

  procedure DoblDobl_Breakup ( infilename,outfilename : in string ) is

    use DoblDobl_Laur_Poly_Convertors;

    infile,outfile : file_type;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    dim : natural32;

  begin
    if infilename /= "" then
      Witness_Sets_io.DoblDobl_Read_Embedding(infilename,lq,sols,dim);
      Create_Output_File(outfile,outfilename);
      if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
        DoblDobl_Breakup(outfile,infilename,lq.all,sols,dim);
      else
        lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                    (Positive_Laurent_Polynomial_System(lq.all));
        DoblDobl_Breakup(outfile,infilename,lp.all,sols,dim);
      end if;
    else
      new_line;
      put_line("Reading the name of the file with embedding...");
      declare
        name : constant string := Read_String;
      begin
        Open_Input_File(infile,name);
        Witness_Sets_io.DoblDobl_Read_Embedding(infile,lq,sols,dim);
        Create_Output_File(outfile,outfilename);
        if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
          DoblDobl_Breakup(outfile,name,lq.all,sols,dim);
        else
          lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                      (Positive_Laurent_Polynomial_System(lq.all));
          DoblDobl_Breakup(outfile,name,lp.all,sols,dim);
        end if;
      end;
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end DoblDobl_Breakup;

  procedure QuadDobl_Breakup ( infilename,outfilename : in string ) is

    use QuadDobl_Laur_Poly_Convertors;

    infile,outfile : file_type;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    dim : natural32;

  begin
    if infilename /= "" then
      Witness_Sets_io.QuadDobl_Read_Embedding(infilename,lq,sols,dim);
      Create_Output_File(outfile,outfilename);
      if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
        QuadDobl_Breakup(outfile,infilename,lq.all,sols,dim);
      else
        lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                    (Positive_Laurent_Polynomial_System(lq.all));
        QuadDobl_Breakup(outfile,infilename,lp.all,sols,dim);
      end if;
    else
      new_line;
      put_line("Reading the name of the file with the embedding...");
      declare
        name : constant string := Read_String;
      begin
        Open_Input_File(infile,name);
        Witness_Sets_io.QuadDobl_Read_Embedding(infile,lq,sols,dim);
        Create_Output_File(outfile,outfilename);
        if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
          QuadDobl_Breakup(outfile,name,lq.all,sols,dim);
        else
          lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                      (Positive_Laurent_Polynomial_System(lq.all));
          QuadDobl_Breakup(outfile,name,lp.all,sols,dim);
        end if;
      end;
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end QuadDobl_Breakup;

end Drivers_to_Breakup_Solutions;
