with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with Degrees_in_Sets_of_Unknowns;
with Standard_Linear_Product_System;
with m_Homogeneous_Bezout_Numbers;
with m_Homogeneous_Start_Systems;        use m_Homogeneous_Start_Systems;
with m_Homogeneous_Permanent_Factors;

procedure ts_permstar is

-- DESCRIPTION :
--   Solving a linear-product start system via the row expansion
--   for the permanent of a degree matrix.

  procedure Solve_m_Homogeneous_Start_System
              ( ind_sols : in Lists_of_Integer_Vectors.List;
                q : out Poly_Sys; qsols : out Solution_List ) is

  -- DESCRIPTION :
  --   Solves the linear-product start system stored in 
  --   Standard_LInear_Product_System, using the list of indices.

    use Lists_of_Integer_Vectors;
    use m_Homogeneous_Permanent_Factors;

    tmp : List := ind_sols;
    lvc : Standard_Integer_Vectors.Link_to_Vector;
    qsols_last : Solution_List;

  begin
    q := Standard_Linear_Product_System.Polynomial_System;
    while not Is_Null(tmp) loop
      lvc := Head_Of(tmp);
      Append(qsols,qsols_last,Solve_Linear_System(lvc.all));
      tmp := Tail_Of(tmp);
    end loop;
  end Solve_m_Homogeneous_Start_System;

  procedure Write_to_File
              ( q : in Poly_Sys; qsols : in Solution_List ) is

  -- DESCRIPTION :
  --   Prompts for a file name and writes the system and its
  --   solutions to file.

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of an output file ...");
    Read_Name_and_Create_File(file);
    Standard_System_and_Solutions_io.put(file,q,qsols);
    Close(file);
  end Write_to_File;

  procedure Partition_and_Degree_Table ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Finds a partition of the set of unknowns of p
  --   and constructs the corresponding degree table.

    n : constant natural32 := natural32(p'last);
    z : Partition(1..n);
    m : natural32;
    b : natural64;
    indsols : Lists_of_Integer_Vectors.List;
    q : Poly_Sys(p'range);
    qsols : Solution_List;

  begin
    m_Homogeneous_Bezout_Numbers.PB(p,b,m,z);
    m_Homogeneous_Bezout_Numbers.Patch(p,z,m,b);
    new_line;
    put("The partition : "); put(z(1..m)); new_line;
    put("with Bezout number : "); put(b,1); put_line(".");
    m_Homogeneous_Permanent_Factors.Permanent_Factors(p,z(1..m),indsols);
    m_Homogeneous_Start_System(p,z(1..m));
    Solve_m_Homogeneous_Start_System(indsols,q,qsols);
    put("Computed "); put(Length_Of(qsols),1); put_line(" solutions.");
    Write_to_File(q,qsols);
  end Partition_and_Degree_Table;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a system and then creates a partition
  --   with corresponding degree table.

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    get(lp);
    Partition_and_Degree_Table(lp.all);
  end Main;

begin
  Main;
end ts_permstar;
