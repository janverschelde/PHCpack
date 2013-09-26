with text_io,integer_io;                 use text_io,integer_io;
with File_Scanning;                      use File_Scanning;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with C_Integer_Arrays,C_Double_Arrays;   use C_Integer_Arrays,C_Double_Arrays;
with C_to_Ada_Arrays;                    use C_to_Ada_Arrays;
with Coefficient_Support_Poly_Systems;   use Coefficient_Support_Poly_Systems;
with Coefficient_Solution_Vectors;       use Coefficient_Solution_Vectors;

-- NOTICE :
--   The following packages are included for use in path_tracks, to make
--   sure that the linker will get these in the .ali files.

with C_Integer_Arrays_io;                use C_Integer_Arrays_io;
with C_Double_Arrays_io;                 use C_Double_Arrays_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Homotopy;
with Increment_and_Fix_Continuation;     use Increment_and_Fix_Continuation;
with Standard_Root_Refiners;             use Standard_Root_Refiners;

procedure ts_path_tracker is

-- DESCRIPTION :
--   Test on the C interface to the path tracking routines.

  function path_tracker ( n,m : integer;
                          p_ms : C_Integer_Array;
                          p_ns : integer; p_s : C_Integer_Array;
                          p_nc : integer; p_c : C_Double_Array;
                          q_ms : C_Integer_Array;
                          q_ns : integer; q_s : C_Integer_Array;
                          q_nc : integer; q_c : C_Double_Array;
                          nbmul : integer; mul : C_Integer_Array;
                          nbcfs : integer; cfs : C_Double_Array ) 
                        return integer;
  pragma Import(C,path_tracker,"path_tracker");

  -- DESCRIPTION :
  --   The function path_tracker calls the Ada path tracking routines.

  procedure Call_Path_Tracker 
              ( p,q : in Poly_Sys; sols : in Solution_List ) is

    n : constant natural := p'length;
    m : constant natural := Number_of_Unknowns(p(p'first));
    p_moncnt : constant C_Integer_Array := Monomial_Count(p);
    p_monsum : constant natural := Sum(p_moncnt);
    p_sup : constant C_Integer_Array := Support(n,p_monsum,p_moncnt,p);
    p_cff : constant C_Double_Array := Coefficients(p_monsum,p_moncnt,p);
    q_moncnt : constant C_Integer_Array := Monomial_Count(q);
    q_monsum : constant natural := Sum(q_moncnt);
    q_sup : constant C_Integer_Array := Support(n,q_monsum,q_moncnt,q);
    q_cff : constant C_Double_Array := Coefficients(q_monsum,q_moncnt,q);
    mul : constant C_Integer_Array := Multiplicities(sols);
    cfs : constant C_Double_Array := Coefficients(sols);
   -- mul : constant Standard_Integer_Vectors.Vector := Multiplicities(sols);
   -- cfs : constant Standard_Floating_Vectors.Vector := Coefficients(sols);
   -- c_mul : constant C_Integer_Array := Convert(mul);
   -- c_cfs : constant C_Double_Array := Convert(cfs);
    return_of_call : integer;

  begin
    return_of_call
      := path_tracker(n,m,p_moncnt,p_sup'length,p_sup,p_cff'length,p_cff,
                          q_moncnt,q_sup'length,q_sup,q_cff'length,q_cff,
                         -- c_mul'length,c_mul,c_cfs'length,c_cfs);
                          mul'length,mul,cfs'length,cfs);
    if return_of_call = 0
     then put_line("Call to C terminated normally.");
     else put("Call to C terminated abnormally with exit code ");
          put(return_of_call,1); new_line;
    end if;
  end Call_Path_Tracker;

  procedure Read_Start_System
               ( q : out Poly_Sys; sols : out Solution_List ) is

  -- DESCRIPTION :
  --   Reads the name of a file for start system q with solutions sols.

    file : file_type;
    n : natural;
    found : boolean;

  begin
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(file);
    get(file,n,q);
    Scan_and_Skip(file,"SOLUTIONS",found);
    if found
     then get(file,sols);
     else new_line;
          Read(sols);
    end if;
  end Read_Start_System;

  procedure Main is

    lp : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Testing the path tracking routines.");
    new_line;
    put_line("Reading the target system...");
    get(lp);
    declare
      q : Poly_Sys(lp'range);
      ans : character;
      newt : Complex_Number;
    begin
      Read_Start_System(q,sols);
      if Head_Of(sols).t = Create(1.0)
       then put("The value for t is one.  Do you want to change ? (y/n) ");
            Ask_Yes_or_No(ans);
            if ans = 'y'
             then put("Give new (complex) value for t : "); get(newt);
                  Set_Continuation_Parameter(sols,newt);
            end if;
      end if;
      Call_Path_Tracker(lp.all,q,sols);
    end;
  end Main;

begin
  Main;
end ts_path_tracker;

