with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Sets_of_Unknowns;                  use Sets_of_Unknowns;
with Sets_of_Unknowns_io;               use Sets_of_Unknowns_io;
with Partitions_of_Sets_of_Unknowns;    use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io; use Partitions_of_Sets_of_Unknowns_io;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with P_Intrinsic_Diagonal_Homotopies;   use P_Intrinsic_Diagonal_Homotopies;
with P_Intrinsic_Diagonal_Solvers;      use P_Intrinsic_Diagonal_Solvers;

package body Drivers_to_Intrinsic_Solvers is

-- SOME AUXILIARIES :

  procedure Read_Polynomials ( n : out natural32; p,q : out Poly ) is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables n, and then two
  --   polynomials in n variables.

  begin
    n := 0;
    new_line;
    put("Give the number of variables : "); get(n);
    put_line("Give a polynomial p (terminated with semicolon) : ");
    Symbol_Table.Init(n);
    get(p);
    put("Your polynomial p : "); put(p); new_line;
    put_line("Give a polynomial q (terminated with semicolon) : ");
    get(q);
    put("Your polynomial q : "); put(q); new_line;
  end Read_Polynomials;

  procedure Write_Polynomials ( file : in file_type; p,q : in Poly ) is

  -- DESCRIPTION :
  --   Writes the "See output file for results..." message to screen
  --   and the polynomials p and q to file.

  begin
    new_line;
    put_line("See the output file for results...");
    new_line;
    put_line(file,"Intersecting the polynomial ");
    put(file,p); new_line(file);
    put_line(file,"with the polynomial ");
    put(file,q); new_line(file);
  end Write_Polynomials;

  function Read_Partition ( n : natural32 ) return Partition is

  -- DESCRIPTION :
  --   Prompts the user for the number of sets and the sets
  --   in the partition of a set of n variables.

  -- REQUIRED :
  --   The symbol table has been initialized.

    m : natural32 := 0;
    variables : Set;

  begin
    variables := Create(n);
    for i in 1..n loop
      Add(variables,i);
    end loop;
    put("The set of all variables : "); put(variables); new_line; 
    put("Give the number of sets in the partition : "); get(m);
    declare
      z : Partition(1..m);
    begin
      Create(z,n);
      put_line("Reading the partition of the set of variables...");
      put_line("Begin and end each set with a curly brace.");
      for i in 1..m loop
        put("Give set "); put(i,1); put(" : "); get(z(i));
      end loop;
      put("The partition : "); put(z); new_line;
      return z;
    end;
  end Read_Partition; 

-- DRIVERS TO SOLVERS FOR TWO POLYNOMIALS :

  procedure Total_Degree_Hypersurface_Intersection
               ( n : in natural32; p,q : in Poly ) is

    n2 : constant natural32 := 2*n;
    b2 : Vector(1..integer32(n2));
    w2 : VecVec(1..2);
    sols : Array_of_Solution_Lists(1..2);
    fail : boolean;

  begin
    Witness_Points(n,p,q,b2,w2,sols,fail);
    if not fail then
      if not Is_Null(sols(1))
       then Combine_Solutions(Standard_Output,n,sols(1),b2,w2(1..1));
      end if;
      if not Is_Null(sols(2))
       then Combine_Solutions(Standard_Output,n,sols(2),b2,w2);
      end if;
    end if;
  end Total_Degree_Hypersurface_Intersection;

  procedure Total_Degree_Hypersurface_Intersection
               ( file : in file_type; n : in natural32; p,q : in Poly ) is

    n2 : constant natural32 := 2*n;
    b2 : Vector(1..integer32(n2));
    w2 : VecVec(1..2);
    sols : Array_of_Solution_Lists(1..2);
    fail : boolean;

  begin
    Witness_Points(file,n,p,q,b2,w2,sols,fail);
    if not fail then
      if not Is_Null(sols(1))
       then Combine_Solutions(file,n,sols(1),b2,w2(1..1));
      end if;
      if not Is_Null(sols(2))
       then Combine_Solutions(file,n,sols(2),b2,w2);
      end if;
    end if;
  end Total_Degree_Hypersurface_Intersection;

  procedure Multihomogeneous_Hypersurface_Intersection
               ( n : in natural32; p,q : in Poly; z : in Partition ) is

    n2 : constant natural32 := 2*n;
    m : constant natural32 := natural32(z'length);
    m2 : constant natural32 := m*m;
    sols : Array_of_Solution_Lists(1..integer32(m2));
    b2 : Vector(1..integer32(n2));
    w2 : Array_of_VecVecs(1..integer32(m2));

  begin
    Witness_Points(n,p,q,z,b2,w2,sols);
    Combine_Solutions(Standard_Output,n,sols,b2,w2);
  end Multihomogeneous_Hypersurface_Intersection;

  procedure Multihomogeneous_Hypersurface_Intersection
               ( file : in file_type;
                 n : in natural32; p,q : in Poly; z : in Partition ) is

    n2 : constant natural32 := 2*n;
    m : constant natural32 := natural32(z'length);
    m2 : constant natural32 := m*m;
    b2 : Vector(1..integer32(n2));
    w2 : Array_of_VecVecs(1..integer32(m2));
    sols : Array_of_Solution_Lists(1..integer32(m2));

  begin
    Witness_Points(file,n,p,q,z,b2,w2,sols);
    Combine_Solutions(file,n,sols,b2,w2);
  end Multihomogeneous_Hypersurface_Intersection;

-- MAIN INTERACTIVE DRIVERS ROUTINES :

  procedure Intersection_of_Two_Hypersurfaces is

    ans : character;
    file : file_type;
    file_output : boolean;
    n : natural32 := 0;
    p,q : Poly;

  begin
    new_line;
    put("Do you wish intermediate output on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      file_output := true;
    else
      file_output := false;
    end if;
    Read_Polynomials(n,p,q);
    put("Is there a multi-homogeneous structure to respect ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        z : constant Partition := Read_Partition(n);
      begin
        if file_output then
          Write_Polynomials(file,p,q);
          put(file,"The multi-homogeneous structure is ");
          put(file,z); new_line(file);
          Multihomogeneous_Hypersurface_Intersection(file,n,p,q,z);
        else
          Multihomogeneous_Hypersurface_Intersection(Standard_Output,n,p,q,z);
        end if;
      end;
    else
      if file_output then
        Write_Polynomials(file,p,q);
        Total_Degree_Hypersurface_Intersection(file,n,p,q);
      else
        Total_Degree_Hypersurface_Intersection(n,p,q);
      end if;
    end if;
  end Intersection_of_Two_Hypersurfaces;

  procedure Driver_to_Hypersurface_Solvers is

    timer : Timing_Widget;
    file : file_type;
    lp : Link_to_Poly_Sys;
    n : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Reading a polynomial system...");
    get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    n := Number_of_Unknowns(lp(lp'first));
    put(file,lp'last,n,lp.all);
    Driver_for_Continuation_Parameters(file);
    if lp'length = 1 then
      declare
        p : Poly := lp(lp'first);
        dp : constant integer32 := Degree(p);
        b : Vector(1..n) := Random_Vector(1,n);
        v : Vector(1..n) := Random_Vector(1,n);
        tp : Vector(1..dp);
        fail : boolean;
      begin
        new_line;
        put_line("See the output file for results...");
        new_line;
        tstart(timer);
        Hypersurface_Witness_Points(file,n,dp,p,b,v,tp,fail);
        tstop(timer);
      end;
    else
      new_line;
      put("Is there a multi-homogeneous structure to respect ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        declare
          z : constant Partition := Read_Partition(n);
        begin
          new_line;
          put_line("See the output file for results...");
          new_line;
          tstart(timer);
          Multihomogeneous_Hypersurface_Solver(file,n,z,lp.all);
          tstop(timer);
        end;
      else
        new_line;
        put_line("See the output file for results...");
        new_line;
        declare
          b : Vector(1..n);
          w : VecVec(lp'range);
          sols : Array_of_Solution_Lists(1..n);
        begin
          tstart(timer);
          Total_Degree_Hypersurface_Solver(file,n,lp.all,b,w,sols);
          tstop(timer);
        end;
      end if;
    end if;
    new_line(file);
    print_times(file,timer,"Computation of Witness Points");
  end Driver_to_Hypersurface_Solvers;

end Drivers_to_Intrinsic_Solvers;
