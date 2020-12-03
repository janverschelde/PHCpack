with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Permutations,Shuffle_Polynomials;   use Permutations,Shuffle_Polynomials; 
with Intrinsic_Witness_Sets_io;          use Intrinsic_Witness_Sets_io;
with Equation_by_Equation_Solvers;       use Equation_by_Equation_Solvers;

package body Main_Eqn_by_Eqn_Solvers is

  function Shuffle ( file : file_type; p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);
    ans : character;

  begin
    new_line;
    put_line("MENU to shuffle polynomials in the system : ");
    put_line("  0. leave the polynomials in their given order;");
    put_line("  1. sort the polynomials in increasing order of degrees;");
    put_line("  2. sort the polynomials in decreasing order of degrees;");
    put_line("  3. apply a random permutation to rearrange the order;");
    put_line("  4. give your own permutation to shuffle the polynomials.");
    put("Type 0, 1, 2, 3, or 4 for your shuffle selection : ");
    Ask_Alternative(ans,"01234");
    new_line(file);
    case ans is
      when '1' =>
        res := Increasing_Degree_Sort(p);
        put_line(file,"The system sorted in increasing order of degrees :");
        put(file,res); new_line(file);
      when '2' =>
        res := Decreasing_Degree_Sort(p);
        put_line(file,"The system sorted in decreasing order of degrees :");
        put(file,res); new_line(file);
      when '3' =>
        declare
          pp : constant Permutation(p'range) := Random_Permutation(p'last);
        begin
          put("Applying the permutation");
          put(Standard_Integer_Vectors.Vector(pp)); new_line;
          res := Permute(p,pp);
          put(file,"The system after random permutation");
          put(file,Standard_Integer_Vectors.Vector(pp));
          put_line(file," :"); put(file,res); new_line(file);
        end;
      when '4' =>
        declare
          pp : Permutation(p'range);
        begin
          put("Reading a permutation of 1 to ");
          put(p'last,1); put_line(" ...");
          get(Standard_Integer_Vectors.Vector(pp));
          put("Applying the permutation ");
          put(Standard_Integer_Vectors.Vector(pp)); new_line;
          res := Permute(p,pp);
          put(file,"The system after given permutation");
          put(file,Standard_Integer_Vectors.Vector(pp));
          put_line(file," :"); put(file,res); new_line(file);
        end;
      when others =>
         put_line(file,"Given order of polynomial equations taken.");
         new_line(file);
         return p;
    end case;
    return res;
  end Shuffle;

  function Maximum ( a,b : integer32 ) return integer32 is
  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Maximum;

  procedure Shuffle_Polynomials_and_Solve
              ( file : in file_type; filename : in string;
                p : in Poly_Sys ) is

    ans : character;
    have_stone,want_stone,report : boolean;
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    q : Poly_Sys(p'range);
    md : constant integer32 := Maximum(p'last,n);
    witset : Array_of_Solution_Lists(1..md);
    planes : VecMat(1..md);
    s : Link_to_Poly_Sys;
    s_sols : Solution_List;
    s_dim,oc : natural32;
    k : integer32;
    s_p : Link_to_Matrix;

  begin
    new_line;
    put("Do you have already a witness set for the first equations ? ");
    Ask_Yes_or_No(ans);
    have_stone := (ans = 'y');
    if have_stone then
      Read_Witness_Stone(s,s_dim,natural32(k),s_sols,s_p);
      planes(k) := s_p;    -- planes(s_dim) := s_p;
      witset(k) := s_sols; -- witset(s_dim) := s_sols;
      q := p;
    else
      q := Shuffle(file,p);
      k := 0;
    end if;
    put(file,natural32(q'last),natural32(n),q);
    new_line;
    put("Do you want stepping stones of indetermediate witness sets ? ");
    Ask_Yes_or_No(ans);
    want_stone := (ans = 'y');
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(file);
    put_line(file,"Running the equation-by-equation solver ...");
    new_line(file);
    P_Solve_Equation_by_Equation
      (file,filename,report,want_stone,q'last,n,k,q,witset,planes);
  end Shuffle_Polynomials_and_Solve;

  procedure Solver_using_Generics
              ( file : in file_type; filename : in string;
                p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Prepares the executable versions to instantiate the
  --   generic equation-by-equation solver.

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    ep : constant Eval_Poly_Sys(p'range) := Create(p);
    jm : constant Jaco_Mat(p'range,1..n) := Create(p);
    jf : constant Eval_Jaco_Mat(p'range,1..n) := Create(jm);

    function Eval ( k : integer32; x : Vector ) return Complex_Number is
    begin
      return Eval(ep(k),x);
    end Eval;

    function Diff ( k : integer32; x : Vector ) return Vector is

      res : Vector(1..n);

    begin
      for i in 1..n loop
        res(i) := Eval(jf(k,i),x);
      end loop;
      return res;
    end Diff;

    procedure Solve is new G_Solve_Equation_by_Equation(Eval,Diff);

    md : constant integer32 := Maximum(p'last,n);
    witset : Array_of_Solution_Lists(1..md);
    planes : VecMat(1..md);
    d : Standard_Natural_Vectors.Vector(p'range);
    ans : character;
    have_stone,want_stone,report : boolean;
    s_dim,oc : natural32;
    k : integer32;
    s_sols : Solution_List;
    s_p : Link_to_Matrix;

  begin
    put(file,natural32(p'last),natural32(n),p);
    new_line;
    put("Do you have already a witness set for the first equations ? ");
    Ask_Yes_or_No(ans);
    have_stone := (ans = 'y');
    if have_stone then
      Read_Witness_Stone(s_dim,natural32(k),s_sols,s_p);
      witset(k) := s_sols; -- witset(s_dim) := s_sols;
      planes(k) := s_p;    -- planes(s_dim) := s_p;
    else
      k := 0;
    end if;
    new_line;
    put("Do you want stepping stones of indetermediate witness sets ? ");
    Ask_Yes_or_No(ans);
    want_stone := (ans = 'y');
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(file);
    put_line(file,"Running the equation-by-equation solver ...");
    new_line(file);
    for i in p'range loop
      d(i) := natural32(Degree(p(i)));
    end loop;
    Solve(file,filename,report,want_stone,p'last,n,k,d,witset,planes);
  end Solver_using_Generics;

  procedure Read_System ( name : in string; p : out Link_to_Poly_Sys ) is

    file : file_type;

  begin
    Open(file,in_file,name);
    get(file,p);
  exception
    when others => new_line;
                   put_line("Exception occurred with file " & name & ".");
                   get(p);
  end Read_System;

  procedure Interactive_Create_Output_File
               ( file : in out file_type; name : in string;
                 new_name : out Link_to_String ) is

    procedure Ask_File is
    begin
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file,new_name);
    end Ask_File;

  begin
    if name = ""
     then Ask_File;
     else Create_Output_File(file,name,new_name); 
    end if;
  end Interactive_Create_Output_File;

  procedure Main ( infilename,outfilename : in string ) is

    file : file_type;
    lp : Link_to_Poly_Sys;
    name : Link_to_String;

  begin
    if infilename /= ""
     then Read_System(infilename,lp);
     else new_line; get(lp);
    end if;
    Interactive_Create_Output_File(file,outfilename,name);
    if name /= null
     then Shuffle_Polynomials_and_Solve(file,name.all,lp.all);
     else Shuffle_Polynomials_and_Solve(file,outfilename,lp.all);
    end if;
  end Main;

end Main_Eqn_by_Eqn_Solvers;
