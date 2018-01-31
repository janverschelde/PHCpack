with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User; 
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Boolean_Matrices;
with Boolean_Matrices_io;                use Boolean_Matrices_io;
with Standard_Random_Matrices;
with Integer_Permanents;                 use Integer_Permanents;

procedure ts_mtperm is

-- DESCRIPTION :
--   Interactive test on the permanent computation of an integer matrix.

  function Random_Boolean_Matrix
             ( dim : natural32 ) return Boolean_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a random Boolean matrix of the given dimension.

    idm : constant integer32 := integer32(dim);
    res : Boolean_Matrices.Matrix(1..idm,1..idm);
    rnd : integer32;

  begin
    for i in 1..idm loop
      for j in 1..idm loop
        rnd := Standard_Random_Numbers.Random(0,1);
        res(i,j) := (rnd = 1);
      end loop;
    end loop;
    return res;
  end Random_Boolean_Matrix;

  procedure Integer_Random_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random integer matrix of the given dimension dim
  --   and computes its permanent.

    idm : constant integer32 := integer32(dim);
    mat : constant Standard_Integer_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,0,9);
    per : integer64;
    ans : character;
    timer : Timing_Widget;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Do you want output of the factors ? (y/n) ");
    Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' 
     then per := Permanent(standard_output,mat);
     else per := Permanent(mat);
    end if;
    tstop(timer);
    put("The permanent : "); put(per,1); new_line;
    new_line;
    print_times(standard_output,timer,"Permanent Computation");
  end Integer_Random_Test;

  procedure Boolean_Random_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random Boolean matrix of the given dimension dim
  --   and computes its permanent.

    idm : constant integer32 := integer32(dim);
    mat : Boolean_Matrices.Matrix(1..idm,1..idm)
        := Random_Boolean_Matrix(dim);
    per : integer64;
    ans : character;
    timer : Timing_Widget;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Do you want output of the factors ? (y/n) ");
    Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' 
     then per := Permanent(standard_output,mat);
     else per := Permanent(mat);
    end if;
    tstop(timer);
    put("The permanent : "); put(per,1); new_line;
    new_line;
    print_times(standard_output,timer,"Permanent Computation");
  end Boolean_Random_Test;

  procedure Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and asks whether the matrix is Boolean or not.

    dim : natural32 := 0;
    ans : character;

  begin
    put("Give the dimension : "); get(dim);
    put("Generate a Boolean matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Boolean_Random_Test(dim);
     else Integer_Random_Test(dim);
    end if;
  end Random_Test;

  procedure Start_Integer_Random_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random integer matrix of the given dimension dim
  --   and computes its permanent, starting with the column indices
  --   for a number of rows first.

    idm : constant integer32 := integer32(dim);
    mat : constant Standard_Integer_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,0,9);
    nbr : integer32 := 0;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Give the number of rows : "); get(nbr);
    declare
      nbstc : constant integer32
            := Number_of_Start_Columns(integer32(dim),nbr);
      stcol : Standard_Integer_VecVecs.VecVec(1..nbstc);
      cols : Standard_Integer_Vectors.Vector(mat'range(1))
           := (mat'range(1) => 0);
      cnts : Standard_Integer_Vectors.Vector(mat'range(1))
           := (mat'range(1) => 1);
      idx : integer32 := 0;
      per : integer64;
      ans : character;
    begin
      Start_Columns(mat'first(1),nbr,mat,cols,cnts,idx,stcol);
      put_line("The start column indices : ");
      for k in 1..idx loop
        put(stcol(k)); new_line;
      end loop;
      put("Do you want output of the factors ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then per := Start_Permanent(standard_output,nbr+1,stcol(1..idx),mat);
       else per := Start_Permanent(nbr+1,stcol(1..idx),mat);
      end if;
      put("The permanent : "); put(per,1); new_line;
      per := Permanent(mat);
      put("The permanent : "); put(per,1); new_line;
    end;
  end Start_Integer_Random_Test;

  procedure Start_Boolean_Random_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random Boolean matrix of the given dimension dim
  --   and computes its permanent, starting with the column indices
  --   for a number of rows first.

    idm : constant integer32 := integer32(dim);
    mat : Boolean_Matrices.Matrix(1..idm,1..idm)
        := Random_Boolean_Matrix(dim);
    nbr : integer32 := 0;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Give the number of rows : "); get(nbr);
    declare
      nbstc : constant integer32
            := Number_of_Start_Columns(integer32(dim),nbr);
      stcol : Standard_Integer_VecVecs.VecVec(1..nbstc);
      cols : Standard_Integer_Vectors.Vector(mat'range(1))
           := (mat'range(1) => 0);
      cnts : Standard_Integer_Vectors.Vector(mat'range(1))
           := (mat'range(1) => 1);
      idx : integer32 := 0;
      per : integer64;
      ans : character;
    begin
      Start_Columns(mat'first(1),nbr,mat,cols,cnts,idx,stcol);
      put_line("The start column indices : ");
      for k in 1..idx loop
        put(stcol(k)); new_line;
      end loop;
      put("Do you want output of the factors ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then per := Start_Permanent(standard_output,nbr+1,stcol(1..idx),mat);
       else per := Start_Permanent(nbr+1,stcol(1..idx),mat);
      end if;
      put("The permanent : "); put(per,1); new_line;
      per := Permanent(mat);
      put("The permanent : "); put(per,1); new_line;
    end;
  end Start_Boolean_Random_Test;

  procedure Start_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and asks whether the matrix is Boolean or not.

    dim : natural32 := 0;
    ans : character;

  begin
    put("Give the dimension : "); get(dim);
    put("Generate a Boolean matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Start_Boolean_Random_Test(dim);
     else Start_Integer_Random_Test(dim);
    end if;
  end Start_Random_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test a permament computation : ");
    put_line("  1. test on a randomly generated matrix;");
    put_line("  2. compute start columns for the first rows.");
    put("Type 1 or 2 to choose a test : ");
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Random_Test;
      when '2' => Start_Random_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtperm;
