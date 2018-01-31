with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer64_Vectors;
with Standard_Integer_VecVecs;
with Integer_Permanents;
with Static_Columns_Queue;
with Multitasking;

package body Multitasking_Integer_Permanents is

  procedure Initialize
              ( mat : in Boolean_Matrices.Matrix;
                nbrows : in integer32; size : out integer32 ) is

    dim : constant integer32 := mat'last(1);
    nbstc : constant integer32
          := Integer_Permanents.Number_of_Start_Columns(dim,nbrows);
    stcol : Standard_Integer_VecVecs.VecVec(1..nbstc);
    cols : Standard_Integer_Vectors.Vector(mat'range(1))
         := (mat'range(1) => 0);
    cnts : Standard_Integer_Vectors.Vector(mat'range(1))
         := (mat'range(1) => 1);
    idx : integer32 := 0;

  begin
    Integer_Permanents.Start_Columns
      (mat'first(1),nbrows,mat,cols,cnts,idx,stcol);
    Static_Columns_Queue.Initialize(stcol(1..idx));
   -- put_line("The start column indices : ");
   -- for k in 1..idx loop
   --   put(stcol(k)); new_line;
   -- end loop;
    size := idx;
  end Initialize;

  procedure Initialize
              ( mat : in Standard_Integer_Matrices.Matrix;
                nbrows : in integer32; size : out integer32 ) is

    dim : constant integer32 := mat'last(1);
    nbstc : constant integer32
          := Integer_Permanents.Number_of_Start_Columns(dim,nbrows);
    stcol : Standard_Integer_VecVecs.VecVec(1..nbstc);
    cols : Standard_Integer_Vectors.Vector(mat'range(1))
         := (mat'range(1) => 0);
    cnts : Standard_Integer_Vectors.Vector(mat'range(1))
         := (mat'range(1) => 1);
    idx : integer32 := 0;

  begin
    Integer_Permanents.Start_Columns
      (mat'first(1),nbrows,mat,cols,cnts,idx,stcol);
    Static_Columns_Queue.Initialize(stcol(1..idx));
   -- put_line("The start column indices : ");
   -- for k in 1..idx loop
   --   put(stcol(k)); new_line;
   -- end loop;
    size := idx;
  end Initialize;

  function Permanent ( mat : Boolean_Matrices.Matrix;
                       nbrows,ntasks : integer32; verbose : boolean )
                     return integer64 is

    nbjobs : integer32;
    factors : Standard_Integer64_Vectors.Vector(1..ntasks)
            := (1..ntasks => 0);

    procedure Next ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Retrieves the next job from the queue and computes
    --   the corresponding factor in the permanent.
 
      nextcols : Standard_Integer_Vectors.Link_to_Vector;
      colidx : Standard_integer_Vectors.Vector(mat'range(2));
      cnts : Standard_integer_Vectors.Vector(mat'range(2));
      per : integer64;

      use Standard_Integer_Vectors;

    begin
      loop
        nextcols := Static_Columns_Queue.Next_Columns;
        exit when (nextcols = null);
        for i in nextcols'range loop
          colidx(i) := nextcols(i);
        end loop;
        cnts := (mat'range(2) => 1);
        for col in 1..nbrows loop
          cnts(colidx(col)) := 0;
        end loop;
        per := 0;
        Integer_Permanents.Permanent(nbrows+1,mat,colidx,cnts,per);
       -- put(colidx); put(" : "); put(per,1); new_line;
        factors(i) := factors(i) + per;
      end loop;
    end Next;
    procedure silent is new Multitasking.Silent_Workers(Next);
    procedure report is new Multitasking.Reporting_Workers(Next);

  begin
    Initialize(mat,nbrows,nbjobs);
    if verbose then
      put("Initialized the queue with "); put(nbjobs,1);
      put_line(" jobs.");
    end if;
    if verbose
     then report(ntasks);
     else silent(ntasks);
    end if;
    return Standard_Integer64_Vectors.Sum(factors);
  end Permanent;

  function Permanent ( mat : Standard_Integer_Matrices.Matrix;
                       nbrows,ntasks : integer32; verbose : boolean )
                     return integer64 is
    nbjobs : integer32;
    factors : Standard_Integer64_Vectors.Vector(1..ntasks)
            := (1..ntasks => 0);

    procedure Next ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Retrieves the next job from the queue and computes
    --   the corresponding factor in the permanent.
 
      nextcols : Standard_Integer_Vectors.Link_to_Vector;
      colidx : Standard_integer_Vectors.Vector(mat'range(2));
      cnts : Standard_integer_Vectors.Vector(mat'range(2));
      per : integer64;

      use Standard_Integer_Vectors;

    begin
      loop
        nextcols := Static_Columns_Queue.Next_Columns;
        exit when (nextcols = null);
        for i in nextcols'range loop
          colidx(i) := nextcols(i);
        end loop;
        cnts := (mat'range(2) => 1);
        for col in 1..nbrows loop
          cnts(colidx(col)) := 0;
        end loop;
        per := 0;
        Integer_Permanents.Permanent(nbrows+1,mat,colidx,cnts,per);
       -- put(colidx); put(" : "); put(per,1); new_line;
        factors(i) := factors(i) + per;
      end loop;
    end Next;
    procedure silent is new Multitasking.Silent_Workers(Next);
    procedure report is new Multitasking.Reporting_Workers(Next);

  begin
    Initialize(mat,nbrows,nbjobs);
    if verbose then
      put("Initialized the queue with "); put(nbjobs,1);
      put_line(" jobs.");
    end if;
    if verbose
     then report(ntasks);
     else silent(ntasks);
    end if;
    return Standard_Integer64_Vectors.Sum(factors);
  end Permanent;

end Multitasking_Integer_Permanents;
