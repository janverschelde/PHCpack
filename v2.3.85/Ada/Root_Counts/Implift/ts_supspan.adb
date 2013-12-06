with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Supports_of_Polynomial_Systems;
with Span_of_Supports;                   use Span_of_Supports;

procedure ts_supspan is

-- DESCRIPTION :
--   Computes the span of the supports of a polynomial system.
--   If the span is less than the dimension of the ambient space,
--   then the generators of the normals to the hyperplanes that
--   contain the supports are computed.

  procedure Rank ( support : in List; dim : in integer32;
                   output : in boolean ) is

  -- DESCRIPTION :
  --   Computes the rank of the span of the points in support.
  --   If output, then additional information is printed.

    rnksup : natural32;
    normals : Link_to_Matrix;

  begin
    new_line;
   -- rnksup := Rank32_of_Support(support,dim,output);
    rnksup := Rank64_of_Support(support,dim,output);
    put("The computed rank : "); put(rnksup,1); new_line;
    if rnksup < natural32(dim) then
      Rank_of_Support(support,dim,output,rnksup,normals);
      put_line("A basis of vectors normal to the span :");
      put(normals.all);
    end if;
  end Rank;

  procedure Generate_Support is

  -- DESCRIPTION :
  --   Prompts the user for the ambient dimension, the number of points,
  --   lower, upper bounds on the random numbers, and a bound on the rank
  --   of the point configuration stored in a list of points.
  --   After the generation of the support, its rank will be determined.

    dim,nbp,low,upp,rnk : integer32 := 0;
    support : List;
    ans : character;

  begin
    new_line;
    put("Give the ambient dimension : "); get(dim);
    put("Give the number of points : "); get(nbp);
    put("Give a lower bound on the coordinates : "); get(low);
    put("Give an upper bound on the coordinates : "); get(upp);
    put("Give a bound on the rank : "); get(rnk);
    support := Random_Support(dim,nbp,low,upp,rnk);
    put_line("The generated support set : "); put(support);
    new_line;
    put("Do you want output ? (y/n) "); Ask_Yes_or_No(ans);
    Rank(support,dim,ans='y');
  end Generate_Support;

  procedure Generate_Tuple_of_Supports is

  -- DESCRIPTION :
  --   Prompts the user for the ambient dimension and the number of
  --   different supports followed by their type of mixture.

    dim,mix,low,upp,rnk : integer32 := 0;

  begin
    new_line;
    put("Give the ambient dimension : "); get(dim);
    put("Give the number of distinct supports : "); get(mix);
    declare
      nbp : Vector(1..mix);
      supports : Array_of_Lists(1..mix);
      cayemb : List;
    begin
      for i in 1..mix loop
         put("Give the number of points in support "); put(i,1);
         put(" : "); nbp(i) := 0; get(nbp(i));
      end loop;
      put("Give a lower bound on the coordinates : "); get(low);
      put("Give an upper bound on the coordinates : "); get(upp);
      put("Give a bound on the rank : "); get(rnk);
      supports := Random_Tuple_of_Supports(dim,mix,low,upp,rnk,nbp);
      put_line("The generated tuple of supports :"); put(supports);
      cayemb := Cayley_Embedding(supports);
      put_line("The Cayley embedding of the supports :"); put(cayemb);
      Rank(cayemb,dim+mix-1,true);
    end;
  end Generate_Tuple_of_Supports;

  procedure Given_Support is

  -- DESCRIPTION :
  --   Prompts the user for the dimension and the number of points
  --   to read from standard input and stored as a list of points.
  --   After the reading of the support, its rank will be determined.

    dim,nbp : integer32 := 0;
    support,support_last : List;
    ans : character;

  begin
    new_line;
    put("Give the ambient dimension : "); get(dim);
    put("Give the number of points : "); get(nbp);
    put("Give the coordinates as a ");
    put(nbp,1); put("-by-"); put(dim,1); put_line(" matrix :");
    declare
      point : Vector(1..dim);
    begin
      for i in 1..nbp loop
        get(point);
        Append(support,support_last,point);
      end loop;
    end;
    put_line("The support set : "); put(support);
    new_line;
    put("Do you want output ? (y/n) "); Ask_Yes_or_No(ans);
    Rank(support,dim,ans='y');
  end Given_Support;

  procedure Span_of_System ( p : in Laur_Sys; n : in integer32 ) is

  -- DESCRIPTION :
  --   Computes the span of the supports of p.
  --   The ambient dimension is n.

    s : Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);
    c : List := Cayley_Embedding(s);

  begin
    new_line;
    put_line("The supports : "); put(s);
    put_line("The Cayley embedding : "); put(c);
    Rank(c,2*n-1,true);
    for i in s'range loop
      Clear(s(i));
    end loop;
    Clear(c);
  end Span_of_System;

  procedure Read_System is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and computes the span of its support.

    lp : Link_to_Laur_Sys;
    nq,nv : integer32 := 0;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("-> read "); put(nq,1); put(" polynomials in ");
    put(nv,1); put_line(" variables");
    Span_of_System(lp.all,nv);
  end Read_System;

  procedure Main is

  -- DESCRIPTION :
  --   Gives the user a choice between reading a system
  --   or generating a random support.

    ans : character;

  begin
    new_line;
    put_line("Computing the rank of a support set...");
    put_line("  1. enter a Laurent polynomial system; or");
    put_line("  2. generate one a random support;");
    put_line("  3. give your own support set;");
    put_line("  4. generate a tuple of support sets.");
    put("Type 1, 2, 3, or 4 to make your choice : ");
    Ask_Alternative(ans,"1234");
    case ans is
       when '1' => Read_System;
       when '2' => Generate_Support;
       when '3' => Given_Support;
       when '4' => Generate_Tuple_of_Supports;
       when others => null;
    end case;
  end Main;

begin
  Main;
end ts_supspan;
