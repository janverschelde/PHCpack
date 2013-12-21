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
with Standard_Exponent_Transformations;  use Standard_Exponent_Transformations;
with Span_of_Supports;                   use Span_of_Supports;
with Transformation_of_Supports;

procedure ts_supspan is

-- DESCRIPTION :
--   Computes the span of the supports of a polynomial system.
--   If the span is less than the dimension of the ambient space,
--   then the generators of the normals to the hyperplanes that
--   contain the supports are computed.
--   The test is to generate a polynomial system with random supports
--   of a given positive corank, write the system to file and then
--   relaunch this test giving the generated system and tranform the
--   system so as to eliminate as many variables as the corank.

  procedure Transform
              ( support : in List; dim : in integer32;
                output : in boolean; rnksup : in natural32;
                normals : in Standard_Integer_Matrices.Matrix;
                fail : out boolean;
                transfo : out Standard_Integer_Matrices.Link_to_Matrix ) is

  -- DESCRIPTION :
  --   Computes the coordinate transformation for a support list of rank 
  --   rnksup < dim, with given normals to the linear span of the support.

  -- REQUIRED : rnksup < dim.

  -- ON ENTRY :
  --   support  a support set with rank of its span less than dim;
  --   dim      all vectors in the support have range 1..dim;
  --   output   true if intermediate output is desired;
  --   rnksup   the rank of the span of the support;
  --   normal   a basis of vector perpendicular to the span of the support.

  -- ON RETURN :
  --   fail     true if no integer valued unimodular coordinate
  --            transformation is possible, false otherwise;
  --   transfo  if not fail, an unimodular coordinate transformation
  --            with in its first rows the columns of normals.

    M : Standard_Integer_Matrices.Matrix(1..dim,1..dim);
    trasup : List;

  begin
    if output
     then Unimodular_Coordinate_Transformation(standard_output,normals,fail,M);
     else Unimodular_Coordinate_Transformation(normals,fail,M);
    end if;
    if fail then
      put_line("no uni-rational coordinate transformation possible ...");
    else
      transfo := new Matrix'(M);
      put_line("The unimodular coordinate transformation : ");
      put(M);
      trasup := Transformation_of_Supports.Transform(support,M);
      put_line("The transformed support : "); put(trasup);
    end if;
  end Transform;

  procedure Rank_of_Span
              ( support : in List; dim : in integer32;
                output : in boolean; suprnk : out natural32;
                normals : out Standard_Integer_Matrices.Link_to_Matrix;
                fail : out boolean;
                transfo : out Standard_Integer_Matrices.Link_to_Matrix ) is

  -- DESCRIPTION :
  --   Computes the rank of the span of the points in support.
  --   If output, then additional information is printed.
  
  -- ON ENTRY :
  --   support  a list of points;
  --   dim      every vector in the support is of range 1..dim;
  --   output   if true, then additional information is printed.

  -- ON RETURN :
  --   suprnk   rank of the support;
  --   normals  if suprnk < dim, then normals contains in its columns a
  --            basis of vectors perpendicular to the span of the supports.
  --   fail     true if no integer unimodular coordinate transformation
  --            could be found, false otherwise;
  --   transfo  if not fail, then transfo is a unimodular transformation
  --            with in its first rows the columns of the normals.

  begin
    new_line;
   -- suprnk := Rank32_of_Support(support,dim,output);
    suprnk := Rank64_of_Support(support,dim,output);
    put("The computed rank : "); put(suprnk,1); new_line;
    if suprnk < natural32(dim) then
      Rank_of_Support(support,dim,output,suprnk,normals);
      put_line("A basis of vectors normal to the span :");
      put(normals.all);
      Transform(support,dim,output,suprnk,normals.all,fail,transfo);
    end if;
  end Rank_of_Span;

  procedure Rank_of_Span
              ( supports : in Array_of_Lists; dim : in integer32;
                output : in boolean; suprnk : out natural32;
                normals : out Standard_Integer_Matrices.Link_to_Matrix;
                fail : out boolean;
                transfo : out Standard_Integer_Matrices.Link_to_Matrix ) is

  -- DESCRIPTION :
  --   Computes the rank of the span of the points in the supports.
  --   If output, then additional information is printed.
  
  -- ON ENTRY :
  --   supports is a tuple of lists of points;
  --   dim      every vector in the supports is of range 1..dim;
  --   output   if true, then additional information is printed.

  -- ON RETURN :
  --   suprnk   rank of the span of the supports;
  --   normals  if suprnk < dim, then normals contains in its columns a
  --            basis of vectors perpendicular to the span of the supports.
  --   fail     true if no integer unimodular coordinate transformation
  --            could be found, false otherwise;
  --   transfo  if not fail, then transfo is a unimodular transformation
  --            with in its first rows the columns of the normals.

    c : List := Cayley_Embedding(supports);

  begin
    Rank_of_Span(c,2*dim-1,true,suprnk,normals,fail,transfo);
    if not fail then
      put_line("The transformation : "); put(transfo.all);
      declare
        M : constant Standard_Integer_Matrices.Matrix(1..dim,1..dim)
          := Remove_Cayley_Embedding(transfo.all,dim);
        N : constant Standard_Integer_Matrices.Matrix(1..dim,normals'range(2))
          := Remove_Cayley_Rows(normals.all,dim);
      begin
        Clear(transfo);
        transfo := new Standard_Integer_Matrices.Matrix'(M);
        Clear(normals);
        normals := new Standard_Integer_Matrices.Matrix'(N);
      end;
    end if;
    Clear(c);
  end Rank_of_Span;

  procedure Generate_Support is

  -- DESCRIPTION :
  --   Prompts the user for the ambient dimension, the number of points,
  --   lower, upper bounds on the random numbers, and a bound on the rank
  --   of the point configuration stored in a list of points.
  --   After the generation of the support, its rank will be determined.

    dim,nbp,low,upp,rnk : integer32 := 0;
    support : List;
    ans : character;
    rnksup : natural32;
    fail : boolean;
    normals,transfo : Link_to_Matrix;

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
    Rank_of_Span(support,dim,ans='y',rnksup,normals,fail,transfo);
  end Generate_Support;

  procedure Write_to_File ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and writes the system p
  --   to that file.

    file : file_type;
 
  begin
    new_line;
    put_line("Reading the name of a file ...");
    Read_Name_and_Create_File(file);
    put(file,natural32(p'last),p);
    close(file);
  end Write_to_File;

  procedure Generate_Support_and_System is

  -- DESCRIPTION :
  --   Prompts the user for the ambient dimension, the number of points,
  --   lower, upper bounds on the random numbers, and a bound on the rank
  --   of the point configuration stored in a list of points.
  --   After the generation of the support, a polynomial system
  --   with random coefficients and the generated support will
  --   be created and written to file.

    dim,nbp,low,upp,rnk : integer32 := 0;
    support : List;

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
    declare
      use Supports_of_Polynomial_Systems;
      p : constant Laur_Sys := Random_Complex_Laurent_System(dim,support);
      ans : character;
    begin
      put_line("A random Laurent system : "); put(p);
      new_line;
      put("Write to file ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_to_File(p);
      end if;
    end;
  end Generate_Support_and_System;

  procedure Generate_Tuple_of_Supports is

  -- DESCRIPTION :
  --   Prompts the user for the ambient dimension and the number of
  --   different supports followed by their type of mixture.

    dim,mix,low,upp,rnk : integer32 := 0;
    rnksup : natural32;
    fail : boolean;
    normals,transfo : Link_to_Matrix;

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
      Rank_of_Span(cayemb,dim+mix-1,true,rnksup,normals,fail,transfo);
    end;
  end Generate_Tuple_of_Supports;

  procedure Generate_Tuple_of_Supports_and_System is

  -- DESCRIPTION :
  --   Prompts the user for the ambient dimension and the number of
  --   different supports followed by their type of mixture.

    dim,mix,low,upp,rnk : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the ambient dimension : "); get(dim);
    put("Give the number of distinct supports : "); get(mix);
    declare
      use Supports_of_Polynomial_Systems;
      nbp,mixtype : Vector(1..mix);
      supports : Array_of_Lists(1..mix);
      p : Laur_Sys(1..dim);
    begin
      for i in 1..mix loop
         put("Give the number of points in support "); put(i,1);
         put(" : "); nbp(i) := 0; get(nbp(i));
         put("Give the number of occurrences of support "); put(i,1);
         put(" : "); mixtype(i) := 0; get(mixtype(i));
      end loop;
      if Sum(mixtype) = dim then
        put_line("Check sum on type of mixture is okay.");
      else
        put("Check sum on mixture type is ");
        put(Sum(mixtype),1); put(" /= "); put(dim,1); put_line(".");
        put_line("Please try again."); return;
      end if;
      put("Give a lower bound on the coordinates : "); get(low);
      put("Give an upper bound on the coordinates : "); get(upp);
      put("Give a bound on the rank : "); get(rnk);
      supports := Random_Tuple_of_Supports(dim,mix,low,upp,rnk,nbp);
      put_line("The generated tuple of supports :"); put(supports);
      p := Random_Complex_Laurent_System(dim,mixtype,supports);
      put_line("A random Laurent system : "); put(p);
      new_line;
      put("Write to file ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_to_File(p);
      end if;
    end;
  end Generate_Tuple_of_Supports_and_System;

  procedure Given_Support is

  -- DESCRIPTION :
  --   Prompts the user for the dimension and the number of points
  --   to read from standard input and stored as a list of points.
  --   After the reading of the support, its rank will be determined.

    dim,nbp : integer32 := 0;
    rnksup : natural32;
    fail : boolean;
    normals,transfo : Link_to_Matrix;
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
    Rank_of_Span(support,dim,ans='y',rnksup,normals,fail,transfo);
  end Given_Support;

  procedure Span_of_System ( p : in Laur_Sys; n : in integer32 ) is

  -- DESCRIPTION :
  --   Computes the span of the supports of p.
  --   The ambient dimension is n.

    s : Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);
    c : List := Cayley_Embedding(s);
    rnksup : natural32;
    fail : boolean;
    normals,transfo : Link_to_Matrix;

  begin
    new_line;
    put_line("The supports : "); put(s);
    put_line("The Cayley embedding : "); put(c);
    Rank_of_Span(c,2*n-1,true,rnksup,normals,fail,transfo);
    if not fail then
      put_line("The transformation : "); put(transfo.all);
      declare
        M : constant Standard_Integer_Matrices.Matrix(1..n,1..n)
          := Remove_Cayley_Embedding(transfo.all,n);
        tp : Laur_Sys(p'range);
      begin
        put_line("-> the Cayley embedding removed : "); put(M);
        tp := Transformation_of_Supports.Transform(p,M);
        put_line("The transformed system : "); put(tp);
      end;
    end if;
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
  --   or generating a random support or systems with random supports.

    ans : character;

  begin
    new_line;
    put_line("Computing the rank of a support set ...");
    put_line("  1. enter a Laurent polynomial system;");
    put_line("  2. generate one random support;");
    put_line("  3. give your own support set;");
    put_line("  4. generate a tuple of random support sets;");
    put_line
      ("  5. make polynomial system with one positive corank support;");
    put_line
      ("  6. make polynomial system with tuple of positive corank supports.");
    put("Type 1, 2, 3, 4, 5, or 6 to make your choice : ");
    Ask_Alternative(ans,"123456");
    case ans is
       when '1' => Read_System;
       when '2' => Generate_Support;
       when '3' => Given_Support;
       when '4' => Generate_Tuple_of_Supports;
       when '5' => Generate_Support_and_System;
       when '6' => Generate_Tuple_of_Supports_and_System;
       when others => null;
    end case;
  end Main;

begin
  Main;
end ts_supspan;
