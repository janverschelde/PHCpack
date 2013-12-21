with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Exponent_Transformations;  use Standard_Exponent_Transformations;
with Supports_of_Polynomial_Systems;
with Span_of_Supports;
with Transformation_of_Supports;

procedure Driver_to_Rank_Supports
            ( outfile : in file_type;
              p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

  procedure Compute_Rank_of_Span
              ( supports : in Array_of_Lists; n : in natural32;
                output : in boolean; suprnk : out natural32;
                normals : out Standard_Integer_Matrices.Link_to_Matrix;
                fail : out boolean;
                transfo : out Standard_Integer_Matrices.Link_to_Matrix ) is

  -- DESCRIPTION :
  --   Computes the rank of the span of the points in the supports.
  --   If output, then additional information is printed.
  
  -- ON ENTRY :
  --   supports is a tuple of lists of points;
  --   n        every vector in the supports is of range 1..n;
  --   output   if true, then additional information is printed.

  -- ON RETURN :
  --   suprnk   rank of the span of the supports;
  --   normals  if suprnk < dim, then normals contains in its columns a
  --            basis of vectors perpendicular to the span of the supports.
  --   fail     true if no integer unimodular coordinate transformation
  --            could be found, false otherwise;
  --   transfo  if not fail, then transfo is a unimodular transformation
  --            with in its first rows the columns of the normals.

    use Span_of_Supports;
    c : List := Cayley_Embedding(supports);
    dim : constant integer32 := 2*integer32(n)-1;

  begin
    suprnk := Rank64_of_Support(c,dim,output);
   -- put("The computed rank : "); put(suprnk,1); new_line;
    if suprnk < natural32(dim) then
      Rank_of_Support(c,dim,output,suprnk,normals);
     -- put_line("A basis of vectors normal to the span :");
     -- put(normals.all);
      declare
        T : Standard_Integer_Matrices.Matrix(1..dim,1..dim);
      begin
        if output
         then Unimodular_Coordinate_Transformation
                (standard_output,normals.all,fail,T);
         else Unimodular_Coordinate_Transformation(normals.all,fail,T);
        end if;
        if fail then
          put_line("No integer valued unimodular coordinate transformation.");
        else
         -- put_line("The transformation : "); put(T);
          declare
            d : constant integer32 := integer32(n);
            M : constant Standard_Integer_Matrices.Matrix(1..d,1..d)
              := Remove_Cayley_Embedding(T,d);
            N : constant Standard_Integer_Matrices.Matrix
                           (1..d,normals'range(2))
              := Remove_Cayley_Rows(normals.all,d);
          begin
            transfo := new Standard_Integer_Matrices.Matrix'(M);
            Clear(normals);
            normals := new Standard_Integer_Matrices.Matrix'(N);
          end;
        end if;
      end;
    end if;
    suprnk := suprnk - n + 1;  -- Cayley embedding !
    Clear(c);
  end Compute_Rank_of_Span;

  procedure Main is

    n : constant natural32
       := Standard_Complex_Laurentials.Number_of_Unknowns(p(p'first));
    s : Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);
    tp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    rnksup : natural32;
    debug_info : constant boolean := false;
    fail : boolean;
    normals,transfo : Link_to_Matrix;

  begin
    new_line;
    Compute_Rank_of_Span(s,n,debug_info,rnksup,normals,fail,transfo);
    put("The rank of the span of the supports : "); put(rnksup,1); new_line;
    if rnksup < n then
      put_line("A basis of vectors normal to the span :");
      put(normals.all);
      if not fail then
        put_line("The unimodular coordinate transformation : ");
        put(transfo.all);
        tp := Transformation_of_Supports.Transform(p,transfo.all);
        put(outfile,natural32(tp'last),1); new_line(outfile);
        put(outfile,tp);
        new_line(outfile);
        put_line(outfile,"Unimodular coordinate transformation :");
        put(outfile,transfo.all);
      end if;
    end if;
    put(outfile,"Supports of system on input have rank ");
    put(outfile,rnksup,1); put_line(outfile,".");
    new_line;
    put_line("See the output file for results ...");
    new_line;
  end Main;

begin
  Main;
end Driver_to_Rank_Supports;
