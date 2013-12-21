with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Inner_Normal_Cones;                 use Inner_Normal_Cones;
with Normal_Cone_Intersections;          use Normal_Cone_Intersections;
with Normal_Cone_Intersections_io;       use Normal_Cone_Intersections_io;

procedure ts_conint is
 
-- DESCRIPTION :
--   Test on the operations of the package Normal_Cone_Intersections.

  procedure Compute_Intersection_Matrix
              ( supports : in Array_of_Lists; g : in List;
                x : Vector; i : in integer32 ) is

  -- DESCRIPTION :
  --   Computes and displays the intersection matrix on screen.
  --   Lists all complementary columns.

    ans : character;
    n1 : constant integer32 := supports'length - 1;
    mg : constant integer32 := integer32(Length_Of(g));
    nc : constant integer32 := Number_of_Cones(supports,i);
    ima : Intersection_Matrix(n1,mg,nc);

    procedure Write_Selection ( cols : in Vector; continue : out boolean ) is
    begin
      put("selected columns : "); put(cols); new_line;
      continue := true;
    end Write_Selection;
    procedure Write_Complementary_Columns is
      new Complementary_Columns(Write_Selection);

    procedure Check_Selection ( cols : in Vector; continue : out boolean ) is

      part : Array_of_Lists(cols'range);

    begin
      put("selected columns : "); put(cols); new_line;
      part := Partition(ima,cols,g);
      put_line("The partition of the set of generators : "); put(part);
      if Partition_in_Union(part,supports,i,cols)
       then put_line("The partition is contained in the union of cones.");
       else put_line("The partition is NOT contained in the union of cones.");
      end if;
     -- Deep_Clear(part);
      put("The point "); put(x); put(" is ");
      if Contained_in_Union(supports,i,g,ima,cols)   -- double check
       then put("contained in the union of columns ");
       else put("not contained in the union of columns ");
      end if;
      put(cols); put_line(".");
      continue := true;
    end Check_Selection;
    procedure Check_Complementary_Columns is
      new Complementary_Columns(Check_Selection);

  begin
    ima := Create(supports,g,i);
    put(ima);
    put_line("The complementary columns :");
    Write_Complementary_Columns(ima);
    put("Do you want to check the complementary columns ? (y/n) ");
    get(ans);
    if ans = 'y'
     then Check_Complementary_Columns(ima);
    end if;
  end Compute_Intersection_Matrix;

  procedure Test_Intersection_Matrix ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Allows the computation of intersection matrix of cones
  --   on the points in the support lists of the polynomial system.

    supp : constant Array_of_Lists(p'range) := Create(p);
    genx : List;
    ind : integer32 := 0;
    x : Vector(p'range) := (p'range => 0);
    ans : character;

  begin
    loop
      new_line;
      put("Give a vector : "); get(x);
      put("Give an index : "); get(ind);
      genx := Generators(supp(ind),x);
      put("The generators of the normal cone at "); put(x); put_line(" :");
      put(genx);
      Compute_Intersection_Matrix(supp,genx,x,ind);
      put("Do you want more tests ? (y/n) "); Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Intersection_Matrix;

  procedure Main is

    p : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the intersection of normal cones.");
    new_line;
    get(p);
    Test_Intersection_Matrix(p.all);
  end Main;

begin
  Main;
end ts_conint;
