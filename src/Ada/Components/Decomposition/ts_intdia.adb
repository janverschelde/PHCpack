with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Plane_Operations;         use Standard_Plane_Operations;
with Standard_Cascading_Planes;         use Standard_Cascading_Planes;
with Intrinsic_Diagonal_Continuation;   use Intrinsic_Diagonal_Continuation;

procedure ts_intdia is

-- DESCRIPTION :
--   Interactive development of an intrinsic diagonal homotopy method
--   to intersect positive dimensional varieties.

  procedure Test_Linear_Cascade ( n,a,b,h0 : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction of linear spaces used to intersect
  --   two components of dimensions a and b in n-space.

    n2 : constant integer32 := 2*n;
    apb : constant integer32 := a+b;
    m : constant integer32 := n2-apb;
    AA : constant Matrix(1..apb,1..n)
       := Standard_Random_Matrices.Random_Matrix(natural32(apb),natural32(n));
    dA : constant Matrix(1..apb,1..n2) := Double_to_Diagonal(AA);
    BB : constant Matrix(1..apb,1..n)
       := Standard_Random_Matrices.Random_Matrix(natural32(apb),natural32(n));
    CC : constant Matrix(1..n,1..n2)
       := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n2));
    d : constant Vector(1..n) := Random_Vector(1,n);
    offset : constant Vector(1..n2) := Compute_Offset(CC,d);
    L1equ,L2equ : Matrix(1..apb,0..n2);
    L1gen,L2gen : Matrix(1..n2,0..m);
    res : double_float;

  begin
    L1equ := Target_Space(n,n2,apb,b,dA,BB,CC,d);
    L1gen := Generators(L1equ);
    put_line("The offset vector :"); put_line(offset);
   -- put_line("Verification of the offset vector at 1st plane: ");
   -- Evaluate(Standard_Output,L1equ,offset,res);
   -- put("  residual : "); put(res,3); new_line;
    Evaluate(L1equ,offset,res);
    put("Value of offset at 1st plane : "); put(res,3); new_line;
    Shift_Offset(L1gen,offset);
    for h in reverse h0..b-1 loop
      new_line;
      put("Generating space at dimension "); put(h,1); put_line("...");
      new_line;
      L2equ := Target_Space(n,n2,apb,h,dA,BB,CC,d);
      L2gen := Generators(L2equ);
     -- put_line("Verification of the offset vector at 2nd plane: ");
     -- Evaluate(Standard_Output,L2equ,offset,res);
     -- put("  residual : "); put(res,3); new_line;
      Evaluate(L2equ,offset,res);
      put("Value of offset at 2nd plane : "); put(res,3); new_line;
      Shift_Offset(L2gen,offset);
      Intersect(L1equ,L2equ,L1gen,L2gen);
     -- put_line("Verification of 1st plane");
     -- Evaluate(Standard_Output,L1equ,L1gen,res);
     -- put("  residual : "); put(res,3); new_line;
      Evaluate(L1equ,L1gen,res);
      put("Value of 1st generators at 1st plane :");
      put(res,3); new_line;
     -- put_line("Verification of 2nd plane");
     -- Evaluate(Standard_Output,L2equ,L2gen,res);
     -- put("  residual : "); put(res,3); new_line;
      Evaluate(L2equ,L2gen,res);
      put("Value of 2nd generators at 2nd plane :");
      put(res,3); new_line;
      put_line("Evaluation of 2nd generators at 1st plane");
      Evaluate(Standard_Output,L1equ,L2gen,res);
      put_line("Evaluation of 1st generators at 2nd plane");
      Evaluate(Standard_Output,L2equ,L1gen,res);
      L1equ := L2equ; L1gen := L2gen;
    end loop;
  end Test_Linear_Cascade;

  procedure Cascade_of_Linear_Spaces is

  -- DESCRIPTION :
  --   Sets up the cascade of linear spaces used to compute witness
  --   sets for all components of the intersection of two varieties.

    n,a,b,h0 : integer32 := 0;

  begin
    new_line;
    put("Give n (dimension of ambient space): "); get(n);
    put("Give a (dimension of 1st component): "); get(a);
    put("Give b (dimension of 2nd component): "); get(b);
    h0 := Minimal_Intersection_Dimension(n,a,b);
    put("The minimal dimension of intersection h0 = ");
    put(h0,1); new_line;
    Test_Linear_Cascade(n,a,b,h0);
  end Cascade_of_Linear_Spaces;

  procedure Main is

   -- ans : character;

  begin
   -- new_line;
   -- put_line("MENU to test intrinsic diagonal homotopies:");
   -- put_line("  1. set up cascade of linear spaces; or");
   -- put_line("  2. intersect two positive dimensional solution components.");
   -- put("Type 1 or 2 to select : ");
   -- Ask_Alternative(ans,"12");
   -- if ans = '1'
   --  then Cascade_of_Linear_Spaces;
   --  else Main_Driver_to_Intersect_Varieties;
   -- end if;
    new_line;
    put_line("Setting up a cascade of linear spaces ...");
    Cascade_of_Linear_Spaces;
  end Main;

begin
  Main;
end ts_intdia;
