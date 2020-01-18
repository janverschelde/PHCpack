with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Standard_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with Newton_Convolutions;

procedure ts_mtnewton is

-- DESCRIPTION :
--   Tests the development of Newton's method on power series
--   with the reverse mode of algorithmic differentation
--   and linearization to solve the matrix series equations,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Standard_Run
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in Standard_Complex_Vectors.Vector;
	        deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double precision on a solution sol
  --   of the system p, with power series of degree deg.

     use Standard_Speelpenning_Convolutions;

     c : constant Convolution_Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
     s : constant Link_to_System := Create(c,p'last,deg);
     dim : constant integer32 := sol'last;
     scf : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Newton_Convolutions.Series_Coefficients(sol,deg);
     deg,maxit,nbt : integer32 := 0;

   begin
     Add_Parameter_to_Constant(s);
     new_line;
     put("Give the number of iterations : "); get(maxit);
     put("Give the number of tasks : "); get(nbt);
   end Standard_Run;

   procedure Standard_Test is

   -- DESCRIPTION :
   --   Prompts the user for a polynomial system with solutions.

     lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
     sols : Standard_Complex_Solutions.Solution_List;
     nbr,dim : natural32;
     ls : Standard_Complex_Solutions.Link_to_Solution;
     deg : integer32;

   begin
     new_line;
     put_line("Reading a polynomial system with solutions ...");
     Standard_System_and_Solutions_io.get(lp,sols);
     nbr := Standard_Complex_Solutions.Length_Of(sols);
     ls := Standard_Complex_Solutions.Head_Of(sols);
     dim := natural32(ls.n);
     new_line;
     put("Read "); put(nbr,1); put(" solutions in dimension ");
     put(dim,1); put_line(".");
     new_line;
     put("Give the degree of the series : "); get(deg);
     Standard_Run(lp,ls.v,deg);
   end Standard_Test;

   procedure Main is

   -- DESCRIPTION :
   --   Launches the test.

   begin
     new_line;
     put_line("Testing Newton's method on power series ...");
     Standard_Test;
   end Main;

begin
  Main;
end ts_mtnewton;
