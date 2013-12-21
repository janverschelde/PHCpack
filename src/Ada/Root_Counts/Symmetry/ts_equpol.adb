with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Symmetry_Group;                     use Symmetry_Group;
with Symmetry_Group_io;
with Symbolic_Symmetry_Group_io;
with Drivers_for_Symmetry_Group_io;      use Drivers_for_Symmetry_Group_io;
with Equivariant_Polynomial_Systems;     use Equivariant_Polynomial_Systems;

procedure ts_equpol is

-- DESCRIPTION :
--   Test on the (G,V,W)-symmetric polynomial systems.

  lp : Link_to_Poly_Sys;
  n : natural32;
  g,v,w : List_of_Permutations;
  allperms,notsym,inva,equi : boolean;

begin
  new_line;
  put_line("Test on the (G,V,W)-symmetric polynomial systems.");
  new_line;
  get(lp);
  n := natural32(lp'last);
  Read_Permutation_Group(n,g,v,allperms);
  put_line("The symmetry group : ");
  Symbolic_Symmetry_Group_io.put(v);
  Act(v,lp.all,w,notsym,inva,equi);
  put_line("w:"); Symmetry_Group_io.put(w);
  if notsym
   then put_line("The system is not (G,V,W)-symmetric.");
  end if;
end ts_equpol;
