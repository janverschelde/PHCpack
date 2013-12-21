with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Symmetric_Polyhedral_Continuation;  use Symmetric_Polyhedral_Continuation;

package body Symmetric_BKK_Bound_Solvers is

  function Symmetric_BKK_Solve ( p : Poly_Sys; sign : boolean )
                               return Solution_List is
  
    n : constant integer32 := p'length;
    points : Array_of_Lists(p'range) := Create(p);
    sols : Solution_List;
    mix,per : Link_to_Vector;
    mixsub : Mixed_Subdivision;
    lp,pp : Laur_Sys(p'range);
    low : constant Vector := (1..n => 0);
    upp : constant Vector := Adaptive_Lifting(points);

  begin
    Compute_Mixture(points,mix,per); Clear(per);
    declare
      lifted : Array_of_Lists(mix'range);
      nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                    := (mix'range => 0.0);
      file : file_type;
    begin
      Mixed_Coherent_Subdivision
        (n,mix.all,points,false,low,upp,lifted,nbsucc,nbfail,mixsub);
      pp := Polynomial_to_Laurent_System(p);
      lp := Perform_Lifting(n,mix.all,lifted,pp);
      Create(file,out_file,"/tmp/brol");
      sols := Symmetric_Mixed_Solve(file,sign,lp,mixsub,n,mix.all);
      Deep_Clear(lifted);
    end;
    Clear(pp); Clear(lp);
    Deep_Clear(points);
    Clear(mixsub);
    return sols;
  end Symmetric_BKK_Solve;

  function Symmetric_BKK_Solve ( p : Poly_Sys; grp : List_of_Permutations;
                                 sign : boolean ) return Solution_List is
  
    n : constant integer32 := p'length;
    points : Array_of_Lists(p'range) := Create(p);
    sols : Solution_List;
    mix,per : Link_to_Vector;
    mixsub : Mixed_Subdivision;
    lp,pp : Laur_Sys(p'range);
    low : constant Vector := (1..n => 0);
    upp : constant Vector := Adaptive_Lifting(points);

  begin
    Compute_Mixture(points,mix,per); Clear(per);
    declare
      lifted : Array_of_Lists(mix'range);
      nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                    := (mix'range => 0.0);
      file : file_type;
    begin
      Mixed_Coherent_Subdivision
         (n,mix.all,points,false,low,upp,lifted,nbsucc,nbfail,mixsub);
      pp := Polynomial_to_Laurent_System(p);
      lp := Perform_Lifting(n,mix.all,lifted,pp);
      Create(file,out_file,"/tmp/brol");
      sols := Symmetric_Mixed_Solve(file,grp,sign,lp,mixsub,n,mix.all);
      Deep_Clear(lifted);
    end;
    Clear(pp); Clear(lp);
    Deep_Clear(points);
    Clear(mixsub);
    return sols;
  end Symmetric_BKK_Solve;

  function Symmetric_BKK_Solve
               ( file : file_type; p : Poly_Sys; sign : boolean)
               return Solution_List is

    n : constant integer32 := p'length;
    points : Array_of_Lists(p'range) := Create(p);
    sols : Solution_List;
    mix,per : Link_to_Vector;
    mixsub : Mixed_Subdivision;
    lp,pp : Laur_Sys(p'range);
    low : constant Vector := (1..n => 0);
    upp : constant Vector := Adaptive_Lifting(points);

  begin
    Compute_Mixture(points,mix,per); Clear(per);
    declare
      lifted : Array_of_Lists(mix'range);
      nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                    := (mix'range => 0.0);
    begin
      Mixed_Coherent_Subdivision
        (n,mix.all,points,false,low,upp,lifted,nbsucc,nbfail,mixsub);
      put_line(file,"****  THE SYMMETRIC MIXED SUBDIVISION  ****"); 
      put(file,natural32(n),mix.all,mixsub);
      pp := Polynomial_to_Laurent_System(p);
      lp := Perform_Lifting(n,mix.all,lifted,pp);
      new_line(file); put_line(file,"****  SYMMETRIC MIXED SOLVE  ****");
      sols := Symmetric_Mixed_Solve(file,sign,lp,mixsub,n,mix.all);
      Deep_Clear(lifted);
    end;
    Clear(pp); Clear(lp);
    Deep_Clear(points);
    Clear(mixsub);
    return sols;
  end Symmetric_BKK_Solve;

  function Symmetric_BKK_Solve ( file : file_type; p : Poly_Sys;
                                 grp :  List_of_Permutations; sign : boolean )
                               return Solution_List is

    n : constant integer32 := p'length;
    points : Array_of_Lists(p'range) := Create(p);
    sols : Solution_List;
    mix,per : Link_to_Vector;
    mixsub : Mixed_Subdivision;
    lp,pp : Laur_Sys(p'range);
    low : constant Vector := (1..n => 0);
    upp : constant Vector := Adaptive_Lifting(points);

  begin
    Compute_Mixture(points,mix,per); Clear(per);
    declare
      lifted : Array_of_Lists(mix'range);
      nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                    := (mix'range => 0.0);
    begin
      Mixed_Coherent_Subdivision
        (n,mix.all,points,false,low,upp,lifted,nbsucc,nbfail,mixsub);
      put_line(file,"****  THE SYMMETRIC MIXED SUBDIVISION  ****"); 
      put(file,natural32(n),mix.all,mixsub);
      pp := Polynomial_to_Laurent_System(p);
      lp := Perform_Lifting(n,mix.all,lifted,pp);
      new_line(file); put_line(file,"****  SYMMETRIC MIXED SOLVE  ****");
      sols := Symmetric_Mixed_Solve(file,grp,sign,lp,mixsub,n,mix.all);
      Deep_Clear(lifted);
    end;
    Clear(pp); Clear(lp);
    Deep_Clear(points);
    Clear(mixsub);
    return sols;
  end Symmetric_BKK_Solve;

end Symmetric_BKK_bound_Solvers;
