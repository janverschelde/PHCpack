with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Generic_Lists;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Permutations,Symmetry_Group;        use Permutations,Symmetry_Group;

package Orbits_of_Solutions is

-- DESCRIPTION :
--   This package contains routines for manipulating the
--   solutions set of a symmetric polynomial system.

-- DATA STRUCTURES :

  type Orbit ( n : integer32 ) is record
    orb : Permutation(1..n);  -- data representing the orbit
     -- e.g.: orb = [ 1 1 2 2 ] => solution is like (a,a,b,b) or (a,b,a,b)
    nbdiff : natural32;         -- number of different values
    nbgen  : natural32;         -- number of generating solutions
    nbsols : natural32;         -- number of solutions in the orbit
  end record;

  type Link_to_Orbit is access Orbit;

  package Lists_of_Orbits is new Generic_Lists(Link_to_Orbit);
  type List_of_Orbits is new Lists_of_Orbits.List;

-- CONSTRUCTORS :

  function Generating ( sols : Solution_List; sign : boolean;
                        tol : double_float ) return Solution_List;

  -- DESCRIPTION :
  --   Returns a list of generating solutions, by permutations, all
  --   other solutions in sols can be derived from the resulting list.
  --   If sign is true, then also permutations who alter the sign of
  --   the components are tried.

  procedure Analyze ( L : in List_of_Permutations; sign : in boolean;
                      tol : in double_float; sols : in out Solution_List );

  -- DESCRIPTION :
  --   the solution list sols will be checked upon symmetry, according
  --   to the list of permutations l

  -- ON ENTRY :
  --   l              a list of permutations;
  --   sign           if true then sign symmetry has to be taken into account;
  --   tol            tolerance to decide wether two solutions are equal;
  --   sols           a list of solutions.

  -- ON RETURN :
  --   sols           a list of solutions, where only one element per orbit
  --                  is in the list; 
  --                  the multiplicity field is used to indicate the number
  --                  of elements in the orbit.

  procedure Orbit_Structure 
                    ( s : in Solution; tol : in double_float;
                      orbit : in out Permutation; nbdiff : out natural32 );

  -- DESCRIPTION :
  --   This procedure returns the structure of the orbit of the solution.

  -- ON ENTRY :
  --   s              a solution;
  --   tol            tolerance to decide wether two solutions are equal.

  -- ON RETURN :
  --   orbit          orbit(i) = j means that x_i has the same value as x_j;
  --   nbdiff         the number of different values in the orbit.

  function Orbits ( sols : Solution_List; tol : double_float )
                  return Permutation;

  -- DESCRIPTION :
  --   Let orb := Orbits(sols,tol), 
  --   then orb(i) gives the number of solutions belonging to an orbit
  --   with i different values.

  -- ON ENTRY :
  --   sols         a list of solutions;
  --   tol          tolerance to decide wether two solutions are equal.

  procedure Orbits ( grp : in List_of_Permutations; tol : in double_float;
                     sols : in out Solution_List;
                     lorb : in out List_of_Orbits );

  -- DESCRIPTION :
  --   The list of solutions will be analyzed according to the
  --   symmetry group.  Together with the information on the orbits,
  --   the list of generating solutions will be produced.

  -- ON ENTRY :
  --   grp           a list of permutations representing the symmetry group;
  --   tol           used to compare two solutions;
  --   sols          the list of solutions.

  -- ON RETURN :
  --   sols          the generating list of solutions;
  --   lorb          the list containing the information of the orbits.

-- SELECTOR :

  function Same_Orbit ( orb1,orb2 : Permutation ) return boolean;

  -- DESCRIPTION :
  --   Returns true when the structures of the orbits are the same.
  --   For example : 1 1 1 2 = 1 2 1 1, but 1 1 2 2 /= 2 2 1 2.

-- DESTRUCTOR :

  procedure Clear ( lorb : in out List_of_Orbits );

  -- DESCRIPTION :
  --   All memory space allocated for lorb will be freed.
		    
end Orbits_of_Solutions;
