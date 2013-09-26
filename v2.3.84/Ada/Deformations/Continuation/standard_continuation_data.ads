with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Generic_Lists;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Continuation_Data is

-- DESCRIPTION :
--   In order to keep the parameters and variables manageble,
--   they have been grouped into records.
--   When path tracking we convert a solution list into an array.
--   For huge solution lists this causes indexing problems.
--   Therefore we use lists of solution arrays.

-- DATASTRUCTURES FOR VARIABLES :

  type Solu_Info is record    -- contains information about the solution

    sol : Link_to_Solution;   -- the solution: vector, t and multiplicity

    corr,cora,resr,resa,rcond : double_float; 
                              -- last correction (cor) and residual (res), 
                              -- once relative (r) and once absolute (a)
                              -- and estimate for inverse condition of jacobian

    length_path : double_float;  -- length of the path

    nstep,nfail,niter,nsyst : natural32;  -- various counters :
                              -- number of steps, failures, corrector
                              -- iterations and number of linear systems solved
  end record;

  type Solu_Info_Array is array ( integer32 range <> ) of Solu_Info;
  type Link_to_Solu_Info_Array is access Solu_Info_Array;
  package Lists_of_Solu_Info_Arrays is
    new Generic_Lists(Link_to_Solu_Info_Array);
  type Solu_Info_Array_List is new Lists_of_Solu_Info_Arrays.List;

-- CREATERS :

  function Shallow_Create ( s : Link_to_Solution ) return Solu_Info;
  function Deep_Create    ( s : Solution ) return Solu_Info;
  function Shallow_Create ( s : Solution_Array ) return Solu_Info_Array;
  function Deep_Create    ( s : Solution_Array ) return Solu_Info_Array;
  function Shallow_Create ( s : Solution_List )  return Solu_Info_Array;
  function Deep_Create    ( s : Solution_List )  return Solu_Info_Array;

  function Shallow_Create ( s : Solu_Info ) return Link_to_Solution;
  function Deep_Create    ( s : Solu_Info ) return Solution;
  function Shallow_Create ( s : Solu_Info_Array ) return Solution_Array;
  function Deep_Create    ( s : Solu_Info_Array ) return Solution_Array;
  function Shallow_Create ( s : Solu_Info_Array ) return Solution_List;
  function Deep_Create    ( s : Solu_Info_Array ) return Solution_List;

  -- DESCRIPTION :
  --   A shallow create copies the pointer to the solution, while
  --   a deep create allocates memory for a copy of the solution.

  function Create ( s : Solu_Info_Array ) return Solu_Info_Array_List;
  function Create ( s : Link_to_Solu_Info_Array ) return Solu_Info_Array_List;

  -- DESCRIPTION :
  --   Returns a list with one element: the given array s.

  function Create ( s : Solution_List; size : natural32 )
                  return Solu_Info_Array_List;

  -- DESCRIPTION :
  --   Returns a list of solution arrays where every array contains
  --   as many solutions as the given size, except perhaps for the 
  --   last array in the list which may contain fewer.

  procedure Append ( first,last : in out Solu_Info_Array_List;
                     sa : in Solu_Info_Array );

  -- DESCRIPTION :
  --   Appends the solution array to the list of solution arrays,
  --   with head in first and last element pointed to by last.

-- CONVERTOR : 

  function Concat ( s : Solu_Info_Array_List ) return Solution_List;

  -- DESCRIPTION :
  --   Concatenates the solution arrays in s into one big list.
  --   The pointers of the solutions are copied.

-- OPERATIONS ON Solu_Info :

  procedure Copy_Info ( s1 : in Solu_Info; s2 : in out Solu_Info );
  procedure Copy_Solu ( s1 : in Solu_Info; s2 : in out Solu_Info );
  procedure Copy      ( s1 : in Solu_Info; s2 : in out Solu_Info );

  -- DESCRIPTION : 
  --   Copies the information, the solution or everything from s1 to s2.

  procedure Init_Info ( s : in out Solu_Info );

  -- DESCRIPTION :
  --   Initializes the information of the solution.

  procedure Add_Info ( s1 : in out Solu_Info; s2 : in Solu_Info );

  -- DESCRIPTION :
  --   Adds the information in the counters of s2 to s1.

  procedure Update_Info ( s1 : in out Solu_Info; s2 : in Solu_Info );

  -- DESCRIPTION :
  --   Adds the information in the counters of s2 to s1 and copies the
  --   other information from s2 to s1.

-- OPERATIONS ON Solu_Info_Array :

  procedure Copy ( s : in Solu_Info_Array; sa : in out Solution_Array );
  procedure Copy ( sa : in Solution_Array; s : in out Solu_Info_Array );

  -- DESCRIPTION : Copies s into sa or vice versa.

-- DESTRUCTORS :

  procedure Clear ( s : in out Solu_Info );
  procedure Clear ( s : in out Solu_Info_Array );

  -- DESCRIPTION :
  --   This is clear is only needed after a deep create.

  procedure Deep_Clear ( s : in out Link_to_Solu_Info_Array );
  procedure Shallow_Clear ( s : in out Link_to_Solu_Info_Array );
  procedure Deep_Clear ( s : in out Solu_Info_Array_List );
  procedure Shallow_Clear ( s : in out Solu_Info_Array_List );

  -- DESCRIPTION :
  --   A shallow clear releases only the pointers, while a deep clear
  --   destroys also all content and nested pointer structures.

end Standard_Continuation_Data;
