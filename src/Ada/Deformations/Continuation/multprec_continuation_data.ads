with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;
with Continuation_Parameters;
with Standard_Continuation_Data;

package Multprec_Continuation_Data is

-- DESCRIPTION :
--   In order to keep the parameters and variables manageble,
--   they have been grouped into records.

-- DATA STRUCTURES FOR PARAMETERS :

  type Pred_Pars is record    -- contains the parameters for the predictor

    minstep,maxstep : Floating_Number;  -- minimum and maximum step size
    expfac,redfac : Floating_Number;    -- expansion and reduction factor
                                        --  for step length control
    success_steps : natural32;   -- number of successful steps before expansion
    predictor_type : natural32;  -- type of predictor used
    dist_target : Floating_Number;  -- distance to target
    power : natural32;           -- power of t in (polyhedral) homotopy

  end record;

  type Corr_Pars is record    -- contains the parameters for the corrector

    epsrx,epsax,epsrf,epsaf : Floating_Number;  
                              -- desired precisions for x and its residual f(x)
                              -- once relative (r) and once absolute (a)

    maxit,maxtot : natural32; -- maximum number of corrector iterations
                              -- for one step and for the whole path
  end record;

-- DATASTRUCTURES FOR VARIABLES :

  type Solu_Info is record    -- contains information about the solution

    sol : Link_to_Solution;   -- the solution: vector, t and multiplicity

    corr,cora,resr,resa,rcond : Floating_Number; 
                              -- last correction (cor) and residual (res), 
                              -- once relative (r) and once absolute (a)
                              -- and estimate for inverse condition of jacobian

    length_path : Floating_Number;  -- length of the path

    nstep,nfail,niter,nsyst : natural32;  -- various counters :
                              -- number of steps, failures, corrector
                              -- iterations and number of linear systems solved
  end record;

  type Solu_Info_Array is array ( integer32 range <> ) of Solu_Info;

-- CONVERTORS :

  function Convert ( p : Continuation_Parameters.Pred_Pars )
                   return Multprec_Continuation_Data.Pred_Pars;
 
  -- DESCRIPTION :
  --   Converts a data structure with predictor parameters from
  --   standard to multi-precision arithmetic.

  function Convert ( c : Continuation_Parameters.Corr_Pars )
                   return Multprec_Continuation_Data.Corr_Pars;

  -- DESCRIPTION :
  --   Converts a data structure with corrector parameters from
  --   standard to multi-precision arithmetic.

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

end Multprec_Continuation_Data;
