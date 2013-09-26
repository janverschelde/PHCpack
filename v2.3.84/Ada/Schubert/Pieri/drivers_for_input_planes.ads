with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Brackets;                          use Brackets;
with Standard_Floating_Vectors;         use Standard_Floating_Vectors;
with Standard_Complex_VecMats;          use Standard_Complex_VecMats;

package Drivers_for_Input_Planes is

-- DESCRIPTION :
--   This package provides menu's to generate random and real input planes
--   for the three types of Pieri homotopies.
--   The open design of the package is due to its use in both the usual
--   as in the quantum Pieri homotopies.

-- READING CO-DIMENSION CONDITIONS :

  function Read_Codimensions ( m,p,q : natural32 ) return Bracket;

  -- DESCRIPTION :
  --   Reads the vector of co-dimensions with a test whether these
  --   co-dimensions condition will lead to a finite number of solutions.

-- PRIMITIVES :

  function Random_Complex_Planes ( m,p : natural32 ) return VecMat;
  function Random_Real_Planes    ( m,p : natural32 ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of range 1..m*p with random complex/real m-planes,
  --   with the columns orthonormalized.

  function Random_Complex_Planes
             ( m,p : natural32; k : Bracket ) return VecMat;
  function Random_Real_Planes 
             ( m,p : natural32; k : Bracket ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of the same range as k with random (m+1-k(i))-planes,
  --   with the columns orthonormalized.

  function Random_Complex_Planes ( m,p,q : natural32 ) return VecMat;
  function Random_Real_Planes    ( m,p,q : natural32 ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of range 1..m*p+q*(m+p) with random complex or real
  --   m-planes, --   with the columns orthonormalized.

  function Equidistant_Interpolation_Points ( n : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Generates n equidistant interpolation points in [-1,+1],
  --   starting at a random value.

  function Read_Interpolation_Points ( n : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Reads n s-values from standard input.

  function Osculating_Input_Planes ( m,p : natural32 ) return VecMat;
  function Osculating_Input_Planes ( m,p : natural32;
                                     s : Vector ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of range 1..m*p with real m-planes osculating
  --   a rational normal curve, sampled at equidistant points in [-1,+1],
  --   or, if specified at the given s-values.

  function Complex_Osculating_Input_Planes 
                 ( m,p : natural32 ) return VecMat;
  function Complex_Osculating_Input_Planes 
                 ( m,p : natural32; s : Vector ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of range 1..m*p with complex m-planes osculating
  --   a rational normal curve, with s values sampled at equidistant points
  --   on the imaginary axis between -i and i, or if specified, at the
  --   points on the imaginary axis in s.

  function Osculating_Input_Planes
             ( m,p : natural32; k : Bracket ) return VecMat;
  function Osculating_Input_Planes
             ( m,p : natural32; k : Bracket; s : Vector ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of the same range as k with real (m+1-k(i))-planes
  --   osculating a rational normal curve, sampled at equidistant points
  --   in [-1,+1], or, if specified at the given s-values.

  function Osculating_Input_Planes ( m,p,q : natural32 ) return VecMat;
  function Osculating_Input_Planes
             ( m,p,q : natural32; s : Vector ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of the range 1..m*p+q*(m+p) with real m-planes
  --   osculating a rational normal curve, sampled at equidistant points
  --   in [-1,+1], or, if specified at the given s-values.

  function Read_Input_Planes ( m,p : natural32 ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of range 1..m*p with real m-planes retrieved from
  --   a file given by the user after calling this function.
  --   The matrices on return have their columns orthonormalized.

  function Read_Input_Planes ( m,p : natural32; k : Bracket ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of the same range as k with real (m+1-k(i))-planes 
  --   retrieved from a file given by the user after calling this function.
  --   The matrices on return have their columns orthonormalized.

  function Read_Input_Planes ( m,p,q : natural32 ) return VecMat;

  -- DESCRIPTION :
  --   Returns a vector of range 1..m*p+q*(m+p) with real m-planes retrieved
  --   from a file given by the user after calling this function.
  --   The matrices on return have their columns orthonormalized.

-- MAIN INTERACTIVE DRIVERS :

  procedure Driver_for_Input_Planes
              ( file : in file_type; m,p : in natural32; planes : out VecMat );

  -- DESCRIPTION :
  --   Generates m-planes as input to the hypersurface Pieri algorithm.

  -- ON ENTRY :
  --   file     to write output logistics on;
  --   m        dimension of the input planes, number of columns;
  --   p        dimension of the output planes, m+p = number of rows.

  -- ON RETURN :
  --   planes   vector of range 1..m*p with m-planes in dimension m+p.

  procedure Driver_for_Input_Planes
              ( m,p : in natural32; k : in Bracket; planes : out VecMat );
  procedure Driver_for_Input_Planes
              ( file : in file_type;
                m,p : in natural32; k : in Bracket; planes : out VecMat );

  -- DESCRIPTION :
  --   Generates m-planes as input to the general Pieri algorithm.

  -- ON ENTRY :
  --   file     to write coordinates of target input planes on;
  --   m        number of columns of the input planes is m+1-k(i);
  --   p        dimension of the output planes, m+p = number of rows;
  --   k        co-dimension conditions of the input planes.

  -- ON RETURN :
  --   planes   vector of same range as k with (m+1-k(i))-planes in
  --            a space of dimension m+p.

  procedure Driver_for_Input_Planes
              ( file : in file_type; m,p,q : in natural32;
                s : out Vector; planes : out VecMat );

  -- DESCRIPTION :
  --   Generates m-planes as input to the quantum Pieri algorithm.

  -- ON ENTRY :
  --   file     to write output logistics on;
  --   m        number of columns of the input planes is m+1-k(i);
  --   p        dimension of the output planes, m+p = number of rows;
  --   q        degree of the output maps.

  -- ON RETURN :
  --   s        vector of range 1..m*p+q*(m+p) of interpolation points;
  --   planes   vector of range 1..m*p+q*(m+p) with m-planes in
  --            a space of dimension m+p, sampled at s-values.

end Drivers_for_Input_Planes;
