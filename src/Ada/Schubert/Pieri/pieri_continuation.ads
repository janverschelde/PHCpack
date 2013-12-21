with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with Localization_Posets;                use Localization_Posets;

package Pieri_Continuation is

-- DESCRIPTION :
--   Invokes the polynomial continuation on the Pieri homotopies.

  procedure Trace_Paths
               ( file : in file_type; homsys : in Poly_Sys;
                 locmap : in Standard_Natural_Matrices.Matrix;
                 report,outlog : in boolean;
                 plane : in out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Traces the paths defined by the Pieri homtopy for one given plane.

  -- ON ENTRY :
  --   file      to write intermediate results on;
  --   homsys    Pieri homotopy, where t is the last variable;
  --   locmap    localization pattern, 0,1 in identity, 2 = free;
  --   report    switch to indicate whether reporting path tracker;
  --   outlog    if switched on, writes homotopies;
  --   plane     solution plane at t = 0.

  -- ON RETURN :
  --   plane     solution plane at t = 1.

  procedure Trace_Paths
               ( file : in file_type; homsys : in Poly_Sys;
                 locmap : in Standard_Natural_Matrices.Matrix;
                 report,outlog : in boolean;
                 planes : in Standard_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Traces the paths defined by the Pieri homotopy to all planes.

  -- ON ENTRY :
  --   file      to write intermediate results on;
  --   homsys    Pieri homotopy, with t as last variable;
  --   locmap    localization pattern, 0,1 in identity, 2 = free;
  --   report    switch to indicate whether reporting path tracker;
  --   outlog    if switched on, writes homotopies;
  --   planes    solution planes at t = 0.

  -- ON RETURN :
  --   planes    solution planes at t = 1.

  procedure Quantum_Trace_Paths
               ( m,p,q : in natural32; nd : in Node;
                 homsys : in Poly_Sys; conpar,s_mode : in natural32;
                 locmap : in Standard_Natural_Matrices.Matrix;
                 planes : in Standard_Complex_VecMats.VecMat );
  procedure Quantum_Trace_Paths
               ( file : in file_type;
                 m,p,q : in natural32; nd : in Node;
                 xpm : in Standard_Complex_Poly_Matrices.Matrix;
                 s : in Standard_Complex_Vectors.Vector; ip : in VecMat;
                 homsys : in Poly_Sys; conpar,s_mode : in natural32;
                 locmap : in Standard_Natural_Matrices.Matrix;
                 report,outlog : in boolean;
                 planes : in Standard_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   The q-analogue of the path tracing in the Pieri homotopy algorithm.

  -- ON ENTRY :
  --   file      to write intermediate results on;
  --   m         dimension of the input planes;
  --   p         dimension of the output planes;
  --   q         degree of the solution maps;
  --   nd        current node;
  --   xpm       symbolic representation of the solution maps;
  --   s         interpolation points;
  --   ip        interpolation conditions;
  --   homsys    Pieri homotopy, with t as last variable;
  --   conpar    index of the continuation parameter;
  --   s_mode    determines starting value of s,
  --             if equal to zero, then s is zero at the start,
  --             otherwise s is one at the start;
  --   locmap    localization pattern, 0,1 in identity, 2 = free;
  --   report    switch to indicate whether reporting path tracker;
  --   outlog    if switched on, writes homotopies;
  --   planes    solution planes at t = 0.

  -- ON RETURN :
  --   cnt       updated counter;
  --   planes    solution planes at t = 1.

end Pieri_Continuation;
