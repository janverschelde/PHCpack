with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Span_of_Component;                 use Span_of_Component;

package Span_of_Component_Creators is

-- DESCRIPTION :
--   To create the span of a component, samples are taken until the
--   span is either found to be the entire space or some subspace.

  procedure Create_Span
               ( spt : in Standard_Sample; tol : in double_float;
                 sps,sps_last : in out Standard_Sample_List;
                 sp : out Standard_Span );
  procedure Create_Span
               ( spt : in Standard_Sample; tol : in double_float;
                 size : in natural32;
                 sps,sps_last : in out Multprec_Sample_List;
                 sp : out Multprec_Span );

  -- DESCRIPTION :
  --   Takes enough samples until the span is determined.

  -- REQUIRED :
  --   The sampling machine is initialized and properly tuned.

  -- ON ENTRY :
  --   spt       a sample on the component;
  --   tol       tolerance to decide whether a number is zero;
  --   size      size of the multi-precision numbers (optional).

  -- ON RETURN :
  --   sps       samples used to create the span;
  --   sps_last  pointer to last element in the list sps;
  --   sp        span of component.

end Span_of_Component_Creators;
