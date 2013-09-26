procedure Sensitivity_of_Factorization;

-- DESCRIPTION :
--   Investigation of the numerical sensitivity of the factorization of
--   polynomials in two variables with approximate complex coefficients.
--   The user is prompted for a number of factors and for the degrees of
--   each factor.  The monodromy breakup algorithm is launched on 
--   several perturbed instances on the product of the factors.
--   The perturbations are of magnitude 10^(-i), with i from 0 to 14.
