with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_c2fac ( job : integer32;
                     a : C_intarrs.Pointer;
		     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32;

-- DESCRIPTION :
--   Provides a gateway from C to the operations in PHCpack
--   to factor pure dimensional solution sets of polynomial systems,
--   using monodromy certified by linear traces.

-- ON ENTRY :
--   job    =  0 : display the menu of all available options;
--   job    =  1 : prompts the user for a witness set,
--                 stores the system in the systems container,
--                 and its solutions in the solutions container,
--                 and returns in a the dimension of the ambient space
--                 and in b the following two numbers:
--                   b[0] : dimension of the solution set;
--                   b[1] : degree of the solution set.
--   job    =  2 : takes the system and the solutions from the containers
--                 and initializes the sampling machine, on input,
--                 the dimension of the witness set must be in a;
--   job    =  3 : assigns the coefficient c[0] + c[1]*I to the
--                 b-th coefficient of the a-th slice;
--   job    =  4 : stores the c[0] + c[1]*I as random gamma constant
--                 for the a-th equation;
--   job    =  5 : compute a new witness set on the new slices;
--   job    =  6 : swaps slices and solution sets to turn back;
--   job    =  7 : copy embedded system from sampler to systems container;
--   job    =  8 : copy first solution list to container;
--   job    =  9 : put solutions with index in a from monodromy grid
--                 in the solutions container (companion to job = 11);
--   job    = 10 : initializes Monodromy_Permutations with two numbers:
--                   a[0] : number of monodromy loops,
--                   b[0] : degree of the solution component to factor;
--   job    = 11 : store solutions in container to Monodromy_Permutations;
--   job    = 12 : compute permutation by last stored solution list,
--                 and return this new permutation in b;
--   job    = 13 : updates decomposition with a new permutation,
--                 a[0] must contain the dimension and b the permutation;
--   job    = 14 : writes the current decomposition;
--   job    = 15 : applies the linear trace to certify the decomposition;
--   job    = 16 : returns in c the diagnostics of the trace grid;
--   job    = 17 : returns in c difference between trace and actual sum.
--   job    = 18 : finds the index of a solution label in a slice,
--                 on entry: a[0] is label to a solution,
--                           a[1] is the number of a slice;
--                 on return: b is index to solution if label occurs,
--                            otherwise, b is zero;
--   job    = 19 : initialize number of slices in Sampling_Operations
--                 with the content of a;
--   job    = 20 : adds a new slice to Sampling_Operations, 
--                 where a[0] = total number of doubles in the slices;
--                       a[1] = dimension of the solution set;
--                       a[2] = the ambient dimension;
--                 the coefficients are in c.
--   job    = 21 : returns in c the coefficients of a slice,
--                 where a[0] = total number of doubles in the slices;
--                       a[1] = dimension of the solution set;
--                       a[2] = the ambient dimension;
--                 the index to the slice is in b;
--   job    = 22 : sets the target slices to the a-th slice stored 
--                 in Sampling_Operations.
--   job    = 23 : completes one loop, sampling from one solution,
--                 where a[0] = index for the starting hyperplane sections;
--                       a[1] = index for the target hyperplanes sections;
--                   and b = label of the start solution;
--                 on return b contains the label of the matching solution
--                 in the list at the target hyperplane sections.
--   job    = 24 : reads a witness set from file,
--                 on input in b is the file name, and in a the
--                          number of characters in the file name,
--                 stores the system in the systems container,
--                 and its solutions in the solutions container,
--                 and returns in a the dimension of the ambient space
--                 and in b the following two numbers:
--                   b[0] : dimension of the solution set;
--                   b[1] : degree of the solution set.
--   job    = 25 : reads a witness set to file,
--                 on input in b is the file name, and in a the
--                          number of characters in the file name,
--                 the system and the solutions are taken from containers.
--   job    = 26 : returns in a the number of irreducible factors in
--                 the current irreducible decomposition.
--   job    = 27 : given in a an index k to an irreducible component,
--                 returns in a the degree of the k-th component and
--                 in b the labels of the points that span the k-th
--                 component in the current irreducible decomposition.
--
-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the right range.
