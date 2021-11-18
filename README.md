# Gillespie Algorithm for Adaptive Dynamics

I'm trying to implement a general version of Gillespies Algorithm for the settup of adaptive dynamics. 
mmon implementations of Gillespies Algorithm force you to know in advance 
  a) which/how many different types there are present through the course of the simulation
  b) how the tranistion rates are between the finite number of types
  
The goal of this implementation is to get rid of these constrains. Adaptive dynamics often work in a setup where individuals are
characterized by some trait value. These traits mutate at a given rate according to some predefined mutation measure. Already in the
easy case wher an individual is characterized by a real trait say in [0,1] and the mutation measure is given by a gaussian kernel on 
the trait space common implementations of Gillespies algorithm fail to work. Since apriori there are possibly infinitely many types,
and one cannot say in advance which types will appear throughout the simulation.
Trying to approximate the trait space and the mutation measure in a discrete manner is possible, but gets computationally intens when
the approximation is good. Moreover in most of the cases mutation rates are small such that only a small number of different traits
is alive at any given time.

This implementation is ment to be versitile and fast. The biggest problem at the moment is memory allocation.

Currently only the usual implementation of Gillespies algorithm is working. Hence where one knows the number of types in advance.
You can find examples of birth and death processes (where there are no traits - only the population size is tracked) in the example folder.
