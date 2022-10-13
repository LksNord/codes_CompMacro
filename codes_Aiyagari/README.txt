
------------------------------------------------------
------------------------------------------------------
Solving the Aiyagari model with different methods         
collected for Alex Monge's Quant Macro Course @ EUI     
 
by Lukas Nord, November 2021                              
------------------------------------------------------
------------------------------------------------------

This code solves for the steady state of an Aiyagari model, using different methods to solve the household problem

The codes are meant to illustrate different solution methods as easily accessible as possible and are not necessarily optimized for maximum speed (especially VFIc and PFI). Suggestions to further improve these codes are very welcome.

Please direct errors or questions to lukas.nord@eui.eu.


----------
FILES:
----------

MAIN.m: script, sets parameters and calls functions for all solution methods

solve_aiyagari.m: function, contains GE loop over the interest rate, calls solution to HH problem and distribution

solveHH_VFId.m: function, solves HH problem by discrete value function iteration
	- iterates on Bellman Equation, discretized asset choice

solveHH_VFIc.m: function, solves HH problem by continuous value function iteration
	- iterates on Bellman Equation, continuous asset choice

solveHH_VFIcheb.m: function, solves HH problem by projection method
	- iterates on Bellman Equation approx. by Chebyshev polynomials, continuous asset choice

solveHH_PFI.m: function, solves HH problem by policy function iteration
	- iterates forward on Euler Equation, continuous asset choice

solveHH_EGM.m: function, solves HH problem by endogenous grid method
	- iterates backward on Euler Equation, continuous asset choice
	- reference: Carroll (2006, Economic Letters)

getDist_discrete.m: function, computes distribution over state-space based on discrete asset policy

getDist_continuous.m: function, computes distribution over state-space based on continuous asset policy
	- reference: Young (2010, Journal of Economic Dynamics and Control)

interpLN.m: function, simple interpolation routine

cheby_vec.m: function, computes for given x values of the Chebyshev polynomials up to order N

markovappr.m: function, discretizes AR(1) process




