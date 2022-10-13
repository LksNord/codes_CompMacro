------------------------------------------------------
------------------------------------------------------
Solving Transitions of the Aiyagari model              
collected for Alex Monge's Quant Macro Course @ EUI    
   
by Lukas Nord, December 2021                                                                                         
------------------------------------------------------
------------------------------------------------------

These codes solve for the dynamics of an Aiyagari model in response to a 
one-off (MIT) productivity shock. The household problem is solved by
endogenous grid method. For the aggregate dynamics we apply:
a) extended path algortihm (as previously)
b) updates of the interest path using the sequence space Jacobian
   (see Auclert et al. (ECMA,2021))

The codes are meant to illustrate solution methods as easily
accessible as possible and are not necessarily optimized for maximum
speed. Suggestions to further improve the codes are very welcome.

Please direct errors or questions to lukas.nord@eui.eu.


----------
FILES:
----------

MAIN.m: script, sets parameters and calls functions to solve transitions

solve_aiyagari.m: function, contains GE loop over the steady state interest rate, calls solution to HH problem and distribution

solveHH_EGM.m: function, solves HH problem by endogenous grid method
	- iterates backward on Euler Equation, continuous asset choice
	- reference: Carroll (2006, Economic Letters)

getDist_continuous.m: function, computes distribution over state-space based on continuous asset policy
	- reference: Young (2010, Journal of Economic Dynamics and Control)

solve_trans.m: function, computes transitional dynamics to an MIT productivity shock
	- updating of interest rate path based on extended path approach

solve_transSSJ.m: function, computes transitional dynamics to an MIT productivity shock
	- updating of interest rate path based on sequence space Jacobian approach

compute_Jac.m: function, computes the sequence space Jacobian of capital market clearing wrt the path of interest rates

stepEGM.m: function, iterates one step on the household problem by endogenous grid method

stepDist.m: function, iterates one step on the distribution across HHs by Young method

markovappr.m: function, discretizes AR(1) process




