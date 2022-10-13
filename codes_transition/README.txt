------------------------------------------------------
------------------------------------------------------
Solving Transitions of the Aiyagari model              
collected for Alex Monge's Quant Macro Course @ EUI    
   
by Lukas Nord, November 2021                                                                                         
------------------------------------------------------
------------------------------------------------------

These codes solve for transitions between two steady states of a baseline
Aiyagari model. The code can produce either a deterministic transition
to a new steady state after a change in the labor tax rate or the IRF of
the economy in response to a one-off (MIT) productivity shock.

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

solve_trans.m: function, computes transitional dynamics between steady states
	- can solve for either transition after changes in labor income tax or for dynamics after transitory (MIT) productivity shock

stepEGM.m: function, iterates one step on the household problem by endogenous grid method

stepDist.m: function, iterates one step on the distribution across HHs by Young method

markovappr.m: function, discretizes AR(1) process




