# Comparing the effect of different priors on simulations of Hodgkin-Huxley dynamics with probabilistic solvers

In this project my aim is to find out ***whether it is useful to customize different priors for a given ODE and how do they influence the solution***. 

The ODE of choice is the Hodgkin-Huxley (HH) neuron model. The simulation are tested for 50 ms, where the constant input current is applied during the whole period.

As a probabilistic baseline I chose an EK1 solver with fixed steps and fixed diffusion, and an IWP prior. 

In the most popular packages [NEURON — backward Euler with ∆t=0.025 ms, CSIM — Exponential Euler ? , GENESIS - Exponential Euler with ∆t=0.01] the HH ODE is solved with Exponential Euler method with the step size 0.01.

- [ ] EK1 different orders — WP diagram with different number of evaluations, also including RK(…)
- [ ] Mean absolute error during the simulation time for EK1( IWP(1…4) )
- [ ] Maybe? include comparison with EK0
- [ ] Also maybe? WP diagram with time vs error but that one will show more kind of efficacy of implementation of EK1

## 1. The first aspect that I’ll study is how choosing a different prior (Matern or IOUP) effects the solution (putting the same prior on all the dimensions of an ODE).
! the num_derivatives = 3, ∆t = 0.01, FixedDiffusion model throughout all the experiments

  - ### Integrated Ornstein-Uhlenbeck process (IOUP)
      - [] describe paramterization of prior 
      - [] one visual on how does this look like
 
    How does total MAE depends on the choice of parameters
      - [] plot showing for IOUP(num_derivatives = 3) the effect of rate parameter on total MAE; as a red dashed line put the performance of the IWP(3);
      - [] given the optimal parameters — plot the prior and pick a suboptimal parametrization and plot priors for them too; initializing with the stationary distribution ! ; discuss the effect

    How do the parameters influence the speed/accuracy
    - [] given the optimal and suboptimal parameters, compare on WP diagrams the performance error vs time (or time vs error) relative to each other and IWP(3)   

 
  - ### Matern prior
      - [] decribe parametrization
      - [] one visual on how does this prior look like

    How does total MAE depends on the choice of parameters
      - [] plot showing the effect of lengthscale on the same total MAE; as a red dashed line put the performance of the IWP(3)
      - [] same suboptimal vs. optimal plots (in the terms of resulting in smaller MAE)
        
    How do the parameters influence the speed/accuracy
    - [] given the optimal and suboptimal parameters, compare on WP diagrams the performance error vs time (or time vs error) relative to each other and IWP(3)   

  - ### Comparing the Matern, IOUP and IWP with optimal parameters between each other
    - [] WP diagram
    - [] comparing the mean absolute errors across the whole simulation

## 2. How does combination of priors influences the solution? 

## 3. How does transforming the priors before entering the ODE influences the solution? 


[NEURON] — Carnevale, N.T. and Hines, M.L. (2006) The Neuron Book. Cambridge University Press, Cambridge, UK. http://dx.doi.org/10.1017/CBO9780511541612

[CSIM] — … 

[GENESIS] — http://genesis-sim.org/GENESIS/iBoG/iBoGpdf/chapt2.pdf
