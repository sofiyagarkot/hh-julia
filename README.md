# Comparing the effect of different priors on simulations of Hodgkin-Huxley dynamics with probabilistic solvers

In this project my aim is to find out ***whether it is useful to customize different priors for a given ODE and how do they influence the solution***. 

The ODE of choice is the Hodgkin-Huxley (HH) neuron model. The simulation are tested for 50 ms, where the constant input current is applied during the whole period.

As a probabilistic baseline I chose an EK1 solver with fixed steps and fixed diffusion, and an IWP prior. 

In the most popular packages [NEURON — backward Euler with ∆t=0.025 ms, CSIM — Exponential Euler , GENESIS - Exponential Euler with ∆t=0.01] the HH ODE is solved with Exponential Euler method with the step size 0.01.
