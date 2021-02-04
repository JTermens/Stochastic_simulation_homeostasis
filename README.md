# NEUTRAL DRIFT AND STEM CELL FATE IN TISSUE HOMEOSTASIS
## STOCHASTIC SIMULATION ALGORITHMS (version 7)
 
**Author:** Joan Térmens

This repository contains the main algoithms to simulate the collective behaviour of a 1-dim tissue of stem cells with periodic boundary conditions. Each cell choses between proliferating or differentiating (leaving the tissue) according a different model (VM,CBD or density) so that homeostasis is achieved and maintained. Each simulation begins with a single labeled cell and tracks its descendants until the clone dies, reaches the tissue's size or the last measure time is reached. The dynamics are simulated by a simple model of linear growth/migration in which each cell grows util the adult size (l0) and migrates to be at the equilibrium distance with its neighbouring cells. Each algorith simulates a different cell fate model of the following:

* **Critical birth-death model (CBD):** each cell choses between differentiation or proliferation stochastically with *prob=1/2* for both. Lifetimes are generated stochastically as *\~ EXP(1/t_cell)*.

* **Voter Model (VM):** each cell stochastically choses one side and evolves with the opposite fate as that side's neighbour cell. Lifetimes are generated stochastically as *\~ EXP(1/t_cell)*.

* **Density model:** cell fate is determined stochastically in function of the local density (at *\[-L,+L\]* around the progenitor). With this density, two stochastic lifetimes are generated, following *t1 \~ EXP(1/w1\*t_cell)* and *t2 \~ EXP(1/w2\*t_cell)*. If *t1 < t2* the cell will proliferate after *t1*, if not, the cell differentiates after *t2*. This model can be simulated with long and short-range (adhesion) interactions between cells, as well as without any interaction.

The code is written in `Fortran 95` and organized into three modules: `ssa_simul7.f95` which contains the simulation algorithms, `mtmod.f95` which contains the algorithms responsible of random number generation and `system_parameters.f95` which contains the parameters that characterize the simulations. To facilitate its use, the program has been thought to be configured by only changing the `system_parameters.f95` file.

The simulations could be easily compiled with any `Fortran 95` compiler. For instance, a simlple way using `gfortran` could be:

```
gfortran -c system_parameters.f95 mtmod.f95 ssa_simul7.f95
gfortran -o executable_name.out main_simul7.f95 system_parameters.o mtmod.o ssa_simul7.o
./executable_name.out
```
**WARNING** Beware with the number of simulations performed (`Nsimul`) as a low number would lead to bad results due to the lack of statistically significant sample, whereas a big number produces intensive executions that could last for a few hours.

To know more about the scientific background and application of this work, check [this paper](http://diposit.ub.edu/dspace/handle/2445/141704).
