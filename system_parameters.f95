!--------------------------------------------------------------------------------------------------
!                    NEUTRAL DRIFT AND STEM CELL FATE IN TISSUE HOMEOSTASIS
!                                    SYSTEM PARAMETERS
!
! This library contains all the parameters that define the type and characteristics of the
! algorithms of the ssa_simul7 library, allowing to easily configure any simulation.
!
! Joan Térmens Cascalló,                                                last revision: May 24, 2019
!--------------------------------------------------------------------------------------------------


module system_parameters

  logical,parameter::TestMode=.false.
! Testmode: similar to a verbose flag. Will output an additional file containing the number of
! measures as a function of the seed, which could be useful for detecting long simulations.
  logical,parameter::CompMeas=.false.
! CompMeas: Complementary measures. For each simulation, will output two files containing the density
! and position of the 1st cell, respectively, as a function of time. Warning: could produce very large files.
  logical,parameter::ShortRange=.true.
! ShortRange: determines if the cell-cell interactions are shortranged, thus simulating adhesion. The
! applied short-range is equal to cell_radius*cutoff_ratio.
  logical,parameter::GrMrg=.true.
! GrMrg: Growth-Migration. While true, the cells grow and migrate and while false, the cells only move
! due to gaussian noise without any interaction.
  integer(4),parameter::algorithm=2
! algorithm: Algorithm type, simply 1 => VM; 2 => CBD and  3 => density
  integer(4),parameter::Nini=20
! Nini: number of initial cells
  integer(4),parameter::SEED_ini=16797642
! SEED_ini: initial seed
  real(8),parameter::l0=1.d0
! l0 is the size of an adult cell.
  real(8),parameter:: L=2.d0
! L is the distance at which cells could sense density.
  integer(4),parameter::Nmeasure=10**3
! Nmeasure: number of measurements
  integer(4),parameter::max_run=10**6
! max_run: maximum number of iterations that a simulation could reach.
  integer(4),parameter::Nsimul=10**4
! Nsimul: number of simulations to be performed
  integer(4),parameter::Nread=200
! Nread: number of simulations to be done to print a screen output of # simulations and computer time.
  real(8),parameter::t_cell=200.d0
! t_cell is the average life of a cell.
  real(8),parameter::dt=0.05d0
! dt is the time integration step.
  real(8),parameter::t_growth=0.1d0*t_cell
! t_growth: time that a cell needs to become adult by growing from l0/2 to l0.
  real(8),parameter::dens0=1.d0/l0
! dens0 is the equilibrium density
  real(8),parameter::dens_ini=0.95d0
! dens_ini is the initial density
  real(8),parameter::r=0.3d0
! r is the multiplicative constant of the reaction rates used to determine the cell fate.
  real(8),parameter::K=1.d0
! K is the multiplicative constant of the quadratic interaction potential.
  real(8),parameter::delta_ini=(1.d0/dens_ini)-l0
! delta ini is the separation between cells needed to start with density = dens_ini
  real(8),parameter::D=0.02d0
! D is the difussion constant, determines the variance of the gaussian noise.
  real(8),parameter::r_c=1.5d0
! r_c is the cutoff ratio, used to determine the cutoff radius for short-ranged interactions
  real(8),parameter::Lsys=Nini*(l0+delta_ini)
! Lsys is the initial perimeter of the system.

end module system_parameters
