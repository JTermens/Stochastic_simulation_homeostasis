!--------------------------------------------------------------------------------------------------
!                    NEUTRAL DRIFT AND STEM CELL FATE IN TISSUE HOMEOSTASIS
!                                    SYSTEM PARAMETERS
!
! This library contains all the parameters that define the type and characteristics of the
! algorithms of the ssa_simul7 library, allowing to easely configurate simulations.
!
! Joan Térmens Cascalló,                                                last revision: May 24, 2019
!--------------------------------------------------------------------------------------------------


module system_parameters

  logical,parameter::TestMode=.false.
  logical,parameter::CompMeas=.false.
  logical,parameter::ShortRange=.true.
  logical,parameter::GrMrg=.true.
  integer(4),parameter::algorithm=2 !1=>VM; 2=>CBD; 3=>density
  integer(4),parameter::Nini=20
  integer(4),parameter::SEED_ini=16797642
  real(8),parameter::l0=1.d0
  real(8), parameter::L=2.d0
  integer(4),parameter::Nmeasure=10**3
  integer(4),parameter::max_run=10**6
  integer(4),parameter::Nsimul=10**4
  integer(4),parameter::Nread=200
  real(8),parameter::t_cell=200.d0,dt=0.05d0
  real(8),parameter::t_growth=0.1d0*t_cell
  real(8),parameter::dens0=1.d0/l0,dens_ini=0.95d0,r=0.3d0,K=1.d0
  real(8),parameter::delta_ini=(1.d0/dens_ini)-l0,D=0.02d0,r_c=1.5d0
  real(8),parameter::Lsys=Nini*(l0+delta_ini)

end module system_parameters
