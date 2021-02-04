!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------

PROGRAM main_simul7
  USE ssa_simul7
  USE system_parameters
  IMPLICIT NONE

  integer(4)::i,j,i_simul
  real(4)::TIME1,TIME2

  integer(4)::SEED,k_measure
  real(8)::time_measure(Nmeasure+1)
  integer(4)::sum_Nsurv(Nmeasure+1),num_sumands(Nmeasure+1),sum_Ncell(Nmeasure+1)
  real(8)::sum_frac_diff(Nmeasure+1)
  real(8)::average_Nsurv(Nmeasure+1),average_frac(Nmeasure+1)
  integer(4)::clone_size_times(7),clone_size_Nsim(7)
  integer(4)::clone_size_sum(2*Nini,7)
  real(8)::N(7),prob(7)

  integer(4)::Nsurv_rec(Nmeasure+1),Ncell_rec(Nmeasure+1)
  real(8)::frac_rec(Nmeasure+1)

  logical::err
  integer(4)::index_err

  character(2)::charact_L

  write(charact_L,fmt="(I2)") int(L)

  if (algorithm.eq.1) then
    print*,"VOTER MODEL"
  else if (algorithm.eq.2) then
    print*,"CRITICAL BIRTH-DEATH MODEL"
  else if (algorithm.eq.3) then
    print*, "DENSITY, L=",L,", GrMrg=",GrMrg,", ShortRange=",ShortRange
  endif

  SEED=SEED_ini

  time_measure(1)=0.d0
  do i=1,Nmeasure
    time_measure(i+1)=i*t_cell
  end do

  clone_size_times=[5,25,50,100,250,500,750]

  index_err=0
  err=.false.

  sum_Nsurv=0
  sum_frac_diff=0.d0
  sum_Ncell=0
  num_sumands=0
  clone_size_sum=0
  clone_size_Nsim=0

  call cpu_time(TIME1)

  if (TestMode) then
    open(unit=14, file="Test_Data.dat")
  endif

  do i_simul=1,Nsimul
    SEED=SEED+1
    if (algorithm.eq.1) then
      call main_VM(SEED,time_measure,Nsurv_rec,k_measure)
    else if (algorithm.eq.2) then
      call main_CBD(SEED,time_measure,Nsurv_rec,frac_rec,k_measure)
    else if (algorithm.eq.3) then
      call main_dens(SEED,time_measure,Ncell_rec,Nsurv_rec,frac_rec,k_measure,err)
    endif

    if (TestMode) then
      write(14,*) "seed=",SEED," k_measure=",k_measure," err=",err
    endif

    if(err.eqv..true.) then
      index_err=index_err+1
    endif

    do j=1,k_measure
      sum_Nsurv(j)=sum_Nsurv(j)+Nsurv_rec(j)
      sum_frac_diff(j)=sum_frac_diff(j)+frac_rec(j)
      sum_Ncell(j)=sum_Ncell(j)+Ncell_rec(j)
      num_sumands(j)=num_sumands(j)+1
    enddo

    do i=1,7
      if(k_measure.gt.(clone_size_times(i)+1)) then
        clone_size_Nsim(i)=clone_size_Nsim(i)+1
        do j=1,Nsurv_rec(clone_size_times(i)+1)
          clone_size_sum(j,i)=clone_size_sum(j,i)+1
        enddo
      end if
    end do

    if(Nread*(i_simul/Nread).eq.i_simul) then !Seguiment per pantalla del
      call cpu_time(TIME2)          !progés del programa
      write(*,*) "# simulations = ",i_simul," CPU time = ",TIME2-TIME1
    endif
  end do

  if (TestMode) then
    write(14,*) " "
    write(14,*) "# simulations with errors/warnings=",index_err
    close(14)
  endif

  if (algorithm.eq.1) then
    open(unit=12,file="VM_average.dat")
    open(unit=13,file="VM_clone_size.dat")

  else if (algorithm.eq.2) then
    open(unit=12,file="CBD_average.dat")
    open(unit=13,file="CBD_clone_size.dat")

  else if (algorithm.eq.3) then
    if (GrMrg) then
      if (ShortRange) then
        open(unit=12,file="L="//charact_L//"_average_GM_SR.dat")
        open(unit=13,file="L="//charact_L//"_clone_size_GM_SR.dat")
      else
        open(unit=12,file="L="//charact_L//"_average_GM_LR.dat")
        open(unit=13,file="L="//charact_L//"_clone_size_GM_LR.dat")
      endif
    else
      open(unit=12,file="L="//charact_L//"_average.dat")
      open(unit=13,file="L="//charact_L//"_clone_size.dat")
    endif
  endif

  do i=1,(Nmeasure+1)
    average_Nsurv(i)=real(sum_Nsurv(i))/real(num_sumands(i))
    if (algorithm.ne.1) then
      average_frac(i)=sum_frac_diff(i)/real(num_sumands(i))
      write(12,*) time_measure(i),average_Nsurv(i),average_frac(i),num_sumands(i)
    else
      write(12,*) time_measure(i),average_Nsurv(i),num_sumands(i)
    endif
  enddo

  do i=1,2*Nini
    do j=1,7
      N(j)=real(i)/average_Nsurv(clone_size_times(j)+1)
      prob(j)=real(clone_size_sum(i,j))/real(clone_size_Nsim(j))
    enddo
    write(13,*) N(1),prob(1),N(2),prob(2),N(3),prob(3),N(4),prob(4),N(5),prob(5),N(6),prob(6),N(7),prob(7)
  enddo

  close(12)
  close(13)

END PROGRAM main_simul7
