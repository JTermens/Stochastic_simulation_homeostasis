!--------------------------------------------------------------------------------------------------
!                    NEUTRAL DRIFT AND STEM CELL FATE IN TISSUE HOMEOSTASIS
!                      STOCHASTIC SIMULATION ALGORITHMS, ssa, (version 7)
!
! This library contains the main algoithms to simulate the collective behaviour of a 1-dim tissue
! of stem cells with periodic boundary conditions. Each cell choses between proliferating or
! differentiating (leaving the tissue) according a different model (VM,CBD or density) so that
! homeostasis is achieved and maintained. Each simulation begins with a single labeled cell and
! tracks its descendants until the clone dies, reaches the tissue's size or the last measure time
! is reached.
! The dynamics are simulated by a simple model of linear growth/migration in which each cell grows
! util the adult size (l0) and migrates to be at the equilibrium distance with its neighbouring
! cells. Each algorith simulates a different cell fate model of the following:
!
!     * Critical birth-death model (CBD): each cell choses between dif/prolif stochastically with
!       prob=1/2 for both. Lifetimes are generated stochastically as ~EXP(1/t_cell).
!
!     * Voter Model (VM): each cell stochastically choses one side and evolves with the opposite
!       fate as that side's neighbour cell. Lifetimes are generated stochastically as~EXP(1/t_cell).
!
!     * Density model: cell fate is determined stochastically in function of the local density
!       (at [-l,+L] around the progenitor). With this density, two stochastic lifetimes are
!       generated, following t1~EXP(1/w1*t_cell) and t2~EXP(1/w2*t_cell). If t1<t2 the cell will
!       proliferate after t1, if not, the cell differentiates after t2. This model can be simulated
!       with long and short-range (adhesion) interactions between cells, as well as without any
!       interaction.
!
! Additionally, this library includes subroutines to compute the local density (number of close
! cells) and gaussian random numbers by the Box-Müller method.
!
! Joan Térmens Cascalló,                                                last revision: May 24, 2019
!--------------------------------------------------------------------------------------------------

module ssa_simul7
  implicit none

contains

!--------------------------------------------------------------------------------------------------
!                                     SUBROUTINE main_dens
!
! Simulates the collective behaviour of 1-dim tissue of stem cells with periodic boundary
! conditions. Each cell choses between proliferating or differentiating (leaving the tissue)
! according a density model.
! Usage: main_dens(TestMode,CompMeas,GrMrg,Adh,L,SEED,time_measure,Ncell_rec,Nsurv_rec,frac_rec,k_measure,err)
!        Where TestMode, CompMeas and Adh are logical values that activate a test mode of the
!        simulation, the recording of complementary measures and short-range the interaction
!        potentials to simulate adhesion, respectively. When GrMrg=1 cell growth and migration are
!        computed, whereas when GrMrg=0, they are not. L is the maximum distance at which cells
!        measure density, SEED is the seed of random numbers, time_measure is an array of size
!        Nmeasure+1 with the times at which measures are taken, Ncell_rec, Nsurv_rec and frac_rec
!        are output arrays with the simulation's results. k_measure is a real output that represents
!        the number of measures and err is an error flag.
!
! * Warning, the recording of complementary measures by setting CompMeas = .true. can output large files.
!--------------------------------------------------------------------------------------------------
  !subroutine main_dens(TestMode,CompMeas,GrMrg,ShortRange,L,SEED,time_measure,Ncell_rec,Nsurv_rec,frac_rec,k_measure,err)
  subroutine main_dens(SEED,time_measure,Ncell_rec,Nsurv_rec,frac_rec,k_measure,err)
    use mtmod
    use system_parameters
    implicit none

    !logical,intent(in)::TestMode,CompMeas,ShortRange,GrMrg !when true, activate the test mode, record
            !complementary measures, short-range potentials (adhesion) or growth/migration of cells
    integer(4),intent(in)::SEED !seed for the random number generation
    real(8),intent(in)::time_measure(Nmeasure+1) !times at which Nsurv is recorded
    !real(8),intent(in)::L !maximum distance at which cells measure density

    real(8),intent(out)::frac_rec(Nmeasure+1) !register of frac_diff at time_measure
    integer(4),intent(out)::Nsurv_rec(Nmeasure+1),Ncell_rec(Nmeasure+1) !register of Nsurv at time_measure
    integer(4),intent(out)::k_measure !# Nsurv measures recorded
    logical,intent(out)::err !error flag

    logical::fate,mode,label !fate of a newborned cell & mode of evolution (true =>prolif), label?

    integer(4)::Ncell,Nclose !# cells, #cells at [-L,+L] around a given cell
    integer(4)::Nsurv,last_Nsurv !# labeled cells, # labeled cells before last dif/prolif

    integer(4)::i,i_run,i_time !generic index, index of main loop runs, index of stak's time loop
    integer(4)::index_cell,index1,index2 !index chosen cell, generic indexes of cells

    real(8)::w1,w2 !reaction rates
    real(8)::time,last_time !time simulated, last time of a prolif/dif
    real(8)::t_life,tmin,t1,t2 !lifetime of a cell, minim-themeum of remaining lifetimes
    real(8)::dens !local density
    real(8)::delta_sz,pst_progenitor,sz_progenitor
    real(8)::dist_eq1,dist_eq2 !distance between a cell and the previous one, equilibrium distance
                               !between a cell and the previous one dist_eq=(sz+prev_sz)/2

    real(8)::xgauss,rnd1,rnd2 !random variables to generate white noise
    real(8)::du1,du2,x1,x2 !energy differentials and positions of a pair of cells

    integer(4)::sum_diff,sum_prolif !sum of differentiating and proliferating cells, respectively
    real(8)::frac_diff !fraction of differentiated cells

    logical,allocatable::cell_logic(:,:) !array of cells'logic values: [label,fate] (true =>prolif)
    real(8),allocatable::cell_real(:,:) !array of cells'real values: [lifetime,size,position]

    logical,allocatable::cell_logic_temp(:,:) !temporal arrays used to save values at de/allocate
    real(8),allocatable::cell_real_temp(:,:)

    real(8),allocatable::prev_pst(:),prev_sz(:) !previous position and previous size

    real(4)::TIME1,TIME2,TIME_STACK,TIME_PROLIF,TIME_DIFF,TIME_RECORD,TIME_INI,TIME_TOTAL,TIME0

    character(2)::charact_L ! L and seed in character form to name the output files
    character(8)::charact_seed

    call cpu_time(TIME1)
    call cpu_time(TIME0)

    call sgrnd(SEED) !initialize random numbers with a given seed

    write(charact_L,fmt="(I2)") int(L)
    write(charact_seed,fmt="(I8)") SEED

    if(CompMeas.eqv..true.) then
      open(unit=1,file="dens_L="//charact_L//"_seed="//charact_seed//".dat")
      open(unit=2,file="posis_L="//charact_L//"_seed="//charact_seed//".dat")
    endif

    if(TestMode.eqv..true.) then
      print*," "
      print*,"TEST RESUME, SEED=",SEED," and L=",L
    endif

    Ncell=Nini !initial number of cells
    k_measure=0 !no measures taken
    time=0.d0

    frac_diff=0.d0
    sum_diff=0
    sum_prolif=0

    TIME_INI=0.0
    TIME_STACK=0.0
    TIME_PROLIF=0.0
    TIME_DIFF=0.0
    TIME_RECORD=0.0

    allocate(cell_logic(Ncell,2),cell_real(Ncell,3))

!-------------INITIALIZE CELL ARRAYS---------------------------------------------------------------
! The following two arrays are initializes with Nini cells:
!   * cell_logic(i)=[label,fate](i), where label=.true. -> labeled and fate=.true. -> proliferation
!   * cell_real(i)=[lifetime,size,position](i), where lifetime>t_growth, size=l0, the adult size
!     and position=l0/2+(Ncell-1)*l0.

    do index_cell=1,Ncell
      cell_logic(index_cell,1)=.false.
      cell_real(index_cell,2)=l0 !all cells initialized as adults
      cell_real(index_cell,3)=l0*0.5d0+dble(index_cell-1)*(l0+delta_ini)
    enddo

    index_cell=int(real(Nini)/2.d0)
    cell_logic(index_cell,1)=.true.

    call Ncell_close(L,index_cell,cell_real(:,3),Lsys,Nclose) !compute numerical cell density

    dens=dble(Nclose)/(2.d0*L+1.d0)

    do index_cell=1,Ncell
!------------------COMPUTE CELL FATE---------------------------------------------------------------
! Computation of fate & lifetime for a cell. w1 and w2 characterise the relation with the density
! by modifying the mean life of a cell. The lower lifetime and its related fate are chosen.

      if(dens.gt.4.33d0) then
        if(TestMode.eqv..true.) then
          print*,"WARNING: density fixed from ",dens
        endif
        dens=4.33d0
      endif

      w1=-r*(dens-dens0)/dens0+1.d0 !compute reaction rates
      w2=r*(dens-dens0)/dens0+1.d0

      t1=-(dlog(grnd())/w1)*t_cell !generate two random lifetimes ~ EXP(1/w*t_cell)
      t2=-(dlog(grnd())/w2)*t_cell

      if(t1.lt.t2) then !chose lower lifetime and its fate
        t_life=t1
        fate=.true.
      else
        t_life=t2
        fate=.false.
      endif

!------------------END CELL FATE COMPUTATION-------------------------------------------------------
      cell_logic(index_cell,2)=fate !assign cell fates and lifetimes
      cell_real(index_cell,1)=t_life

    enddo

!-------------END INITIALIZE CELL ARRAYS-----------------------------------------------------------

    Nsurv=1 !start with 1 labeled survivor

    call cpu_time(TIME2)

    TIME_INI=TIME2-TIME1

!-------------MAIN LOOP----------------------------------------------------------------------------
! The main loop is executed until Nsurv=0, Nsurv=Ncell, time>max(time_measure) or max_run is
! reached. Firstly the cell with lower remaining time is chosen (index_cell) and are executed all
! processes in the stak, i.e. all growth and migration processes. Next, the index_cell ends its
! life either proliferation or diferentiating and indexes are updated accordingly. Finally, if
! index_cell was a labeled cell, some output measures of Nsurv, density and time are recorded.

    do i_run=1,max_run

      call cpu_time(TIME1)

      index_cell=minloc(cell_real(:,1),dim=1) !determine next cell to dif/prolif
      tmin=cell_real(index_cell,1) !time remaining to next dif/prolif
      label=cell_logic(index_cell,1)
      mode=cell_logic(index_cell,2)

      last_time=time
      time=time+tmin
      last_Nsurv=Nsurv

!------------------STACK---------------------------------------------------------------------------
! In the stack all growth and migration processes are executed and cell sizes and positions are
! updated accordingly. Left extreme of the fisrt cell is fixed at 0.d0 and the other positions
! change. First, the first cell is analised. If its size isn't the adult size or it distance to the
! previous cell (i.e. the last one) isn't the equilibrium distance, growth and migration processes,
! respectively, are executed. Cell's size and position are updated and the changes are accumulated
! due to it's effect on the next cell. This loop is repeated for each cell. Lastly Lsys is updated.

      allocate(prev_pst(Ncell),prev_sz(Ncell))

      do i_time=1,int(tmin/dt)
        time=time+dt

        prev_sz=cell_real(:,2)
        prev_pst=cell_real(:,3)

        !evolution for the 1st cell

        dist_eq1=0.5d0*(prev_sz(Ncell)+prev_sz(1))
        dist_eq2=0.5d0*(prev_sz(1)+prev_sz(2))

        x1=prev_pst(1)-(prev_pst(Ncell)-Lsys)
        x2=prev_pst(2)-prev_pst(1)

        if (ShortRange) then !short-range the potentials to simulate adhesion
          du1=-K*(x1-dist_eq1)*0.5d0*(dsign(1.d0,(r_c*prev_sz(1)-abs(x1)))+1.d0)
          du2=K*(x2-dist_eq2)*0.5d0*(dsign(1.d0,(r_c*prev_sz(1)-abs(x2)))+1.d0)
        else
          du1=-K*(x1-dist_eq1) !potentials are not short-ranged
          du2=K*(x2-dist_eq2)
        endif

        rnd1=grnd()
        rnd2=grnd()
        call Box_Muller(0.d0,dt,rnd1,rnd2,xgauss)

        if (GrMrg) then !apply migration
          cell_real(1,3)=prev_pst(1)+(du1+du2)*dt+dsqrt(2.d0*D)*xgauss
        else
          cell_real(1,3)=prev_pst(1)+dsqrt(2.d0*D)*xgauss
        endif

        do i=2,(Ncell-1) !migration processes solved by Euler-Maruyama method

          dist_eq1=0.5d0*(prev_sz(i-1)+prev_sz(i))
          dist_eq2=0.5d0*(prev_sz(i)+prev_sz(i+1))

          x1=prev_pst(i)-prev_pst(i-1)
          x2=prev_pst(i+1)-prev_pst(i)

          if (ShortRange) then !short-range the potentials to simulate adhesion
            du1=-K*(x1-dist_eq1)*0.5d0*(dsign(1.d0,(r_c*prev_sz(i)-abs(x1)))+1.d0)
            du2=K*(x2-dist_eq2)*0.5d0*(dsign(1.d0,(r_c*prev_sz(i)-abs(x2)))+1.d0)
          else
            du1=-K*(x1-dist_eq1) !potentials are not short-ranged
            du2=K*(x2-dist_eq2)
          endif

          rnd1=grnd()
          rnd2=grnd()
          call Box_Muller(0.d0,dt,rnd1,rnd2,xgauss) !generate gaussian random numbers

          if (GrMrg) then !apply migration process
            cell_real(i,3)=prev_pst(i)+(du1+du2)*dt+dsqrt(2.d0*D)*xgauss
          else
            cell_real(i,3)=prev_pst(i)+dsqrt(2.d0*D)*xgauss
          endif
        enddo

        !evolution for the N_{cell}th cell

        dist_eq1=0.5d0*(prev_sz(Ncell)+prev_sz(Ncell-1))
        dist_eq2=0.5d0*(prev_sz(1)+prev_sz(Ncell))

        x1=prev_pst(Ncell)-prev_pst(Ncell-1)
        x2=prev_pst(1)-(prev_pst(Ncell)-Lsys)

        if (ShortRange.eqv..true.) then !short-range the potentials to simulate adhesion
          du1=-K*(x1-dist_eq1)*0.5d0*(dsign(1.d0,(r_c*prev_sz(Ncell)-abs(x1)))+1.d0)
          du2=K*(x2-dist_eq2)*0.5d0*(dsign(1.d0,(r_c*prev_sz(Ncell)-abs(x2)))+1.d0)
        else
          du1=-K*(x1-dist_eq1)*0.5d0
          du2=K*(x2-dist_eq2)*0.5d0
        endif

        rnd1=grnd()
        rnd2=grnd()
        call Box_Muller(0.d0,dt,rnd1,rnd2,xgauss)

        if (GrMrg) then
          cell_real(Ncell,3)=prev_pst(Ncell)+(du1+du2)*dt+dsqrt(2.d0*D)*xgauss
        else
          cell_real(Ncell,3)=prev_pst(Ncell)+dsqrt(2.d0*D)*xgauss
        endif
        !Growth processes

        if(GrMrg) then !apply cell growth
          do i=1,Ncell
            if(prev_sz(i).lt.l0) then
              delta_sz=(0.5d0*l0/t_growth)*dt
              if((prev_sz(i)+delta_sz).gt.l0) then
                cell_real(i,2)=l0
              else
                cell_real(i,2)=prev_sz(i)+delta_sz
              endif
            endif
          enddo
        endif

      enddo

      deallocate(prev_pst,prev_sz)

!------------------END STACK-----------------------------------------------------------------------
      call cpu_time(TIME2)

      TIME_STACK=TIME_STACK+(TIME2-TIME1)

      call cpu_time(TIME1)

      if (mode.eqv..true.) then

!------------------PROLIFERATION-------------------------------------------------------------------
! Execution of a proliferation process at cell_index. The number of cells increase by 1. Firstly,
! cell density is computed and cell indexes are updated. All indexes higher than index_cell are
! incremented by one to leave free space to the two descendants. Next, descendants' fates and
! lifetimes are computed. Finally, the descendants are saved on cell arrays, its size is half the
! size of an adult cell (l0/2) and its positions are -l0/2 and +l0/2 from pregenitor cell's
! position, respectively

        sum_prolif=sum_prolif+1

        call Ncell_close(L,index_cell,cell_real(:,3),Lsys,Nclose) !compute cell density
        dens=dble(Nclose)/(2.d0*L+1.d0)

        pst_progenitor=cell_real(index_cell,3)
        sz_progenitor=cell_real(index_cell,2)

!-----------------------UPDATE CELL INDEXES (proliferation)----------------------------------------
! Number of cells increase by 1. Indexes lower than index_cell remain the same and indexes higher
! are increased by 1. Positions at index_cell & index_cell+1 are left free to the descendants.

        Ncell=Ncell+1 !proliferation generates a new cell

        allocate(cell_logic_temp(Ncell,2),cell_real_temp(Ncell,3))

        if(index_cell.ne.1) then
          cell_logic_temp(1:(index_cell-1),:)=cell_logic(1:(index_cell-1),:) !indexes in previous
          cell_real_temp(1:(index_cell-1),:)=cell_real(1:(index_cell-1),:) !cells remain the same
        endif

        if(index_cell.ne.(Ncell-1)) then
          cell_logic_temp((index_cell+2):Ncell,:)=cell_logic((index_cell+1):(Ncell-1),:) !increase in 1 indexes
          cell_real_temp((index_cell+2):Ncell,:)=cell_real((index_cell+1):(Ncell-1),:) !posteriors to cell_index
        endif

        deallocate(cell_real,cell_logic) !add one cell to cell arrays
        allocate(cell_logic(Ncell,2),cell_real(Ncell,3))
        cell_logic=cell_logic_temp
        cell_real=cell_real_temp
        deallocate(cell_logic_temp,cell_real_temp)

!-----------------------END UPDATE CELL INDEXES (proliferation)------------------------------------

        index1=index_cell !indexes for newborned cells
        index2=index_cell+1

        cell_real(:,1)=cell_real(:,1)-tmin !update lifetimes

        do i=index1,index2
!------------------COMPUTE CELL FATE---------------------------------------------------------------
! Computation of fate & lifetime for a cell. w1 and w2 characterise the relation with the density
! by modifying the mean life of a cell. The lower lifetime and its related fate are chosen.

          if(dens.gt.4.33d0) then
            if(TestMode.eqv..true.) then
              print*,"WARNING: density fixed from ",dens
            endif
            dens=4.33d0
          endif

          w1=-r*(dens-dens0)/dens0+1.d0 !compute reaction rates
          w2=r*(dens-dens0)/dens0+1.d0

          t1=-(dlog(grnd())/w1)*t_cell !generate two random lifetimes ~ EXP(1/w*t_cell)
          t2=-(dlog(grnd())/w2)*t_cell

          if(t1.lt.t2) then !chose lower lifetime and its fate
            t_life=t1
            fate=.true.
          else
            t_life=t2
            fate=.false.
          endif

!------------------END CELL FATE COMPUTATION-------------------------------------------------------
          cell_logic(i,2)=fate !assign fate and lifetime to newborned cells
          cell_real(i,1)=t_life
        enddo

        if(label.eqv..true.) then !if the progenitor cell is labeled, descendants will be too
          cell_logic(index1,1)=.true.
          cell_logic(index2,1)=.true.
          Nsurv=Nsurv+1 !number of survivors increase in 1
        else
          cell_logic(index1,1)=.false.
          cell_logic(index2,1)=.false.
        endif

        if (GrMrg) then
          cell_real(index1,2)=sz_progenitor*0.5d0 !newborned cells have half the size of the progenitor
          cell_real(index2,2)=sz_progenitor*0.5d0
        else
          cell_real(index1,2)=sz_progenitor !newborned cells have the size of the progenitor, l0
          cell_real(index2,2)=sz_progenitor
        endif

        cell_real(index1,3)=pst_progenitor-l0*0.25d0 !assign positions to newborned cells
        cell_real(index2,3)=pst_progenitor+l0*0.25d0

!------------------END PROLIFERATION---------------------------------------------------------------

        call cpu_time(TIME2)

        TIME_PROLIF=TIME_PROLIF+(TIME2-TIME1)

      else

!------------------DIFFERENTIATION-----------------------------------------------------------------
! Execution of a differentiation at cell_index. The number of cells decrease by 1 due to one cell
! differentiating and leaving the tissue. The indexes of cells posteriors to index_cell are reduced
! by one, leaving empty the space before occupied by index_cell.

        sum_diff=sum_diff+1

!-----------------------UPDATE CELL INDEX (differentiation)----------------------------------------
! The indexes of cells posteriors to index_cell are reduced by one, the indexes previous to it
! remain the same.

        Ncell=Ncell-1 !differentiation removes a cell from the tissue

        allocate(cell_logic_temp(Ncell,2),cell_real_temp(Ncell,3))

        if(index_cell.ne.1) then
          cell_logic_temp(1:(index_cell-1),:)=cell_logic(1:(index_cell-1),:)
          cell_real_temp(1:(index_cell-1),:)=cell_real(1:(index_cell-1),:)
        endif

        if(index_cell.ne.(Ncell+1)) then
          cell_logic_temp(index_cell:Ncell,:)=cell_logic((index_cell+1):(Ncell+1),:)
          cell_real_temp(index_cell:Ncell,:)=cell_real((index_cell+1):(Ncell+1),:)
        endif

        deallocate(cell_real,cell_logic)
        allocate(cell_logic(Ncell,2),cell_real(Ncell,3))
        cell_logic=cell_logic_temp
        cell_real=cell_real_temp
        deallocate(cell_logic_temp,cell_real_temp)

!-----------------------END UPDATE CELL INDEX (differentiation)------------------------------------
        if(label.eqv..true.) then
          Nsurv=Nsurv-1
        endif

        cell_real(:,1)=cell_real(:,1)-tmin !update lifetimes
      endif

      frac_diff=dble(sum_diff)/(dble(sum_diff)+dble(sum_prolif))

!------------------END DIFFERENTIATION-------------------------------------------------------------
      call cpu_time(TIME2)

      TIME_DIFF=TIME_DIFF+(TIME2-TIME1)

      call cpu_time(TIME1)

!------------------OUTPUT RECORD-------------------------------------------------------------------
! If Nsurv has changed, measures of Nsurv, density and time are taken. Nsurv measures are taken at
! the times determined by time_measure and are stored in a fixed size array. Density and time
! measures are taken each time a labeled cell dif/prolif and are stored in allocatable arrays.

      if(label.eqv..true.) then
        do i=(k_measure+1),(Nmeasure+1)
          if(time_measure(i).lt.time) then
            Nsurv_rec(i)=last_Nsurv
            frac_rec(i)=frac_diff
            Ncell_rec(i)=Ncell
          else
            k_measure=i-1
            exit
          endif
        enddo

        if(CompMeas.eqv..true.) then
          write(1,*) time,dens
          write(2,*) Ncell,time,(merge(1,0,cell_logic(i,1)),cell_real(i,3), i=1,Ncell)
        endif
!------------------END OUTPUT RECORD---------------------------------------------------------------


        if((i.eq.(Nmeasure)).or.(Nsurv.eq.0).or.(Nsurv.eq.Ncell).or.(Ncell.gt.(4*Nini))) then
          if(Ncell.gt.(4*Nini)) then
            err=.true.
          else
            err=.false.
          endif


          call cpu_time(TIME_TOTAL)
          TIME_TOTAL=TIME_TOTAL-TIME0

          if(TestMode.eqv..true.) then
            print*,"  TIME INI=",TIME_INI," frac=",(TIME_INI/TIME_TOTAL)*100.0,"%"
            print*,"  TIME STACK=",TIME_STACK," frac=",(TIME_STACK/TIME_TOTAL)*100.0,"%"
            print*,"  TIME PROLIF=",TIME_PROLIF," frac=",(TIME_PROLIF/TIME_TOTAL)*100.0,"%"
            print*,"  TIME DIFF=",TIME_DIFF," frac=",(TIME_DIFF/TIME_TOTAL)*100.0,"%"
            print*,"  TIME REC=",TIME_RECORD," frac=",(TIME_RECORD/TIME_TOTAL)*100.0,"%"
            print*,"  # MEASURES=",k_measure
            print*," "
          endif
          if(CompMeas.eqv..true.) then
            close(1)
            close(2)
          endif
          return
        endif
      endif
!-------------END MAIN LOOP------------------------------------------------------------------------
      call cpu_time(TIME2)
      TIME_RECORD=TIME_RECORD+(TIME2-TIME1)

      if((TestMode.eqv..true.).and.(Ncell.gt.4*Nini)) then
        print*, "WARNING: Ncell=",Ncell,">",4*Nini
      endif
    enddo

    if(i_run.eq.max_run) then

      call cpu_time(TIME_TOTAL)

      if(TestMode.eqv..true.) then
        print*, "max_run reached, program finished"
        print*," "
        print*,"  TIME INI=",TIME_INI," frac=",(TIME_INI/TIME_TOTAL)*100.0,"%"
        print*,"  TIME STACK=",TIME_STACK," frac=",(TIME_STACK/TIME_TOTAL)*100.0,"%"
        print*,"  TIME PROLIF=",TIME_PROLIF," frac=",(TIME_PROLIF/TIME_TOTAL)*100.0,"%"
        print*,"  TIME DIFF=",TIME_DIFF," frac=",(TIME_DIFF/TIME_TOTAL)*100.0,"%"
        print*,"  TIME REC=",TIME_RECORD," frac=",(TIME_RECORD/TIME_TOTAL)*100.0,"%"
        print*,"  # MEASURES=",k_measure
        print*," "
      endif
      if(CompMeas.eqv..true.) then
        close(1)
        close(2)
      endif
      return
    endif

  end subroutine main_dens


!--------------------------------------------------------------------------------------------------
!                                     SUBROUTINE main_CBD
!
! Simulates the collective behaviour of 1-dim tissue of stem cells with periodic boundary
! conditions. Each cell choses between proliferating or differentiating (leaving the tissue)
! according to the Critical Birth-Death Model (CBD).
! Usage: main_CBD(SEED,time_measure,Nsurv_rec,frac_rec,k_measure)
!        Where SEED is the seed of random numbers, time_measure is an array of size Nmeasure+1 with
!        the times at which measures are taken, Nsurv_rec and frac_rec are the output arrays with
!        the simulation's results and k_measure is the number of measures recorded.
!
! * In the CBD model the position and size of the cells isn't relevant, so descendants born with
! adult size and stack simulation is omited.
!--------------------------------------------------------------------------------------------------
  subroutine main_CBD(SEED,time_measure,Nsurv_rec,frac_rec,k_measure)
    use mtmod
    use system_parameters
    implicit none

    integer(4),intent(in)::SEED
    real(8),intent(in)::time_measure(Nmeasure+1) !times at which Nsurv is recorded

    integer(4),intent(out)::Nsurv_rec(Nmeasure+1) !register of Nsurv at time_measure
    real(8),intent(out)::frac_rec(Nmeasure+1)
    integer(4),intent(out)::k_measure !# Nsurv measures recorded

    logical::fate,mode,label !fate of a newborned cell & mode of evolution (true =>prolif), label

    integer(4)::Ncell !# cells
    integer(4)::Nsurv,last_Nsurv !# labeled cells, # labeled cells before last dif/prolif

    integer(4)::i,i_run !generic index, index of main loop runs
    integer(4)::index_cell,index1,index2 !index chosen cell, generic indexes of cells

    real(8)::time,last_time !time simulated, last time of a prolif/dif
    real(8)::t_life,tmin !lifetime of a cell, minimum of remaining lifetimes
    real(8)::rnd !a random number

    integer(4)::sum_diff,sum_prolif
    real(8)::frac_diff

    logical,allocatable::cell_logic(:,:) !array of cells'logic values: [label,fate] (true =>prolif)
    real(8),allocatable::cell_real(:) !array of cells'real values: [lifetime]

    logical,allocatable::cell_logic_temp(:,:) !temporal arrays used tosize of the system, save values at de/allocate
    real(8),allocatable::cell_real_temp(:)

    call sgrnd(SEED) !initialize random numbers with a given seed

    Ncell=Nini !initial number of cells
    k_measure=0 !no measures taken
    time=0.d0

    sum_diff=0
    sum_prolif=0
    frac_diff=0.d0

    allocate(cell_logic(Ncell,2),cell_real(Ncell))

!-------------INITIALIZE CELL ARRAYS---------------------------------------------------------------
! The following two arrays are initializes with Nini cells:
!   * cell_logic(i)=[label,fate](i), where label=.true. -> labeled and fate=.true. -> proliferation
!   * cell_real(i)=[lifetime,size,position](i), where lifetime>t_growth, size=l0, the adult size
!     and position=l0/2+(Ncell-1)*l0.

    index_cell=1
    cell_logic(index_cell,1)=.true. !cell 1 begins as labeled
    do index_cell=2,Ncell !other cells
      cell_logic(index_cell,1)=.false.
    enddo


    do index_cell=1,Ncell
!------------------COMPUTE CELL FATE---------------------------------------------------------------
! Computation of fate & lifetime for a cell. Cell lifetime is stochastic exponentially distributed
! and cell fate is stochastic with 1/2 chances of being differentiation, 1/2 of being proliferation.

      t_life=(-dlog(grnd())*t_cell)+t_growth

      rnd=grnd()

      if(rnd.gt.0.5d0) then
        fate=.true.
      else
        fate=.false.
      endif

!------------------END CELL FATE COMPUTATION-------------------------------------------------------
      cell_logic(index_cell,2)=fate !assign cell fates and lifetimes
      cell_real(index_cell)=t_life
    enddo
!-------------END INITIALIZE CELL ARRAYS-----------------------------------------------------------

    Nsurv=1 !start with 1 labeled survivor

!-------------MAIN LOOP----------------------------------------------------------------------------
! The main loop is executed until Nsurv=0, Nsurv=Ncell, time>max(time_measure) or max_run is
! reached. Firstly the cell with lower remaining time is chosen (index_cell) and are executed all
! processes in the stak, i.e. all growth and migration processes. Next, the index_cell ends its
! life either proliferation or diferentiating and indexes are updated accordingly. Finally, if
! index_cell was a labeled cell, some output measures of Nsurv, density and time are recorded.

    do i_run=1,max_run
      index_cell=minloc(cell_real,dim=1) !determine next cell to dif/prolif
      tmin=cell_real(index_cell) !time remaining to next dif/prolif
      label=cell_logic(index_cell,1)
      mode=cell_logic(index_cell,2)

      last_time=time
      time=time+tmin
      last_Nsurv=Nsurv

      if (mode.eqv..true.) then

!------------------PROLIFERATION-------------------------------------------------------------------
! Execution of a proliferation process at cell_index. The number of cells increase by 1. Firstly,
! cell density is computed and cell indexes are updated. All indexes higher than index_cell are
! incremented by one to leave free space to the two descendants. Next, descendants' fates and
! lifetimes are computed. Finally, the descendants are saved on cell arrays, its size is half the
! size of an adult cell (l0/2) and its positions are -l0/2 and +l0/2 from pregenitor cell's
! position, respectively

        sum_prolif=sum_prolif+1

!-----------------------UPDATE CELL INDEXES (proliferation)----------------------------------------
! Number of cells increase by 1. Indexes lower than index_cell remain the same and indexes higher
! are increased by 1. Positions at index_cell & index_cell+1 are left free to the descendants.

        Ncell=Ncell+1 !proliferation generates a new cell

        allocate(cell_logic_temp(Ncell,2),cell_real_temp(Ncell))

        if(index_cell.ne.1) then
          cell_logic_temp(1:(index_cell-1),:)=cell_logic(1:(index_cell-1),:) !indexes in previous
          cell_real_temp(1:(index_cell-1))=cell_real(1:(index_cell-1)) !cells remain the same
        endif

        if(index_cell.ne.(Ncell-1)) then
          cell_logic_temp((index_cell+2):Ncell,:)=cell_logic((index_cell+1):(Ncell-1),:) !increase in 1 indexes
          cell_real_temp((index_cell+2):Ncell)=cell_real((index_cell+1):(Ncell-1)) !posteriors to cell_index
        endif

        deallocate(cell_real,cell_logic) !add one cell to cell arrays
        allocate(cell_logic(Ncell,2),cell_real(Ncell))
        cell_logic=cell_logic_temp
        cell_real=cell_real_temp
        deallocate(cell_logic_temp,cell_real_temp)

!-----------------------END UPDATE CELL INDEXES (proliferation)------------------------------------

        index1=index_cell !indexes for newborned cells
        index2=index_cell+1

        do i=index1,index2
!------------------COMPUTE CELL FATE---------------------------------------------------------------
! Computation of fate & lifetime for a cell. Cell lifetime is stochastic exponentially distributed
! and cell fate is stochastic with 1/2 chances of being differentiation, 1/2 of being proliferation.

      t_life=(-dlog(grnd())*t_cell)+t_growth

      rnd=grnd()

      if(rnd.gt.0.5d0) then
        fate=.true.
      else
        fate=.false.
      endif

!------------------END CELL FATE COMPUTATION-------------------------------------------------------
          cell_logic(i,2)=fate !assign fate and lifetime to newborned cells
          cell_real(i)=t_life
        enddo

        if(label.eqv..true.) then !if the progenitor cell is labeled, descendants will be too
          cell_logic(index1,1)=.true.
          cell_logic(index2,1)=.true.
          Nsurv=Nsurv+1 !number of survivors increase in 1
        else
          cell_logic(index1,1)=.false.
          cell_logic(index2,1)=.false.
        endif

!------------------END PROLIFERATION---------------------------------------------------------------
      else
!------------------DIFFERENTIATION-----------------------------------------------------------------
! Execution of a differentiation at cell_index. The number of cells decrease by 1 due to one cell
! differentiating and leaving the tissue. The indexes of cells posteriors to index_cell are reduced
! by one, leaving empty the space before occupied by index_cell.

        sum_diff=sum_diff+1

!-----------------------UPDATE CELL INDEX (differentiation)----------------------------------------
! The indexes of cells posteriors to index_cell are reduced by one, the indexes previous to it
! remain the same.

        Ncell=Ncell-1 !differentiation removes a cell from the tissue

        allocate(cell_logic_temp(Ncell,2),cell_real_temp(Ncell))

        if(index_cell.ne.1) then
          cell_logic_temp(1:(index_cell-1),:)=cell_logic(1:(index_cell-1),:)
          cell_real_temp(1:(index_cell-1))=cell_real(1:(index_cell-1))
        endif

        if(index_cell.ne.(Ncell+1)) then
          cell_logic_temp(index_cell:Ncell,:)=cell_logic((index_cell+1):(Ncell+1),:)
          cell_real_temp(index_cell:Ncell)=cell_real((index_cell+1):(Ncell+1))
        endif

        deallocate(cell_real,cell_logic)
        allocate(cell_logic(Ncell,2),cell_real(Ncell))
        cell_logic=cell_logic_temp
        cell_real=cell_real_temp
        deallocate(cell_logic_temp,cell_real_temp)

!-----------------------END UPDATE CELL INDEX (differentiation)------------------------------------
        if(label.eqv..true.) then
          Nsurv=Nsurv-1
        endif

      endif
!------------------END DIFFERENTIATION-------------------------------------------------------------

      frac_diff=dble(sum_diff)/(dble(sum_diff)+dble(sum_prolif))

!------------------OUTPUT RECORD-------------------------------------------------------------------
! If Nsurv has changed, measures of Nsurv, density and time are taken. Nsurv measures are taken at
! the times determined by time_measure and are stored in a fixed size array. Density and time
! measures are taken each time a labeled cell dif/prolif and are stored in allocatable arrays.

      if(label.eqv..true.) then
        do i=(k_measure+1),(Nmeasure+1)
          if(time_measure(i).lt.time) then
            Nsurv_rec(i)=last_Nsurv
            frac_rec(i)=frac_diff
          else
            k_measure=i-1
            exit
          endif
        enddo

        if((i.eq.(Nmeasure+1)).or.(Nsurv.eq.0).or.(Nsurv.eq.Ncell)) then
          deallocate(cell_logic,cell_real)
          return
        endif
      endif
!------------------END OUTPUT RECORD---------------------------------------------------------------

      cell_real=cell_real-tmin !update lifetimes


!-------------END MAIN LOOP------------------------------------------------------------------------
    enddo

    if(i_run.eq.max_run) then
      deallocate(cell_logic,cell_real)
      print*, "max_run reached, program finished"
      return
    endif

  end subroutine main_CBD


!--------------------------------------------------------------------------------------------------
!                                     SUBROUTINE main_VM
!
! Simulates the collective behaviour of 1-dim tissue of stem cells with periodic boundary
! conditions. Each cell choses between proliferating or differentiating (leaving the tissue)
! according to the Voter Model (VM).
! Usage: main_CBD(SEED,time_measure,Nsurv_rec,k_measure)
!        Where SEED is the seed of random numbers, time_measure is an array of size Nmeasure+1 with
!        the times at which measures are taken, Nsurv_rec is the output array with the simulation's
!        results and k_measure is the number of measures recorded.
!--------------------------------------------------------------------------------------------------
  subroutine main_VM(SEED,time_measure,Nsurv_rec,k_measure)
    use mtmod
    use system_parameters
    implicit none

    integer(4),intent(in)::SEED
    real(8),intent(in)::time_measure(Nmeasure+1) !times at which Nsurv is recorded

    integer(4),intent(out)::Nsurv_rec(Nmeasure+1) !register of Nsurv at time_measure
    integer(4),intent(out)::k_measure !# Nsurv measures recorded

    logical::label1,label2

    integer(4)::Ncell,Ntime !# cells
    integer(4)::Nsurv,last_Nsurv !# labeled cells, # labeled cells before last dif/prolif

    integer(4)::i,i_run, index_time !generic index, index of main loop runs
    integer(4)::index_cell,index_progenitor !index chosen cell, index of the progenitor cell

    real(8)::time, tmin(2) !time simulated, times of the next differentiation/proliferation
    real(8)::rnd !a random number

    logical::cell_logic(Nini) !array of cells'logic values [label]
    real(8)::time_array(Nini)

    call sgrnd(SEED) !initialize random numbers with a given seed

    Ncell=Nini !initial number of cells
    k_measure=0 !no measures taken
    time=0.d0
    Ntime=2*Ncell

    index_cell=1
    cell_logic(index_cell)=.true.
    do index_cell=2,Ncell
      cell_logic(index_cell)=.false.
    enddo

    do i=1,Ncell
      time_array(i)=-dlog(grnd())*t_cell+t_growth
    enddo

    Nsurv=1

!-------------MAIN LOOP----------------------------------------------------------------------------
    do i_run=1,max_run
      last_Nsurv=Nsurv

      do i=1,2
        index_time=minloc(time_array, dim=1)
        tmin(i)=time_array(index_time)

        time_array=time_array-tmin(i)
        time_array(index_time)=-dlog(grnd())*t_cell+t_growth
      enddo

      time=time+tmin(1)

      index_cell=1+floor(Ncell*grnd())

      rnd=grnd()
      if(rnd.ge.0.5d0) then
        if(index_cell.eq.Ncell) then
          index_progenitor=1
        else
          index_progenitor=index_cell+1
        endif
      else
        if(index_cell.eq.1) then
          index_progenitor=Ncell
        else
          index_progenitor=index_cell-1
        endif
      endif

      label1=cell_logic(index_cell)

      if(label1.eqv..true.) then
        Nsurv=Nsurv-1
!------------------OUTPUT RECORD-------------------------------------------------------------------
! If Nsurv has changed, measures of Nsurv, density and time are taken. Nsurv measures are taken at
! the times determined by time_measure and are stored in a fixed size array. Density and time
! measures are taken each time a labeled cell dif/prolif and are stored in allocatable arrays.

        do i=(k_measure+1),(Nmeasure+1)
          if(time_measure(i).lt.time) then
            Nsurv_rec(i)=last_Nsurv
          else
            k_measure=i-1
            exit
          endif
        enddo

        if((i.eq.(Nmeasure+1)).or.(Nsurv.eq.0).or.(Nsurv.eq.Ncell)) then
          return
        endif
!------------------END OUTPUT RECORD---------------------------------------------------------------
      endif

      last_Nsurv=Nsurv
      time=time+tmin(2)

      cell_logic(index_cell)=cell_logic(index_progenitor)

      label2=cell_logic(index_cell)

      if(label2.eqv..true.) then
        Nsurv=Nsurv+1
!------------------OUTPUT RECORD-------------------------------------------------------------------
! If Nsurv has changed, measures of Nsurv, density and time are taken. Nsurv measures are taken at
! the times determined by time_measure and are stored in a fixed size array. Density and time
! measures are taken each time a labeled cell dif/prolif and are stored in allocatable arrays.

        do i=(k_measure+1),(Nmeasure+1)
          if(time_measure(i).lt.time) then
            Nsurv_rec(i)=last_Nsurv
          else
            k_measure=i-1
            exit
          endif
        enddo

        if((i.eq.(Nmeasure+1)).or.(Nsurv.eq.0).or.(Nsurv.eq.Ncell)) then
          return
        endif
!------------------END OUTPUT RECORD---------------------------------------------------------------
      endif

!-------------END MAIN LOOP------------------------------------------------------------------------
    enddo

    if(i_run.eq.max_run) then
      print*, "max_run reached, program finished"
      return
    endif

  end subroutine main_VM


!--------------------------------------------------------------------------------------------------
!                                    SUBROUTINE Nclose
!
! Computes the number of cells in the interval [-L,+L] around a given cell, index_cell. Returns its
! value af Nclose.
! Usage: Ncell_close(L,index_cell,cell_positions,Lsys,Nclose)
!        Where L is the maximum distance at which cells measure density, index_cell is the index of
!        a given cell, cell_positions=cell_real(:,3), Lsys is the system's size and Nclose
!        (output) is the numer of cells in [-L,+L] around index_cell.
!--------------------------------------------------------------------------------------------------
  subroutine Ncell_close(L,index_cell,cell_positions,Lsys,Nclose)
    implicit none

    real(8),intent(in)::L,Lsys
    integer(4),intent(in)::index_cell
    real(8),intent(in)::cell_positions(:)

    integer(4),intent(out)::Nclose

    integer(4)::i,Ncell,sum,index_process
    real(8)::xleft,xright ![xleft,xright]=[position(index_cell)-L,position(index_cell)+L]

    sum=0 !initialize sum
    Ncell=size(cell_positions)

    xleft=cell_positions(index_cell)-L !calculation of xleft and xright
    xright=cell_positions(index_cell)+L

    if(xleft.lt.0.d0) then !if index_cell is sufficiently close to 1 xleft could be <0
      index_process=1
      xleft=xleft+Lsys !so xleft is normalized to fit in [0,Lsys]

      do i=index_cell,Ncell
        if(cell_positions(i).le.xright) then !sum of cells equal or posterior to index_cell
          sum=sum+1
        else
          exit
        endif
      enddo

      sum=sum+(index_cell-1) !if xleft<0, all cells with positions lower than index cell are
                             !are included
      do i=Ncell,(index_cell+1),-1 !sum of remaining cells positioned after Lsys+xleft
        if(cell_positions(i).ge.xleft) then
          sum=sum+1
        else
          exit
        endif
      enddo

    else if(xright.gt.Lsys) then !analogous loop in the situation where xright>Lsys
      index_process=2
      xright=xright-Lsys

      sum=Ncell-(index_cell-1)

      do i=1,index_cell
        if(cell_positions(i).le.xright) then
          sum=sum+1
        else
          exit
        endif
      enddo

      do i=(index_cell-1),1,-1
        if(cell_positions(i).ge.xleft) then
          sum=sum+1
        else
          exit
        endif
      enddo

    else
      index_process=3
      do i=index_cell,Ncell !normal loop where both extrema fit into [0,Lsys]
        if(cell_positions(i).le.xright) then
          sum=sum+1
        else
          exit
        endif
      enddo

      do i=(index_cell-1),1,-1
        if(cell_positions(i).ge.xleft) then
          sum=sum+1
        else
          exit
        endif
      enddo
    endif

    Nclose=sum

    return
  end subroutine Ncell_close

!--------------------------------------------------------------------------------------------------
!                                    SUBROUTINE Box_Muller
!
! This subroutine generates a random gaussian number xgauss, with mean xmu and
! variance xsigma^2 with the Box-Müller method.
!--------------------------------------------------------------------------------------------------
  subroutine Box_Muller(xmu,xsigma,rnd1,rnd2,xgauss)

    implicit none

    real(8),intent(in)::xmu,xsigma !input statistical moments
    real(8),intent(in)::rnd1,rnd2 !input uniform random numbers
    real(8),intent(out)::xgauss !output gaussian random number

    real(8)::r,phi !random radius and random phase
    real(8),parameter::pi=4.D0*datan(1.D0)

    r=dsqrt(-2.d0*log(rnd1)) ! random radius
    phi=2.d0*pi*rnd2 ! random phase
    xgauss=xmu+xsigma*(r*dcos(phi))
    return
  end subroutine Box_Muller

end module ssa_simul7
