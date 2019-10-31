program main
use global_variables
use initialize_update
use compute_energy_ewald
use input_output
implicit none

  !##########Data Dictionary############!
    integer :: i, j, k
    real*8  :: EE=0, EE1=0, EE2=0, st, fn
    real*8  :: DeltaE, time(3)
  !#####################################!

  !################begin#################!
  call cpu_time(started)
  call random_seed()
  !
  !input and initialize system, timing and histogram parameters.
  call initialize_parameters
  !
  !initialize energy and parameters of potential
  call initialize_energy_parameters_Ewald
  !
  ! restart or continue
  if (restart_or_continue == 0) then
    i=1
    !
    !initialize position  
    call initialize_position
    call Build_Charge_Ewald
    call write_pos
    call write_pos1(i)
    call write_hist    
    !
    !initialize energy array related to pos array
    call initialize_energy_arrays_ewald
    !
    !error analysis
    call error_analysis(0,EE)
  else if (restart_or_continue /= 0) then
    !
    !read position and histogram data
    call continue_read_data(i)
    !
    !initialize energy array related to pos array
    call initialize_energy_arrays_ewald
    !
    !error analysis
    call error_analysis(0,EE)
  end if

  call SPME_Ewald(energy_long0)
  call real_energy_spme(energy_short)

  EE2 = EE

  !##############preheating##############!
  if ( i <= StepNum0 ) then
    do step=i, StepNum0
      call monte_carlo_move(EE, DeltaE)
      if ( mod(step, multistep) == 0 ) then
       call update_multistep(EE, EE2)
      end if
      if ( mod(step,DeltaStep1) == 0 ) then
        call compute_physical_quantities
        call write_physical_quantities( step )
      end if
      if ( mod(step,DeltaStep3) == 0 ) then
        call error_analysis(0, EE1)
        call write_energy(step,EE,EE1)
        call write_pos1(step)
        EE = EE1
        EE2 = EE
      end if
    end do
    i=step
  end if

  !################running###############!
  do step=i,StepNum0+StepNum            
    call monte_carlo_move(EE, DeltaE)
    if ( mod(step, multistep) == 0 ) then
      call update_multistep(EE, EE2)
    end if
    if ( mod(step,DeltaStep1) == 0 ) then
      call compute_physical_quantities
      call write_physical_quantities( step )
    end if
    if ( mod(step,DeltaStep2) == 0 ) then     
      call histogram
    end if
    if ( mod(step,DeltaStep3) == 0 ) then
      call write_pos1(step)               
      call write_hist
    end if
  end do

  !##################end#################!
  write(*,*)
  write(*,*)
  write(*,*) '************************************************************'
  write(*,*) '*************************Finished***************************'
  write(*,*) 'Finished!     '
  write(*,*) 'time:(minutes)', real(total_time/60), 'minutes'
  write(*,*) 'time:(hours)  ', real(total_time/3600), 'hours'
  write(*,*) 'time:(days)   ', real(total_time/86400), 'days'
  write(*,*) '************************************************************'
  write(*,*)
  write(*,*)

end program main



