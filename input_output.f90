module input_output
implicit none 
  ! 1. Initialize parameters and allocate arrays
  ! 2. Continue read data after power off
  ! 3. Calculate physical quantities and histogram
  ! 4. Write physical quantities and histogram to file 

save
  !mass center height of the star
  real*8 :: hh
  real*8 :: h_end
  real*8 :: h_branch
  !height distribution of star, star branch and end
  integer, allocatable, dimension(:,:), private :: phi_s
  integer, allocatable, dimension(:,:), private :: phi_sb
  integer, allocatable, dimension(:,:), private :: phi_se
  integer, allocatable, dimension(:,:), private :: phi_a    ! anions (p+salt)
  integer, allocatable, dimension(:,:), private :: phi_as   ! anions (salt)
  integer, allocatable, dimension(:,:), private :: phi_ap   ! anions (p)
  integer, allocatable, dimension(:,:), private :: phi_i    ! ions(polymer)
  integer, allocatable, dimension(:,:), private :: phi_is   ! salt ions
  integer, allocatable, dimension(:,:), private :: phi_xy   !xy的分布图
  integer, allocatable, dimension(:,:), private :: phi_zx   !xz的分布图

contains

subroutine initialize_parameters
  !--------------------------------------!
  !Initialize system parameters
  !
  !Input
  !  none
  !Output
  !  none
  !External Variables
  !  none
  !Routine Referenced:
  !  none
  !--------------------------------------!
  implicit none
  !
  !read parameters from file
  call read_data
  !
  !obtain other parameters from the read parameters
  call data_operation
  !
  !write parameters to screen
  call write_data
  !
  ! allocate arrays
  call data_allocate
  !
  ! initialize 108 bond vectors
  call initialize_bond_vector
  !
  ! initialize new direction array
  call initialize_new_direction
  !
  ! bonds after move
  call initialize_move

end subroutine initialize_parameters


subroutine read_data
  !--------------------------------------!
  !read system parameters
  !
  !Input
  !  none
  !Output
  !  none
  !External Variables
  !  parameters in global variables
  !Routine Referenced:
  !  none
  !--------------------------------------!
  use global_variables
  implicit none
  !
  !read system data
  open(10,file='./system_data.txt')
    read(10,*) pH_or_not
    read(10,*) Lz2
    read(10,*) Z_empty  
    read(10,*) sigmag   
    read(10,*) Beta       
    read(10,*) qq
    read(10,*) qqi
    read(10,*) ion_ratio
    read(10,*) arm
    read(10,*) Nma
    read(10,*) Nga
    read(10,*) man_s
    read(10,*) multistep 
    read(10,*) StepNum0           
    read(10,*) StepNum   
    read(10,*) DeltaStep        
    read(10,*) DeltaStep1
    read(10,*) DeltaStep2     
    read(10,*) DeltaStep3 
    read(10,*) pH_pKa        
  close(10)
end subroutine read_data


subroutine data_operation
  !--------------------------------------!
  !Initialize system parameters
  !and judge whether restarted or continue
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  !
  !--------------------------------------!
  use global_variables
  implicit none
  logical alive
  integer :: i, nx, dLx, Charge_ions

  !
  !Monomers of each star
  Ns = Nma*arm + 1
  !
  !The total monomers of star brushes
  Npe = Ns * Nga
  !
  !The total charged monomers on PE
  if ( abs(qq) == 0 ) then
    Nq_PE = 0
  else
    if ( man_s /= 0) then
      ! the anchored particle without charge,
      ! each arm is allocated charge independently
      Nq_PE = Nma / man_s * arm * Nga     
    else
      Nq_PE = 0 
    end if
  end if
  !
  !the number of Salt monomers
  if ( abs(qqi) == 0 ) then
    Nq_salt_ions = 0
  else
    Charge_ions  = nint( ion_ratio * Nq_PE * abs(qq) )
    Nq_salt_ions = Charge_ions / abs(qqi)
  end if
  !
  !The total charged particles Nq and total particles NN in system
  Nq = Nq_PE * ( abs(qq)+1 ) + Nq_salt_ions * ( abs(qqi) + 1 )
  NN = Npe + Nq_PE * abs(qq) + Nq_salt_ions * ( abs(qqi) + 1 )
  Nq_net = Nq
  !
  !System size, keep mod(Lx,nx)=0
  Lx2 = nint(sqrt( Nga / sigmag )) * 2
  nx = nint( sqrt(1.*Nga) )    
  Lx2 = Lx2 - mod( Lx2, nx )
  Ly2 = Lx2
  Z_empty = ( 1.*Lz2 + 1.*Lz2 * Z_empty ) / (1.*Lz2)
  Lx = Lx2/2.D0
  Ly = Ly2/2.D0
  Lz = Lz2/2.D0
  sigmag1 = 1.D0 * Nga / Lx / Ly  ! true grafting density
  !
  !number of bonds in system
  N_bond = Nma * arm * Nga
  !
  !whether continue or restart
  Inquire( file='start_time.txt', exist=alive )
  if (alive) then
    open(11,file='./start_time.txt')
      read(11,*) restart_or_continue
    close(11)
  else
    restart_or_continue = 0
  end if

  accept = 0
  accept_pH = 0
  accept_long = 0
  num_move = 0
  num_pH = 0
  num_long = 0
  total_accept = 0
  total_accept_long = 0
  total_accept_pH = 0
end subroutine data_operation


subroutine write_data
  !--------------------------------------!
  !Write parameters to screen
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  !
  !--------------------------------------!
  use global_variables
  implicit none

  write(*,*)
  write(*,*)
  write(*,*) '****************************************************************'
  write(*,*) '*************************system_data****************************'
  write(*,*) 'Incert space in z direction,       Z_empty:', Z_empty-1
  write(*,*) 'Grafting density,                   sigmag:', sigmag
  write(*,*) 'True grafting density,             sigmag1:', sigmag1
  write(*,*) 'Beta=1/kT,                            Beta:', Beta
  write(*,*) 'Charges of polymer,                     qq:', qq
  write(*,*) 'Charges of salt ions,                  qqi:', qqi
  write(*,*) 'Arms of star brushes                   Arm:', arm
  write(*,*) 'Monomers of each arm                   Nma:', Nma
  write(*,*) 'Number of grafted star chains          Nga:', Nga
  write(*,*) 'Each man_s monomers with one charge, man_s:', man_s
  write(*,*) 'total particles,                        NN:', NN
  write(*,*) 'total charged particles,                Nq:', Nq
  write(*,*) 'total charged particles in polymer,  Nq_PE:', Nq_PE
  write(*,*) 'total brushes particles,               Npe:', Npe
  write(*,*) 'total charged salt particles: Nq_salt_ions:', Nq_salt_ions
  write(*,*) 'Lattice number in x direction,         Lx2:', Lx2
  write(*,*) 'Lattice number in y direction,         Ly2:', Ly2
  write(*,*) 'Lattice number in z direction,         Lz2:', Lz2
  write(*,*) 'Number of bonds in system,          N_bond:', N_bond
  write(*,*) 'pH-pKa,                             pH-pKa:', pH_pKa
  write(*,*) '****************************************************************'

  write(*,*)
  write(*,*) '****************************************************************'
  write(*,*) '************************running_steps***************************'
  write(*,*) 'restart (0), continue (0), restart_continue:',restart_or_continue
  write(*,*) 'Preheating steps                   StepNum0:', StepNum0
  write(*,*) 'Running steps                       StepNum:', StepNum
  write(*,*) 'Total steps                StepNum0+StepNum:', (StepNum0+StepNum)
  write(*,*) 'Step inteval, physical quantum,   DeltaStep:', DeltaStep
  write(*,*) 'Step inteval, physical quantum,  DeltaStep1:', DeltaStep1
  write(*,*) 'Step inteval, histogram,         DeltaStep2:', DeltaStep2
  write(*,*) 'Step inteval, output data,       DeltaStep3:', DeltaStep3
  write(*,*) '****************************************************************'
  write(*,*)
  write(*,*)

end subroutine write_data


subroutine data_allocate
  !--------------------------------------!
  !Allocate pos, vel, acc and histogram arrays
  !Initialize histogram arrays
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  Arrays in Global Variables
  !Routine Referenced:
  !  David P. Landau, Kurt Binder. A Guide to Monte Carlo Simulation in
  !  Statistical Physics, 3rd (Cambridge, UK: Cambridge University Press, 
  !  2009), 453.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i
  !
  !position, velocity, and acceleration
  allocate( pos(NN,4) )
  allocate( latt(Lx2,Ly2,Lz2+1) )
  allocate( pos0(NN,4) )
  allocate( latt0(Lx2,Ly2,Lz2+1) )
  pos = 0
  latt = 0
  pos0 = 0
  latt0 = 0

  allocate( ipx(Lx2) )
  allocate( ipy(Ly2) )
  allocate( ipz(Lz2) )

  allocate( ip2x(Lx2) )
  allocate( ip2y(Ly2) )
  allocate( ip2z(Lz2) )

  allocate( imx(Lx2) )
  allocate( imy(Ly2) ) 
  allocate( imz(Lz2) )

  allocate( bonds(108,3) )
  allocate( move(108,6) )
  allocate( b1(108) )
  allocate( b12(108) )
  allocate( monbd(Npe, arm+1) )
  allocate( bond_numb(N_Bond) )
  allocate( bond_numb0(N_Bond) )
  allocate( phi_xy(SizeHist,SizeHist))
  allocate( phi_zx(SizeHist,SizeHist))
  phi_xy = 0
  phi_zx = 0
  monbd = 0

  do i = 1, Lx2
    ipx(i) = i+1
    ip2x(i) = i+2
    imx(i) = i-1
  end do
  ipx(Lx2) = 1       !periodical boundary condition
  ip2x(Lx2-1) = 1
  ip2x(Lx2) = 2
  imx(1) = Lx2

  do i = 1, Ly2
    ipy(i) = i+1
    ip2y(i) = i+2
    imy(i) = i-1
  end do
  ipy(Ly2) = 1       !periodical boundary condition
  ip2y(Ly2-1) = 1
  ip2y(Ly2) = 2
  imy(1) = Ly2

  do i = 1, Lz2      !finite boundary condition
    ipz(i) = i+1
    ip2z(i) = i+2
    imz(i) = i-1
  end do 
  ip2z(Lz2) = Lz2 + 1
  imz(1) = 1   

  !
  !Allocate arrays and initialize them
  allocate( phi_s(Lz2, 2) )
  allocate( phi_sb(Lz2,2) )
  allocate( phi_se(Lz2,2) )  
  allocate( phi_a(Lz2, 2) )
  allocate( phi_as(Lz2,2) )
  allocate( phi_ap(Lz2,2) ) 
  allocate( phi_i(Lz2, 2) )
  allocate( phi_is(Lz2,2) )
  phi_s = 0
  phi_sb = 0
  phi_se = 0 
  phi_a = 0
  phi_as = 0
  phi_ap = 0
  phi_i = 0
  phi_is = 0

end subroutine data_allocate


subroutine initialize_bond_vector
  !--------------------------------------!
  ! The bonds.txt file are generated by MATLAB according to the subroutine
  ! bdibfl in pp. 450-452 in Landau's book.
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  Arrays in Global Variables
  !Routine Referenced:
  !  David P. Landau, Kurt Binder. A Guide to Monte Carlo Simulation in
  !  Statistical Physics, 3rd (Cambridge, UK: Cambridge University Press, 
  !  2009), see Subroutine bdibfl in pp. 450-452.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j

  !
  !bonds.txt is generated by matlab by referring Landau's book.
  open(101,file='bonds.txt')
    read(101,*) ((bonds(i,j),j=1,3),i=1,108)
  close(101)
  !
  !testify bonds array
  !   do i = 1, 108
  !    write(*,*) bonds(i,:)
  !   end do

  do i = 1, 108
    b12(i) = dot_product(bonds(i,:), bonds(i,:))
    b1(i)  = sqrt(b12(i))
  end do
  !
  !testify b1, b12 array
  !   do i = 1, 108
  !     write(*,*) b1(i), b12(i)
  !   end do

end subroutine initialize_bond_vector


subroutine initialize_new_direction
  !
  ! 6 probable move direction
  use global_variables
  implicit none

  new_direction(1,1) = 1
  new_direction(1,2) = 0
  new_direction(1,3) = 0

  new_direction(2,1) = 0
  new_direction(2,2) = 1
  new_direction(2,3) = 0

  new_direction(3,1) = 0
  new_direction(3,2) = 0
  new_direction(3,3) = 1  

  new_direction(4,1) = 0
  new_direction(4,2) = 0
  new_direction(4,3) = -1

  new_direction(5,1) = 0
  new_direction(5,2) = -1 
  new_direction(5,3) = 0   

  new_direction(6,1) = -1
  new_direction(6,2) = 0
  new_direction(6,3) = 0  

end subroutine initialize_new_direction


subroutine initialize_move
  !--------------------------------------!
  ! To generate the array move(108,6) whose element move(i,j) is a bond number 
  ! belonging to bonds array after the original bond, i-th bond in bonds array,
  ! with the j-th move that have 6 probable choice.
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  ! David P. Landau, Kurt Binder. A Guide to Monte Carlo Simulation in
  ! Statistical Physics, 3rd (Cambridge, UK: Cambridge University Press, 
  ! 2009), see Subroutine inimove in pp. 461.
  !Note:
  ! There are little errors in Landau's book, and are corrected here.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k
  integer, dimension(6,3) :: new
  logical :: test

  move = 0     ! move(i,j)==0, means the new bond is broken
  do i = 1, 108
    new(1,1) = bonds(i,1) + 1
    new(1,2) = bonds(i,2)
    new(1,3) = bonds(i,3)

    new(2,1) = bonds(i,1)
    new(2,2) = bonds(i,2) + 1
    new(2,3) = bonds(i,3)

    new(3,1) = bonds(i,1)
    new(3,2) = bonds(i,2)
    new(3,3) = bonds(i,3) + 1  

    new(4,1) = bonds(i,1)
    new(4,2) = bonds(i,2)
    new(4,3) = bonds(i,3) - 1

    new(5,1) = bonds(i,1)
    new(5,2) = bonds(i,2) - 1 
    new(5,3) = bonds(i,3)    

    new(6,1) = bonds(i,1) - 1
    new(6,2) = bonds(i,2)
    new(6,3) = bonds(i,3)

    do j = 1, 6
      test = .false.
      do k = 1, 108 
        test = ( new(j,1) .eq. bonds(k,1) ) .and. &
  &            ( new(j,2) .eq. bonds(k,2) ) .and. &
  &            ( new(j,3) .eq. bonds(k,3) )
        if (test) then
          move(i,j) = k ! move(i,j)==k, means the new bond 
                        ! is k bond in bonds array
          exit
        end if
      end do
    end do
  end do

  !testify the move array
  !   do i = 1, 108
  !     write(*,*) move(i,:)
  !   end do

end subroutine initialize_move


subroutine continue_read_data(l)
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
!   use compute_energy_ewald, only : charge, inv_charge
  implicit none
  integer, intent(out) :: l
  integer :: i, j 
  integer, allocatable, dimension(:,:) :: phi
  integer :: xi,yi,zi,xp,yp,zp

  allocate(phi(Lz2,9))

  open(20,file='./data/pos1.txt')
    read(20,*) ((pos(i,j),j=1,4),i=1,NN)
  close(20)
  open(199,file='./data/phi_xy.txt')
  open(200,file='./data/phi_zx.txt')
    read(199,*) ((phi_xy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(200,*) ((phi_zx(i,j),j=1,SizeHist),i=1,SizeHist)
  close(199)
  close(200)
  open(19,file='./start_time.txt')
    read(19,*) restart_or_continue
    read(19,*) l
    read(19,*) Nq_net
    read(19,*) total_time
  close(19)
  open(22,file='./data/phi.txt')
    read(22,*) ((phi(i,j),j=1,9),i=1,Lz2)
      phi_s(:,2) = phi(:,2)
      phi_sb(:,2) = phi(:,3)
      phi_se(:,2) = phi(:,4)
      phi_a(:,2) = phi(:,5)
      phi_as(:,2) = phi(:,6)
      phi_ap(:,2) = phi(:,7)
      phi_i(:,2) = phi(:,8)
      phi_is(:,2) = phi(:,9)
  close(22)

  open(23,file='./data/monbd.txt')
    read(23,*) ((monbd(i,j),j=1,arm+1),i=1,Npe)
  close(23)

  open(24,file='./data/bond_numb1.txt')
    do i = 1, N_Bond
      read(24,*) bond_numb(i)
    end do
  close(24)

  allocate(charge(Nq))

  allocate(inv_charge(NN))

  open(24,file='./data/charge.txt')
    do i = 1, Nq
      read(24,*) charge(i)
    end do
  close(24)

  open(24,file='./data/inv_charge.txt')
    do i = 1, NN
      read(24,*) inv_charge(i)
    end do
  close(24)

  do i = 1, NN
    xi = pos(i,1)
    yi = pos(i,2)
    zi = pos(i,3)
    xp = ipx(xi)
    yp = ipy(yi)
    zp = ipz(zi)
    latt(xi,yi,zi) = 1
    latt(xi,yi,zp) = 1
    latt(xi,yp,zi) = 1
    latt(xi,yp,zp) = 1
    latt(xp,yi,zi) = 1
    latt(xp,yi,zp) = 1
    latt(xp,yp,zi) = 1
    latt(xp,yp,zp) = 1
  end do

  pos0 = pos
  latt0 =latt
  bond_numb0 = bond_numb

end subroutine continue_read_data


subroutine compute_physical_quantities
  !----------------------------------------!
  !
  !input:
  !  pos
  !output:
  !  Rg, Rgz, RR2, RR2z, hh, hh_max
  !External Variables:
  !  Ngl, Nml, Npe, NN,
  !----------------------------------------!
  use global_variables
  implicit none
  integer i,j,k
  
  hh = 0
  do i = 1, Npe
    hh = hh + pos(i,3)
  end do
  hh = hh / Npe
  h_end = 0
  h_branch = 0
  do i = 1, Nga
    do j = 2, arm
      k = (i-1)*Ns+j*Nma+1
      h_end = h_end + pos(k,3)
    end do
    k = (i-1)*Ns+Nma+1
    h_branch = h_branch + pos(k,3)
  end do
  h_end = h_end/(Nga*(arm-1))
  h_branch = h_branch/Nga
  
end subroutine compute_physical_quantities


subroutine histogram
  !----------------------------------------!
  !input:
  !  pos
  !output:
  !  hist1(distribution hisotgram from PE to rod)
  !External Variants:
  !  Npe, Ngl, NN, Nml, Sizehist, Lz 
  !----------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k
  real*8, dimension(3) :: rij
  real*8 rsqr, max_h, theta, rr
  integer :: x, y, z
  
  !
  !star
  do i=1, Npe
    if ( mod(i,(arm*Nma+1))==1 .and. i<=Npe ) cycle
    k = pos(i,3)
    phi_s(k,2) = phi_s(k,2) + 1
    if (pos(i,4)/=0) then
      phi_i(k,2) = phi_i(k,2) + 1
    end if
  end do
  do i = 1, Nga
    j = pos((i-1)*(arm*Nma+1)+Nma+1,3)
    phi_sb(j,2) = phi_sb(j,2) + 1
    do j = 2, arm
      k = pos((i-1)*(arm*Nma+1)+Nma+1+(arm-1)*Nma,3)
      phi_se(k,2) = phi_se(k,2) + 1
    end do
  end do
  do i = Npe+1, Npe+Nq_PE*abs(qq)
    if (pos(i,4)/=0) then
      j = pos(i,3)
      phi_a(j,2) = phi_a(j,2) + 1
      phi_ap(j,2) = phi_ap(j,2) + 1
    end if
  end do
  do i = Npe+Nq_PE*abs(qq)+1, NN - Nq_salt_ions
    j = pos(i,3)
    phi_a(j,2) = phi_a(j,2) + 1
    phi_as(j,2) = phi_as(j,2) + 1  
  end do
  do i = NN - Nq_salt_ions+1, NN
    j = pos(i,3)
    phi_is(j,2) = phi_is(j,2) + 1
  end do
  !2d diatribution
  do i=1,Npe
    x=ceiling((pos(i,1)/2.)/(Lx/SizeHist))
    y=ceiling((pos(i,2)/2.)/(Ly/SizeHist))
    z=ceiling((pos(i,3)/2.)/(Lz/SizeHist))
    if (x==0 .or. y==0 .or. z==0) cycle
    phi_xy(x,y)=phi_xy(x,y)+1
    phi_zx(x,z)=phi_zx(x,z)+1
  end do

end subroutine histogram


subroutine write_pos
  use global_variables
  implicit none
  integer :: i, j, k

  open(102,file='./data/pos.txt')
    do i = 1, NN
      write(102,'(4I10)') (pos(i,j),j=1,4)
    end do
  close(102)

  open(103,file='./data/bond_numb.txt')
    do i = 1, N_Bond
      write(103,*) bond_numb(i)
    end do
  close(103)

  open(104,file='./data/monbd.txt')
    do i = 1, Npe
      write(104,*) monbd(i,:)
    end do
  close(104)

  open(105,file='./data/move.txt')
    do i = 1, 109
      write(105,*) move(i,:)
    end do
  close(105)  

  open(106, file='./data/ipxyz.txt')
    do i = 1, Lx2
      write(106,*) imx(i),ipx(i),ip2x(i),imy(i),ipy(i),ip2y(i),imz(i),ipz(i),ip2z(i)
    end do
  close(106)

!   open(100,file='./data/latt.txt')
!     do i = 1, Lx
!       do j = 1, Ly
!         do k = 1, Lz2+1
!           write(100,'(4I6)') i,j,k,latt(i,j,k)
!         end do
!       end do
!     end do
!   close(100)

end subroutine write_pos

subroutine write_pos1(l)
  use global_variables
!   use compute_energy_ewald, only : charge, inv_charge
  implicit none
  integer :: i, j, k
  integer, intent(in) :: l

  open(107,file='./data/pos1.txt')
    do i = 1, NN
      write(107,'(4I10)') (pos(i,j),j=1,4)
    end do
  close(107)

  open(108,file='./data/bond_numb1.txt')
    do i = 1, N_Bond
      write(108,*) bond_numb(i)
    end do
  close(108)

  open(108,file='./data/charge.txt')
    do i = 1, Nq
      write(108,*) charge(i)
    end do
  close(108)  

  open(108,file='./data/inv_charge.txt')
    do i = 1, NN
      write(108,*) inv_charge(i)
    end do
  close(108) 

  open(109,file='./start_time.txt')
    write(109,*) 1
    write(109,*) l
    write(109,*) Nq_net
    call cpu_time(finished)
    total_time=total_time+finished-started
    call cpu_time(started)
    write(109,*) total_time
    write(109,*) 'time:(minutes):', real(total_time/60)
    write(109,*) 'time:(hours)  :', real(total_time/3600)
    write(109,*) 'time:(days)   :', real(total_time/86400)
    write(109,*) 'Lx2           :', real(Lx2)
    write(109,*) 'Ly2           :', real(Ly2)
    write(109,*) 'Lz2           :', real(Lz2)
    write(109,*) 'Lx            :', Lx
    write(109,*) 'Ly            :', Ly
    write(109,*) 'Lz            :', Lz
    write(109,*) 'Nga           :', Nga
    write(109,*) 'Nma           :', Nma
    write(109,*) 'Nq            :', Nq
    write(109,*) 'NN            :', NN
    write(109,*) 'Z_empty       :', Z_empty-1
    write(109,*) 'sigmag        :', sigmag
    write(109,*) 'true sigmag   :', sigmag1
    write(109,*) 'Beta          :', Beta
    write(109,*) 'qq            :', qq
    write(109,*) 'qqi           :', qqi
    write(109,*) 'Arm           :', arm
    write(109,*) 'Nma           :', Nma
    write(109,*) 'Nga           :', Nga
    write(109,*) 'man_s         :', man_s
    write(109,*) 'NN            :', NN
    write(109,*) 'Nq            :', Nq
    write(109,*) 'Nq_PE         :', Nq_PE
    write(109,*) 'Npe           :', Npe
    write(109,*) 'Nq_salt_ions  :', Nq_salt_ions
    write(109,*) 'Lx2           :', Lx2
    write(109,*) 'Ly2           :', Ly2
    write(109,*) 'Lz2           :', Lz2
    write(109,*) 'N_bond        :', N_bond
    write(109,*) 'pH-pKa        :', pH_pKa
    write(109,*) 'restart_continue:',restart_or_continue
    write(109,*) 'StepNum0      :', StepNum0
    write(109,*) 'StepNum       :', StepNum
    write(109,*) 'StepNum0+StepNum:', (StepNum0+StepNum)
    write(109,*) 'DeltaStep     :', DeltaStep
    write(109,*) 'DeltaStep1    :', DeltaStep1
    write(109,*) 'DeltaStep2    :', DeltaStep2
    write(109,*) 'DeltaStep3    :', DeltaStep3
  close(109)

!   open(100,file='./data/latt1.txt')
!     do i = 1, Lx2
!       do j = 1, Ly2
!         do k = 1, Lz2+1
!           write(100,'(4I6)') i,j,k,latt(i,j,k)
!         end do
!       end do
!     end do
!   close(100)

end subroutine write_pos1


subroutine write_hist
  !----------------------------------------!
  !Write distribution histogram to the file hist1.txt ... hist4.txt
  !input:
  !  hist1, hist2, hist3, hist4,...hist9
  !External Variants: 
  !  SizeHist, Lz, 
  !----------------------------------------!
  use global_variables
  implicit none
  integer i,j
  
  open(34,file='./data/phi.txt')
    do i=1,Lz2
      write(34,340) i, phi_s(i,2), phi_sb(i,2), phi_se(i,2), phi_a(i,2),&
       phi_as(i,2), phi_ap(i,2), phi_i(i,2), phi_is(i,2)
    end do
    340 format(9I10)
  close(34)
  open(199,file='./data/phi_xy.txt')
  open(200,file='./data/phi_zx.txt')
    do i=1,SizeHist
      write(199,'(500I10)') (phi_xy(i,j),j=1,SizeHist)
      write(200,'(500I10)') (phi_zx(i,j),j=1,SizeHist)
    end do
  close(199)
  close(200)

end subroutine write_hist


subroutine write_physical_quantities(j)
  !----------------------------------------!
  !
  !----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: j
  
  open(110,position='append', file='./data/height.txt')
    write(110,1100) 1.*j, hh, h_end, h_branch
    1100 format(4F15.6)
  close(110)

end subroutine write_physical_quantities


subroutine write_energy(j,EE,EE1)
  use global_variables
  implicit none
  integer, intent(in) :: j
  real*8, intent(in) :: EE
  real*8, intent(in) :: EE1

  open(36,position='append', file='./data/energy.txt')
    write(36,361) 1.*j, EE, EE1, Nq_net*1.D0/Nq_PE*100, &
           NN-(Npe/man_s-Nq_net)*2.D0, accept_ratio, &
           accept_ratio_pH, accept_ratio_long, rmse*100
    361 format(9F17.5)
  close(36)  

end subroutine write_energy


end module input_output


