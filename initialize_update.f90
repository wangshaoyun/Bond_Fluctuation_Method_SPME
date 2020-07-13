module initialize_update
implicit none
! compute_energy
! Delta_energy
! Delta_energy_add
! Delta_energy_delete

contains

subroutine initialize_position
  use global_variables
  implicit none
  integer :: i

  !
  !Graft star brushes, the arms are all vertical to the plate.
  !The start point on plate can be uniform of random
  call uniform_star_brushes
  !
  !random distribution in the box
  if ( qq /= 0 ) then
    call initialize_ions
  end if

  pos0 = pos
  latt0 = latt
  bond_numb0 = bond_numb

end subroutine initialize_position


subroutine uniform_star_brushes
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n
  integer :: nx, dLx
  integer :: base, base1, base2
  integer :: xi, yi, zi, xp, yp, zp
  integer, dimension(6,4) :: bond_vector

  nx = nint( sqrt(1.*Nga) )
  dLx = Lx2 / nx
  !up, bond_vector is bond number in bonds array by looking up bonds.txt
  bond_vector(1,1) = 0
  bond_vector(1,2) = 0
  bond_vector(1,3) = 2
  bond_vector(1,4) = 3     !bond number in bonds.txt
  !right
  bond_vector(2,1) = 2
  bond_vector(2,2) = 0
  bond_vector(2,3) = 0 
  bond_vector(2,4) = 1     
  !forward
  bond_vector(3,1) = 0
  bond_vector(3,2) = 2
  bond_vector(3,3) = 0 
  bond_vector(3,4) = 5     
  !left
  bond_vector(4,1) = -2
  bond_vector(4,2) = 0
  bond_vector(4,3) = 0   
  bond_vector(4,4) = 2     
  !backward
  bond_vector(5,1) = 0
  bond_vector(5,2) = -2
  bond_vector(5,3) = 0
  bond_vector(5,4) = 6    
  !down 
  bond_vector(6,1) = 0
  bond_vector(6,2) = 0
  bond_vector(6,3) = -2
  bond_vector(6,4) = 4     !bond number in bonds.txt
  do i = 1, nx
    do j = 1, nx
      l = (i-1)*nx + j
      base = (l-1)*(arm*Nma+1)
      pos(base+1, 1) = (i-1)*dLx + 3
      pos(base+1, 2) = (j-1)*dLx + 3
      pos(base+1, 3) = 1
      monbd(base+1,1) = base + 2 - l    
      monbd(base+1,arm+1) = 1
      do k = 2, Nma + 1
        pos(base+k, 1:3) = pos(base+k-1,1:3) + bond_vector(1,1:3)
        bond_numb(base+k-l) = bond_vector(1,4)
        monbd(base+k,1) = base + k - l       
        monbd(base+k,2) = base + k + 1 - l 
        monbd(base+k,arm+1) = 2
      end do
      do m = 2, arm
        base1 = base + Nma + 1      
        base2 = base1 + (m-2) * Nma 
        pos(base2+1, 1:3) = pos(base1, 1:3) + bond_vector(m-1,1:3)
        if (pos(base2+1, 1)>Lx2) then
          pos(base2+1,1) = pos(base2+1,1)-Lx2
        end if
        if (pos(base2+1, 2)>Ly2) then
          pos(base2+1,2) = pos(base2+1,2)-Ly2
        end if
        bond_numb(base2+1-l) = bond_vector(m,4)
        monbd(base1,m) = base2 + 1 - l     
        monbd(base1,arm+1) = arm
        monbd(base2+1,1) = base2 + 1 - l   
        monbd(base2+1,2) = base2 + 2 - l    
        monbd(base2+1,arm+1) = 2
        do n = 2, Nma
          if ( m == 2 ) then
            pos(base2+n,1:3) = pos(base2+n-1,1:3) + bond_vector(1,1:3)
            bond_numb(base2+n-l) = bond_vector(1,4)
          else
            pos(base2+n,1:3) = pos(base2+n-1,1:3) + bond_vector(1,1:3)
            bond_numb(base2+n-l) = bond_vector(1,4)
          end if
          monbd(base2+n,1) = base2 + n - l                
          monbd(base2+n,2) = base2 + n + 1 - l            
          monbd(base2+n,arm+1) = 2
        end do
        monbd(base2+Nma,2) = 0               
        monbd(base2+Nma,arm+1) = 1
      end do 
    end do
  end do

  !the vortexes of the lattice with particle are occupied
  do j = 1, Nga
    k = 1
    do l = 1, Nma + 1
      i = (j-1)*Ns + l
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
      if ( mod(l-1,man_s)==0 .and. l/=1 ) then
        pos(i,4) = qq
      end if
    end do
    do k = 2, arm
      do l = 1, Nma
        i = (j-1)*Ns + (Nma + 1) + ((k-2) * Nma) + l
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
        if ( mod(l,man_s)==0 ) then
          pos(i,4) = qq
        end if
      end do
    end do
  end do

  !the above and bottom plates are occupied
  do i = 1, Lx2
    do j = 1, Ly2
      latt(i,j,1) = 1
      latt(i,j,Lz2+1) = 1
    end do
  end do

end subroutine uniform_star_brushes


subroutine initialize_ions
  use global_variables
  implicit none
  integer :: i, j, k
  integer :: total
  integer :: xi, yi, zi, xp, yp, zp
  real*8, dimension(3) :: rnd
  logical :: test

  do i = Npe+1, NN
    test=.true.
    do while (test)
      call random_number(rnd)
      xi = floor(rnd(1)*Lx2) + 1
      yi = floor(rnd(2)*Ly2) + 1
      zi = floor(rnd(3)*Lz2) + 1
      xp = ipx(xi)
      yp = ipy(yi)
      zp = ipz(zi)
      if (zp==0) then
        write(*,*) zp
      end if
      total = latt(xi,yi,zi)+latt(xi,yi,zp)+latt(xi,yp,zi)+latt(xi,yp,zp) +   &
              latt(xp,yi,zi)+latt(xp,yi,zp)+latt(xp,yp,zi)+latt(xp,yp,zp)
      if ( total == 0 ) then
        pos(i,1) = xi
        pos(i,2) = yi
        pos(i,3) = zi
        latt(xi,yi,zi) = 1
        latt(xi,yi,zp) = 1
        latt(xi,yp,zi) = 1
        latt(xi,yp,zp) = 1
        latt(xp,yi,zi) = 1
        latt(xp,yi,zp) = 1
        latt(xp,yp,zi) = 1
        latt(xp,yp,zp) = 1        
        test = .false.
      end if
    end do
    if (i<=Npe+Nq_PE*abs(qq)) then
      pos(i,4) = - qq !/ abs(qq)
    elseif (i>NN-Nq_salt_ions) then
      pos(i,4) = qqi
    else
      pos(i,4) = - qqi / abs(qqi)
    end if
  end do

end subroutine initialize_ions


subroutine error_analysis(n, EE)
  use global_variables
  use compute_energy_ewald
  implicit none
  integer, intent(in) :: n
  real*8, intent(out) :: EE

  if (n==0) then
    call error_analysis_ewald(EE)
  else
    call total_energy_spme(EE)
  end if

!   write(*,*) 'initial energy',EE

end subroutine error_analysis


subroutine monte_carlo_move( EE, DeltaE )
  use global_variables
  use compute_energy_ewald
  implicit none
  real*8, intent(inout) :: EE
  real*8, intent(out) :: DeltaE
  integer :: i
  real*8 :: st,fn
  real*8 :: real_time, fourier_time, EE1, EE2

  call cpu_time(st)
    do i = 1, NN - Nga
!       if ( mod(i,DeltaStep) == 0 .and. qq/=0 .and. pH_or_not/=0 ) then
      if ( mod(i,(NN - Nga)) == 0 .and. qq/=0 .and. pH_or_not/=0 ) then
        num_pH = num_pH + 1
        call choose_particle_pH
        if (pos(ip,4)==0) then
          call add_particle(EE,DeltaE)
        else
          call delete_particle(EE,DeltaE)
        end if
      else
        call choose_particle
        call new_position(EE,DeltaE)
        num_move = num_move + 1
      end if
    end do
  call cpu_time(fn)

end subroutine monte_carlo_move


subroutine update_multistep(EE,EE2)
  use global_variables
  use compute_energy_ewald
  implicit none
  real*8 :: rnd, DeltaE
  real*8, intent(inout) :: EE
  real*8, intent(inout) :: EE2
  
  call compute_Nq_net

  if (Nq_net /=0 ) then
    call SPME_Ewald(energy_long1)
  end if

  call random_number(rnd)

  DeltaE = energy_long1 - energy_long0
!   write(*,*) DeltaE
  EE = EE + DeltaE
  num_long = num_long + 1
  total_accept = total_accept + accept
  total_accept_pH = total_accept_pH + accept_pH
  total_accept_long = total_accept_long + 1

  if ( DeltaE > 0 ) then
    if (rnd>exp(-beta*DeltaE)) then
      pos = pos0
      latt = latt0
      bond_numb = bond_numb0
      EE = EE2
      energy_long1 = energy_long0
      call compute_Nq_net
      call initialize_energy_arrays_ewald
      total_accept = total_accept - accept
      total_accept_pH = total_accept_pH - accept_pH
      total_accept_long = total_accept_long - 1
    end if
  end if
  pos0 = pos
  latt0 = latt
  bond_numb0 = bond_numb
  EE2 = EE
  energy_long0 = energy_long1
  accept = 0
  accept_pH = 0
  accept_ratio = 1.*total_accept / num_move
  accept_ratio_pH = 1.*total_accept_pH / num_pH
  accept_ratio_long = 1.*total_accept_long / num_long

!   write(*,*) EE

!   total_accept = 0
!   total_accept_long = 0
!   total_accept_pH = 0
!   num_move = 0
!   num_long = 0
!   num_pH = 0
!   write(*,*) accept_ratio,accept_ratio_long,accept_ratio_pH

end subroutine update_multistep


subroutine choose_particle
  use global_variables
  implicit none
  real*8 :: rnd

  call random_number(rnd)
  ip = floor(rnd*NN) + 1

  !
  !The monomer anchored on the plate can't move, so we need to choose again.
  do while( (mod(ip,Ns)==1 .and. ip<=Npe) .or. (ip>Npe .and. pos(ip,4)==0) )
    call random_number(rnd)
    ip = int(rnd*NN) + 1
  end do

end subroutine choose_particle


subroutine choose_particle_pH
  use global_variables
  use compute_energy_ewald
  implicit none
  real*8 :: rnd
  integer :: i

  call random_number( rnd )
  i = floor(rnd*Nq_PE) + 1
  ip = charge(i)
  ip1 = Npe + i

end subroutine choose_particle_pH


subroutine add_particle(EE,DeltaE)
  use global_variables
  use compute_energy_ewald
  implicit none
  real*8, intent(inout) :: EE
  real*8, intent(out) :: DeltaE
  real*8 :: rnd(3), U_prot, rnd1, Delta_Prob
  integer :: xi,yi,zi,xp,yp,zp,total

  U_prot = log(10.D0)/beta*pH_pKa !+: add, -:delete

  pos_ip0 = pos(ip,:)         !old polymer
  pos_ip0i = pos(ip1,:)       !old ions

  !
  !new position
  call random_number(rnd)
  pos_ip1i(1) = floor(rnd(1)*Lx2)+1  !new ions
  pos_ip1i(2) = floor(rnd(2)*Ly2)+1
  pos_ip1i(3) = floor(rnd(3)*Lz2)+1
  pos_ip1i(4) = -qq/abs(qq)
  pos_ip1 = pos_ip0                  !new polymer
  pos_ip1(4) = qq

  !
  !judge the excluded volume condition
  xi = pos_ip1i(1)
  yi = pos_ip1i(2)
  zi = pos_ip1i(3)
  xp = ipx(xi)
  yp = ipy(yi)
  zp = ipz(zi)
  total = latt(xi,yi,zi)+latt(xi,yi,zp)+latt(xi,yp,zi)+latt(xi,yp,zp) +   &
          latt(xp,yi,zi)+latt(xp,yi,zp)+latt(xp,yp,zi)+latt(xp,yp,zp)
  if (total == 0) then
    call Delta_Energy_Ewald_add(DeltaE)
    if (Nq_net/=0) then
      Delta_Prob = 1. * ( Nq_PE - Nq_net ) / Nq_net * exp(-(DeltaE-U_prot)*beta)
    else
      Delta_Prob = 1
    end if
!     write(*,*) DeltaE
    if (Delta_Prob>=1) then
      latt(xi,yi,zi) = 1
      latt(xi,yi,zp) = 1
      latt(xi,yp,zi) = 1
      latt(xi,yp,zp) = 1
      latt(xp,yi,zi) = 1
      latt(xp,yi,zp) = 1
      latt(xp,yp,zi) = 1
      latt(xp,yp,zp) = 1
      pos(ip,:) = pos_ip1
      pos(ip1,:) = pos_ip1i
      call update_cell_list_pH_add_Ewald
      EE = EE + DeltaE
      accept_pH = accept_pH + 1
    else
      call random_number(rnd1)
      if (rnd1<Delta_Prob) then
        latt(xi,yi,zi) = 1
        latt(xi,yi,zp) = 1
        latt(xi,yp,zi) = 1
        latt(xi,yp,zp) = 1
        latt(xp,yi,zi) = 1
        latt(xp,yi,zp) = 1
        latt(xp,yp,zi) = 1
        latt(xp,yp,zp) = 1
        pos(ip,:) = pos_ip1
        pos(ip1,:) = pos_ip1i
        call update_cell_list_pH_add_Ewald
        EE = EE + DeltaE 
        accept_pH = accept_pH + 1
      end if 
    end if
  end if

end subroutine add_particle


subroutine delete_particle(EE, DeltaE)
  use global_variables
  use compute_energy_ewald
  implicit none
  real*8, intent(inout) :: EE
  real*8, intent(out) :: DeltaE
  real*8 :: U_prot, rnd, Delta_Prob
  integer :: xi, yi, zi, xp, yp, zp

  U_prot = log(10.D0)/beta*pH_pKa !+: add, -:delete

  pos_ip0 = pos(ip,:)             !polymer
  pos_ip1 = pos_ip0
  pos_ip1(4) = 0
  pos_ip0i = pos(ip1,:)           !ions
  pos_ip1i = pos_ip0i
  pos_ip1i(4) = 0
  xi = pos_ip0i(1)
  yi = pos_ip0i(2)
  zi = pos_ip0i(3)
  xp = ipx(xi)
  yp = ipy(yi)
  zp = ipz(zi)

  call Delta_Energy_Ewald_delete(DeltaE)
  if ((Nq_PE - Nq_net)/=0) then
    Delta_Prob = 1. * Nq_net / ( Nq_PE - Nq_net ) * exp(-(DeltaE+U_prot)*beta)
  else
    Delta_Prob = 1
  end if
!   write(*,*) DeltaE
  if (Delta_Prob>=1) then
    latt(xi,yi,zi) = 0
    latt(xi,yi,zp) = 0
    latt(xi,yp,zi) = 0
    latt(xi,yp,zp) = 0
    latt(xp,yi,zi) = 0
    latt(xp,yi,zp) = 0
    latt(xp,yp,zi) = 0
    latt(xp,yp,zp) = 0    
    pos(ip,:) = pos_ip1
    pos(ip1,:) = pos_ip1i
    call update_cell_list_pH_delete_Ewald
    EE = EE + DeltaE
    accept_pH = accept_pH + 1
  else
    call random_number(rnd)
    if (rnd<Delta_Prob) then
      latt(xi,yi,zi) = 0
      latt(xi,yi,zp) = 0
      latt(xi,yp,zi) = 0
      latt(xi,yp,zp) = 0
      latt(xp,yi,zi) = 0
      latt(xp,yi,zp) = 0
      latt(xp,yp,zi) = 0
      latt(xp,yp,zp) = 0    
      pos(ip,:) = pos_ip1
      pos(ip1,:) = pos_ip1i
      call update_cell_list_pH_delete_Ewald
      EE = EE + DeltaE
      accept_pH = accept_pH + 1
    end if
  end if

end subroutine delete_particle


subroutine new_position(EE, DeltaE)
  use global_variables
  use compute_energy_ewald
  implicit none
  real*8,  intent(out)   :: DeltaE
  real*8,  intent(inout) :: EE
  integer :: dir  ! direction, with 6 choice
  integer :: bn   ! number of bonds connect to the particle
  integer, allocatable, dimension(:) :: new_bonds
  integer :: i, ix, iy, iz, testlat, summ
  integer :: xm1, xp1, xp2, ym1, yp1, yp2, zm1, zp1, zp2 
  real*8  :: rnd
  logical :: test

  call random_number(rnd)
  dir = floor(6*rnd)+1

  pos_ip0 = pos(ip,:)
  pos_ip1(1:3) = pos_ip0(1:3) + new_direction(dir,:)
  pos_ip1(4)   = pos_ip0(4)
  call periodic_condition( pos_ip1(1:2) )

  ix = pos_ip0(1)
  iy = pos_ip0(2)
  iz = pos_ip0(3)

  test = .false.
  if ( ip <= Npe ) then
    bn = monbd(ip,arm+1)
    if ( allocated(new_bonds) ) deallocate(new_bonds)
    allocate(new_bonds(bn))
    new_bonds(1) = move( bond_numb( monbd( ip, 1 ) ), dir )
    do i = 2, bn
      new_bonds(i) = move( bond_numb( monbd( ip, i ) ), 7-dir )
    end do
    test = .true.
    do i = 1, bn
      if ( new_bonds(i) == 0 ) then
        test = .false.
        exit
      end if
    end do
  end if

  if ( test .or. (ip>Npe) )then
    select case (dir)
    case ( 1 )
      xp2 = ip2x(ix)
      yp1 = ipy(iy)
      zp1 = ipz(iz)
      testlat = latt(xp2,iy,iz)  + latt(xp2,yp1,iz) + &
                latt(xp2,iy,zp1) + latt(xp2,yp1,zp1)
      if( testlat == 0 ) then
        if (pos_ip1(4)/=0) then
          call Delta_Energy_Ewald(DeltaE)
        else
          DeltaE = 0
        end if
        if ( DeltaE <= 0 ) then
          pos(ip,1:3) = pos_ip1(1:3)
          EE = EE + DeltaE
          if ( ip <= Npe ) then
            do i = 1, bn
              bond_numb( monbd( ip, i ) )= new_bonds(i)
            end do
          end if
          latt(xp2,iy,iz)   = 1 
          latt(xp2,yp1,iz)  = 1 
          latt(xp2,iy,zp1)  = 1 
          latt(xp2,yp1,zp1) = 1 
          latt(ix,iy,iz)    = 0 
          latt(ix,yp1,iz)   = 0 
          latt(ix,iy,zp1)   = 0 
          latt(ix,yp1,zp1)  = 0 
          if (pos_ip1(4)/=0) then
            call update_real_cell_list_Ewald
          end if
          accept = accept + 1
        else
          call random_number(rnd)
          if ( rnd < Exp(-Beta*DeltaE) ) then
            pos(ip,1:3) = pos_ip1(1:3) 
            EE = EE + DeltaE
            if ( ip <= Npe ) then
              do i = 1, bn
                bond_numb( monbd( ip, i ) )= new_bonds(i)
              end do
            end if
            latt(xp2,iy,iz)   = 1 
            latt(xp2,yp1,iz)  = 1 
            latt(xp2,iy,zp1)  = 1 
            latt(xp2,yp1,zp1) = 1 
            latt(ix,iy,iz)    = 0 
            latt(ix,yp1,iz)   = 0 
            latt(ix,iy,zp1)   = 0 
            latt(ix,yp1,zp1)  = 0 
            if (pos_ip1(4)/=0  ) then
              call update_real_cell_list_Ewald
            end if
            accept = accept + 1
          end if
        end if
      end if
    case ( 6 ) 
      xm1 = imx(ix)
      xp1 = ipx(ix)
      yp1 = ipy(iy)
      zp1 = ipz(iz)
      testlat = latt(xm1,iy,iz)  + latt(xm1,yp1,iz) + &
         &      latt(xm1,iy,zp1) + latt(xm1,yp1,zp1)
      if( testlat == 0 ) then
        if (pos_ip1(4)/=0 ) then
          call Delta_Energy_Ewald(DeltaE)
        else
          DeltaE = 0
        end if
        if ( DeltaE <= 0 ) then
          pos(ip,1:3) = pos_ip1(1:3)
          EE = EE + DeltaE
          if ( ip <= Npe ) then
            do i = 1, bn
              bond_numb( monbd( ip, i ) )= new_bonds(i)
            end do
          end if
          latt(xm1,iy,iz)   = 1 
          latt(xm1,yp1,iz)  = 1 
          latt(xm1,iy,zp1)  = 1 
          latt(xm1,yp1,zp1) = 1 
          latt(xp1,iy,iz)   = 0 
          latt(xp1,yp1,iz)  = 0 
          latt(xp1,iy,zp1)  = 0
          latt(xp1,yp1,zp1) = 0 
          if (pos_ip1(4)/=0  ) then
            call update_real_cell_list_Ewald
          end if
          accept = accept + 1
        else
          call random_number(rnd)
          if ( rnd < Exp(-Beta*DeltaE) ) then
            pos(ip,1:3) = pos_ip1(1:3) 
            EE = EE + DeltaE
            if ( ip <= Npe ) then
              do i = 1, bn
                bond_numb( monbd( ip, i ) )= new_bonds(i)
              end do
            end if
            latt(xm1,iy,iz)   = 1 
            latt(xm1,yp1,iz)  = 1 
            latt(xm1,iy,zp1)  = 1 
            latt(xm1,yp1,zp1) = 1 
            latt(xp1,iy,iz)   = 0 
            latt(xp1,yp1,iz)  = 0 
            latt(xp1,iy,zp1)  = 0 
            latt(xp1,yp1,zp1) = 0
            if (pos_ip1(4)/=0 ) then
              call update_real_cell_list_Ewald
            end if
            accept = accept + 1
          end if
        end if
      end if
    case (2)
      xp1 = ipx(ix)
      yp2 = ip2y(iy)
      zp1 = ipz(iz)
      testlat = latt(ix,yp2,iz)  + latt(xp1,yp2,iz) + &
        &       latt(ix,yp2,zp1) + latt(xp1,yp2,zp1)
      if( testlat == 0 ) then
        if (pos_ip1(4)/=0 ) then
          call Delta_Energy_Ewald(DeltaE)
        else
          DeltaE = 0
        end if
        if ( DeltaE < 0 ) then
          pos(ip,1:3) = pos_ip1(1:3)
          EE = EE + DeltaE
          if ( ip <= Npe ) then
            do i = 1, bn
              bond_numb( monbd( ip, i ) ) = new_bonds(i)
            end do
          end if
          latt(ix,yp2,iz)   = 1 
          latt(xp1,yp2,iz)  = 1 
          latt(ix,yp2,zp1)  = 1 
          latt(xp1,yp2,zp1) = 1 
          latt(ix,iy,iz)    = 0 
          latt(xp1,iy,iz)   = 0 
          latt(ix,iy,zp1)   = 0 
          latt(xp1,iy,zp1)  = 0 
          if (pos_ip1(4)/=0  ) then
            call update_real_cell_list_Ewald
          end if
          accept = accept + 1
        else
          call random_number(rnd)
          if ( rnd < Exp(-Beta*DeltaE) ) then
            pos(ip,1:3) = pos_ip1(1:3) 
            EE = EE + DeltaE
            if ( ip <= Npe ) then
              do i = 1, bn
                bond_numb( monbd( ip, i ) ) = new_bonds(i)
              end do
            end if
            latt(ix,yp2,iz)   = 1 
            latt(xp1,yp2,iz)  = 1 
            latt(ix,yp2,zp1)  = 1 
            latt(xp1,yp2,zp1) = 1 
            latt(ix,iy,iz)    = 0 
            latt(xp1,iy,iz)   = 0 
            latt(ix,iy,zp1)   = 0 
            latt(xp1,iy,zp1)  = 0 
            if (pos_ip1(4)/=0  ) then
              call update_real_cell_list_Ewald
            end if
            accept = accept + 1
          end if
        end if
      end if
    case (5)
      xp1 = ipx(ix)
      ym1 = imy(iy)
      yp1 = ipy(iy)
      zp1 = ipz(iz)
      testlat = latt(ix,ym1,iz)  + latt(xp1,ym1,iz) + &
       &        latt(ix,ym1,zp1) + latt(xp1,ym1,zp1)
      if( testlat == 0 ) then
        if (pos_ip1(4)/=0 ) then
          call Delta_Energy_Ewald(DeltaE)
        else
          DeltaE = 0
        end if
        if ( DeltaE < 0 ) then
          pos(ip,1:3) = pos_ip1(1:3)
          EE = EE + DeltaE
          if ( ip <= Npe ) then
            do i = 1, bn
              bond_numb( monbd( ip, i ) ) = new_bonds(i)
            end do
          end if
          latt(ix,ym1,iz)   = 1 
          latt(xp1,ym1,iz)  = 1 
          latt(ix,ym1,zp1)  = 1 
          latt(xp1,ym1,zp1) = 1 
          latt(ix,yp1,iz)   = 0 
          latt(xp1,yp1,iz)  = 0 
          latt(ix,yp1,zp1)  = 0 
          latt(xp1,yp1,zp1) = 0 
          if (pos_ip1(4)/=0  ) then
            call update_real_cell_list_Ewald
          end if
          accept = accept + 1
        else
          call random_number(rnd)
          if ( rnd < Exp(-Beta*DeltaE) ) then
            pos(ip,1:3) = pos_ip1(1:3) 
            EE = EE + DeltaE
            if ( ip <= Npe ) then
              do i = 1, bn
                bond_numb( monbd( ip, i ) ) = new_bonds(i)
              end do
            end if
            latt(ix,ym1,iz)   = 1 
            latt(xp1,ym1,iz)  = 1 
            latt(ix,ym1,zp1)  = 1 
            latt(xp1,ym1,zp1) = 1 
            latt(ix,yp1,iz)   = 0 
            latt(xp1,yp1,iz)  = 0 
            latt(ix,yp1,zp1)  = 0 
            latt(xp1,yp1,zp1) = 0 
            if (pos_ip1(4)/=0  ) then
              call update_real_cell_list_Ewald
            end if
            accept = accept + 1
          end if
        end if
      end if
    case (3)
      xp1 = ipx(ix)
      yp1 = ipy(iy)
      zp2 = ip2z(iz)
      testlat = latt(ix,iy,zp2)  + latt(xp1,iy,zp2) + &
         &      latt(ix,yp1,zp2) + latt(xp1,yp1,zp2)
      if( testlat == 0 ) then
        if (pos_ip1(4)/=0) then
          call Delta_Energy_Ewald(DeltaE)
        else
          DeltaE = 0
        end if
        if ( DeltaE < 0 ) then
          pos(ip,1:3) = pos_ip1(1:3)
          EE = EE + DeltaE
          if ( ip <= Npe ) then
            do i = 1, bn
              bond_numb( monbd( ip, i ) ) = new_bonds(i)
            end do
          end if
          latt(ix,iy,zp2)   = 1 
          latt(xp1,iy,zp2)  = 1 
          latt(ix,yp1,zp2)  = 1 
          latt(xp1,yp1,zp2) = 1 
          latt(ix,iy,iz)    = 0 
          latt(xp1,iy,iz)   = 0 
          latt(ix,yp1,iz)   = 0 
          latt(xp1,yp1,iz)  = 0 
          if (pos_ip1(4)/=0   ) then
            call update_real_cell_list_Ewald
          end if
          accept = accept + 1
        else
          call random_number(rnd)
          if ( rnd < Exp(-Beta*DeltaE) ) then
            pos(ip,1:3) = pos_ip1(1:3) 
            EE = EE + DeltaE
            if ( ip <= Npe ) then
              do i = 1, bn
                bond_numb( monbd( ip, i ) ) = new_bonds(i)
              end do
            end if
            latt(ix,iy,zp2)   = 1 
            latt(xp1,iy,zp2)  = 1 
            latt(ix,yp1,zp2)  = 1 
            latt(xp1,yp1,zp2) = 1 
            latt(ix,iy,iz)    = 0 
            latt(xp1,iy,iz)   = 0 
            latt(ix,yp1,iz)   = 0 
            latt(xp1,yp1,iz)  = 0 
            if (pos_ip1(4)/=0   ) then
              call update_real_cell_list_Ewald
            end if
            accept = accept + 1
          end if
        end if
      end if
    case (4)
      xp1 = ipx(ix)
      yp1 = ipy(iy)
      zm1 = imz(iz)
      zp1 = ipz(iz)
      testlat = latt(ix,iy,zm1)  + latt(xp1,iy,zm1) + &
       &        latt(ix,yp1,zm1) + latt(xp1,yp1,zm1)
      if( testlat == 0 ) then
        if (pos_ip1(4)/=0) then
          call Delta_Energy_Ewald(DeltaE)
        else
          DeltaE = 0
        end if
        if ( DeltaE < 0 ) then
          pos(ip,1:3) = pos_ip1(1:3)
          EE = EE + DeltaE
          if ( ip <= Npe ) then
            do i = 1, bn
              bond_numb( monbd( ip, i ) ) = new_bonds(i)
            end do
          end if
          latt(ix,iy,zm1)   = 1 
          latt(xp1,iy,zm1)  = 1 
          latt(ix,yp1,zm1)  = 1 
          latt(xp1,yp1,zm1) = 1 
          latt(ix,iy,zp1)   = 0 
          latt(xp1,iy,zp1)  = 0 
          latt(ix,yp1,zp1)  = 0 
          latt(xp1,yp1,zp1) = 0 
          if (pos_ip1(4)/=0   ) then
            call update_real_cell_list_Ewald
          end if
          accept = accept + 1
        else
          call random_number(rnd)
          if ( rnd < Exp(-Beta*DeltaE) ) then
            pos(ip,1:3) = pos_ip1(1:3) 
            EE = EE + DeltaE
            if ( ip <= Npe ) then
              do i = 1, bn
                bond_numb( monbd( ip, i ) ) = new_bonds(i)
              end do
            end if
            latt(ix,iy,zm1)   = 1 
            latt(xp1,iy,zm1)  = 1 
            latt(ix,yp1,zm1)  = 1 
            latt(xp1,yp1,zm1) = 1 
            latt(ix,iy,zp1)   = 0 
            latt(xp1,iy,zp1)  = 0 
            latt(ix,yp1,zp1)  = 0 
            latt(xp1,yp1,zp1) = 0 
            if (pos_ip1(4)/=0   ) then
              call update_real_cell_list_Ewald
            end if
            accept = accept + 1
          end if
        end if
      end if
    end select
!     summ = 0
!     do i = 1, 401
!       summ = summ + sum(latt(:,:,i))
!     end do
!     write(*,*) summ
  end if

end subroutine new_position


end module initialize_update





