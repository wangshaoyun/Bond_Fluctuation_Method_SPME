module compute_energy_ewald
  !--------------------------------------!
  ! 
  !--------------------------------------!
  implicit none

  save

!############coefficient in potential function#############!
  !
  !coulomb
  !
  !coulomb potential
  real*8,  private :: lb          !Bjerrum length
  real*8,  private :: EF          !electric field
  real*8,  private :: tol         !tolerance
  real*8,  private :: tau_rf      !time ratio of real space and fourier space
  real*8,  private :: alpha       !Ewald screening parameter alpha
  real*8,  private :: alpha2      !alpha2=alpha*alpha
  real*8,  private :: Mz          !total dipole moment
  real*8,  private :: dMz         !delta Mz
  real*8,  private :: Mz_coef     !Coefficients in correction coulomb energy
  real*8,  private :: err
  !
  !real space
  real*8,  private :: rcc         !Cut off radius of real space
  real*8,  private :: rcc2        !rcc2=rcc*rcc  !
  real*8,  private :: clx         !length of cell in x direction
  real*8,  private :: cly         !length of cell in y direction
  real*8,  private :: clz         !length of cell in z direction
  integer, private :: nclx        !number of cell in x direction
  integer, private :: ncly        !number of cell in y direction
  integer, private :: nclz        !number of cell in z direction 
  !reciprocal space
  integer, private :: Kmax1       !max wave number of x direction
  integer, private :: Kmax2       !max wave number of y direction 
  integer, private :: Kmax3       !max wave number of z direction 
  integer, private :: K_total     !Total wave number in reciprocal space
  integer, dimension(3) :: ordr   !Order of spline function
  integer, dimension(3) :: ng     !wave number
  real*8,  dimension(3) :: gdim   !length, width and height of the box
!##########end coefficient in potential function###########!


!##########################arrays##########################!
  !
  !
  real, allocatable, dimension( :, : )            :: posq
  !
  !neighbor cells of the center cell
  integer, allocatable, dimension(:,:,:), private :: cell_near_list 
  !
  !cell list of charge
  integer, allocatable, dimension( : ), private :: cell_list_q
  !
  !inverse cell list of charge
  integer, allocatable, dimension( : ), private :: inv_cell_list_q
  !
  !cell list in real space
  integer, allocatable, dimension( : ), private :: cell_list_r
  !
  !inverse cell list in real space
  integer, allocatable, dimension( : ), private :: inv_cell_list_r
  !
  ! head of chains, cell list
  integer, allocatable, dimension(:,:,:), private :: hoc_r     
  !
  ! head of chains, inverse cell list
  integer, allocatable, dimension(:,:,:), private :: inv_hoc_r
  !
  ! Periodic condition
  integer, allocatable, dimension( : ), private :: periodic_x
  !
  !Periodic condition  
  integer, allocatable, dimension( : ), private :: periodic_y
  !
  !Periodic condition  
  integer, allocatable, dimension( : ), private :: periodic_z
  !
  !Coulomb energy of i,j in real space
  real,  allocatable, dimension(:,:,:), private :: real_ij 
  !
  !coefficients in Fourier space
  real*8,  allocatable, dimension( : ), private :: exp_ksqr
  !
  !structure factor
  complex(kind=8), allocatable, dimension( : ), private :: rho_k
  !
  !difference of structure factor
  complex(kind=8), allocatable, dimension( : ), private :: delta_rhok
  !
  !difference of structure factor
  complex(kind=8), allocatable, dimension( : ), private :: delta_cosk
  !
  !wave vector ordinal number
  integer, allocatable, dimension(:,:), private :: totk_vectk
  !
  real*8,  allocatable, dimension(:,:,:), private :: bspln_cof
                        !Coefficients of b-spline basis
  real*8,  allocatable, dimension(:,:,:), private :: dspln_cof
                        !Coefficients of derivatives of b-spline basis
  complex (kind=8), allocatable, dimension(:,:,:), private :: Q_PME
                        !Q in SPME
  complex (kind=8), allocatable, dimension(:,:,:), private :: BC_PME
                        !B*C in SPME
!########################end arrays########################!


contains


subroutine initialize_energy_parameters_Ewald
  !--------------------------------------!
  !Initial parameters are not inputted from file and compute
  !the total energy of the system.
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !Reference:
  !The computation of alpha, rc_real et al are refered to
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.304-306.
  !--------------------------------------!
  use global_variables
  implicit none
  !
  !read energy parameters from file
  call read_energy_parameters_Ewald
  !
  !
  if ( qq /= 0 ) then
    !
    !Initialize ewald parameters and array allocate.
    call Initialize_ewald_parameters
    !
    !Construct the array totk_vectk(K_total,3), and allocate
    !rho_k(K_total), delta_rhok(K_total).
!     call build_totk_vectk
!     !
!     !Construct the coefficients vector in Fourier space
!     call build_exp_ksqr
    !
    !
    call Periodic_array
    !
    !
    call pre_calculate_real_space
    !
    !build coefficient array of b-spline
    call bspln_coeffs
    !
    !build B*C in SPME
    call PME_BC

    call write_energy_parameters_Ewald

  end if

  allocate(posq(Nq,4))

end subroutine initialize_energy_parameters_Ewald


subroutine initialize_energy_arrays_ewald
  use global_variables
  implicit none
  !
  !
  if ( qq /= 0 ) then
!     !
!     !Initialize charge with lind list. From this subroutine, pos array is needed.
!     call Build_Charge_Ewald
    !
    !Initialize cell list of charge
    call Initialize_cell_list_q_Ewald
    !
    !Initialize real cell list
    call Initialize_real_cell_list_Ewald
    !
    !Construct the structure factor rho_k
!     call build_rho_k

  end if

end subroutine initialize_energy_arrays_ewald


subroutine Periodic_array
  use global_variables
  implicit none
  integer :: i

  allocate(periodic_x(-Lx2:Lx2))
  allocate(periodic_y(-Ly2:Ly2))
  allocate(periodic_z(-Lz2:Lz2))
  Periodic_x = 0
  Periodic_y = 0
  periodic_z = 0

  if (mod(Lx2,2) == 0) then
    do i = -Lx2, Lx2
      if (i<-Lx2/2) then
        periodic_x(i) = i + Lx2
      elseif (i>=Lx2/2) then
        periodic_x(i) = i - Lx2
      else
        periodic_x(i) = i
      end if
    end do
  else
    do i = -Lx2, Lx2
      if (i<-Lx2/2) then
        periodic_x(i) = i + Lx2
      elseif (i>Lx2/2) then
        periodic_x(i) = i - Lx2
      else
        periodic_x(i) = i
      end if
    end do
  end if
  do i = -Lx2, Lx2
    if (Periodic_x(i)<0) then
      periodic_x(i) = -Periodic_x(i)
    end if
  end do

  if (mod(Ly2,2) == 0) then
    do i = -Ly2, Ly2
      if (i<-Ly2/2) then
        periodic_y(i) = i + Ly2
      elseif (i>=Ly2/2) then
        periodic_y(i) = i - Ly2
      else
        Periodic_y(i) = i
      end if
    end do
  else
    do i = -Ly2, Ly2
      if (i<-Ly2/2) then
        periodic_y(i) = i + Ly2
      elseif (i>Ly2/2) then
        periodic_y(i) = i - Ly2
      else
        Periodic_y(i) = i
      end if
    end do
  end if 
  do i = -Ly2, Ly2
    if (Periodic_y(i)<0) then
      periodic_y(i) = -Periodic_y(i)
    end if
  end do

  do i = -Lz2, Lz2
    if (i>=0) then
      Periodic_z(i) = i;
    else
      Periodic_z(i) = -i;
    end if 
  end do


  open(120,file='./data/periodic_xy.txt')
    do i = -Lx2, Lx2
      write(120,*) i,periodic_x(i),Periodic_y(i)
    end do
  close(120)

  open(121,file='./data/periodic_z.txt')
    do i = -Lz2, Lz2
      write(121,*) i,periodic_z(i)
    end do
  close(121)

end subroutine Periodic_array


subroutine total_energy_ewald(EE)
  !--------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1. 
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(out) :: EE
  integer :: i, j, k, l, m, n, x1, y1, z1, x, y, z, t
  integer :: icelx, icely, icelz, ncel
  real*8 :: EE1, EE2,rr(4),q_total,Ec

  EE = 0
  Ec = 0
  !
  !real space
  do i = 1, Nq
    m = charge(i)
    if (pos(m,4)==0) cycle
    EE1 = 0
    icelx = int((pos(m,1)-1)/clx)+1
    icely = int((pos(m,2)-1)/cly)+1
    icelz = int((pos(m,3)-1)/clz)+1
    ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
    do j = 1, cell_near_list(ncel,28,1)
      icelx = cell_near_list(ncel,j,1)
      icely = cell_near_list(ncel,j,2)
      icelz = cell_near_list(ncel,j,3)
      k = hoc_r(icelx,icely,icelz)
      do while(k/=0)
        l = charge(k)
        if (l/=m) then
          x = pos(m,1)-pos(l,1)
          y = pos(m,2)-pos(l,2)
          z = pos(m,3)-pos(l,3)
          x = periodic_x(x)
          y = periodic_y(y)
          z = Periodic_z(z)
          if ((x*x+y*y+z*z)<rcc2) then
            EE1=EE1+pos(l,4)*real_ij(x,y,z)
          end if
        end if
        k = cell_list_r(k)
      end do
    end do
    EE = EE + EE1*pos(m,4)/2.D0 
    !
    !external field
    EE = EE - EF*pos(m,4)*pos(m,3)/sigma_unit
    !
    !self corrected term
    Ec = Ec - sqrt(alpha2/pi) * pos(m,4) * pos(m,4)
  end do
  !
  !fourier space
  EE = EE + Ec/Beta*lb + sum( exp_ksqr * real( conjg(rho_k) * rho_k ) )/2.D0 
  !
  !modified term of slab geometry
  EE = EE + Mz_coef * Mz**2

!   EE = sum( exp_ksqr * real( conjg(rho_k) * rho_k ) )/2.D0
  
end subroutine total_energy_ewald


subroutine Ewald_long(energy_long)
  use global_variables
  implicit none
  real*8 :: EE0,tol1,st,fn,tau_rf1
  real*8, intent(out) :: energy_long

!   tol1 = tol
!   tau_rf1 = tau_rf
!   tol = 2.1
!   tau_rf = 50                
  !
  !Initialize ewald parameters and array allocate.
  call Initialize_ewald_parameters
!   Kmax1=ceiling(Kmax1/3.)
!   Kmax2=ceiling(Kmax2/3.)
!   Kmax3=ceiling(Kmax3/3.)
  !
  !Construct the array totk_vectk(K_total,3), and allocate
  !rho_k(K_total), delta_rhok(K_total).
  call build_totk_vectk
  !
  !Construct the coefficients vector in Fourier space
  call build_exp_ksqr
  !
  !
  call pre_calculate_real_space
  !
  !Initialize cell list of charge
  call Initialize_cell_list_q_Ewald
  !
  !Initialize real cell list
  call Initialize_real_cell_list_Ewald
  !
  !Construct the structure factor rho_k
  call build_rho_k
  !
  !
  energy_long = sum( exp_ksqr * real( conjg(rho_k) * rho_k ) )/2.D0 

  deallocate(rho_k)
  deallocate(totk_vectk)

end subroutine Ewald_long


subroutine total_energy_spme(EE)
  use global_variables
  implicit none
  real*8, intent(out) :: EE
  real*8 :: EE1, EE2

  EE1 = 0
  EE2 = 0

  if (qq/=0) then
    call real_energy_spme(EE1)

    call compute_Nq_net

    if (Nq_net/=0) then

      call SPME_Ewald(EE2)

    end if
  end if

  EE = EE1 + EE2

end subroutine total_energy_spme


subroutine real_energy_spme(EE)
  use global_variables
  implicit none
  real*8, intent(out) :: EE
  integer :: i, j, k, l, m, n, x1, y1, z1, x, y, z, t
  integer :: icelx, icely, icelz, ncel
  real*8 :: EE1, EE2,rr(4),q_total,Ec,st,fn
  EE = 0
  Ec = 0
  !
  !real space
  do i = 1, Nq
    m = charge(i)
    if (pos(m,4)==0) cycle
    EE1 = 0
    icelx = int((pos(m,1)-1)/clx)+1
    icely = int((pos(m,2)-1)/cly)+1
    icelz = int((pos(m,3)-1)/clz)+1
    ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
    do j = 1, cell_near_list(ncel,28,1)
      icelx = cell_near_list(ncel,j,1)
      icely = cell_near_list(ncel,j,2)
      icelz = cell_near_list(ncel,j,3)
      k = hoc_r(icelx,icely,icelz)
      do while(k/=0)
        l = charge(k)
        if (l/=m) then
          x = pos(m,1)-pos(l,1)
          y = pos(m,2)-pos(l,2)
          z = pos(m,3)-pos(l,3)
          x = periodic_x(x)
          y = periodic_y(y)
          z = Periodic_z(z)
          if ((x*x+y*y+z*z)<rcc2) then
            EE1=EE1+pos(l,4)*real_ij(x,y,z)
          end if
        end if
        k = cell_list_r(k)
      end do
    end do
    EE = EE + EE1*pos(m,4)/2.D0 
    !
    !external field
    EE = EE - EF*pos(m,4)*pos(m,3)/sigma_unit
    !
    !self corrected term
    Ec = Ec - sqrt(alpha2/pi) * pos(m,4) * pos(m,4)
  end do

  EE = EE + Ec/Beta*lb
  !
  !modified term of slab geometry
  EE = EE + Mz_coef * Mz**2

end subroutine real_energy_spme


subroutine Delta_Energy_Ewald(DeltaE)
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
	real*8,  intent(out) :: DeltaE

  DeltaE = 0
  !
  !Compute coulomb energy change in real space
  call Delta_real_Energy(DeltaE)   

end subroutine Delta_Energy_Ewald


subroutine Delta_Energy_Ewald_add(DeltaE)
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
  real*8,  intent(out) :: DeltaE

  DeltaE = 0
  !
  !Compute Coulomb energy
  !
  !Compute coulomb energy change in real space
  call Delta_real_Energy_add(DeltaE)

end subroutine Delta_Energy_Ewald_add


subroutine Delta_Energy_Ewald_delete(DeltaE)
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
  real*8,  intent(out) :: DeltaE

  DeltaE = 0
  !
  !Compute coulomb energy change in real space
  call Delta_real_Energy_delete(DeltaE)

end subroutine Delta_Energy_Ewald_delete


subroutine Delta_real_Energy(DeltaE)
  !--------------------------------------!
  !
  !   
  !Input
  !   DeltaE
  !Output
  !   DeltaE
  !External Variables
  !   real_pair_list, real_point
  !   pos_ip0, pos_ip1, 
  !   Lx, Ly, alpha, lb, Beta
  !Reference:
  !1.The correction of Coulomb energy in slab geometry is
  !  referred to:
  !  In-Chul Yeh, Max L. Berkowitz, 'Ewald summation for 
  !  systems with slab geometry', Journal of Chemical 
  !  Physics, vol 111, number 7, (1999).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8  :: rij(3)
  real*8  :: EE1, EE2
  integer :: i,j,k,x,y,z,x1,y1,z1,t,icelx,icely,icelz,ncel

  EE1 = 0
  icelx = int((pos_ip0(1)-1)/clx)+1
  icely = int((pos_ip0(2)-1)/cly)+1
  icelz = int((pos_ip0(3)-1)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        x = pos_ip0(1) - pos(k,1)
        y = pos_ip0(2) - pos(k,2)
        z = pos_ip0(3) - pos(k,3)
        x = periodic_x(x)
        y = periodic_y(y)
        z = Periodic_z(z)
        if ((x*x+y*y+z*z)<rcc2) then
          EE1=EE1+pos(k,4)*real_ij(x,y,z)
        end if
      end if
      j = cell_list_r(j)
    end do
  end do
  EE1 = EE1 * pos_ip1(4)

  EE2 = 0
  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        x = pos_ip1(1) - pos(k,1)
        y = pos_ip1(2) - pos(k,2)
        z = pos_ip1(3) - pos(k,3)
        x = periodic_x(x)
        y = periodic_y(y)
        z = Periodic_z(z)
        if ((x*x+y*y+z*z)<rcc2) then
          EE2=EE2+pos(k,4)*real_ij(x,y,z)
        end if
      end if
      j = cell_list_r(j)
    end do
  end do
  EE2 = EE2 * pos_ip1(4)
  DeltaE = DeltaE + EE2 - EE1
  !
  !Change of correction energy in slab geometry
  dMz = pos_ip0(4) * (pos_ip1(3) - pos_ip0(3))/sigma_unit
  DeltaE = DeltaE + Mz_coef * (2*Mz*dMz + dMz*dMz)
  !
  ! External field energy
  DeltaE=DeltaE-EF*pos_ip0(4)*(pos_ip1(3)-pos_ip0(3))/sigma_unit   !simga unit

end subroutine Delta_real_Energy


subroutine Delta_real_Energy_add(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,x1,y1,z1,t,qq1,qq2,icelx,icely,icelz,ncel
  real*8 :: EE1,EE2
  real*8, dimension(3) :: rij

  DeltaE = 0
  qq1 = pos_ip1(4)
  qq2 = pos_ip1i(4)
  EE1 = 0
  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      x = pos_ip1(1) - pos(k,1)
      y = pos_ip1(2) - pos(k,2)
      z = pos_ip1(3) - pos(k,3)
      x = periodic_x(x)
      y = periodic_y(y)
      z = Periodic_z(z)
      if ((x*x+y*y+z*z)<rcc2) then
        EE1=EE1+pos(k,4)*real_ij(x,y,z)
      end if
      j = cell_list_r(j)
    end do
  end do
  EE1 = EE1 * qq1

  EE2 = 0
  icelx = int((pos_ip1i(1)-1)/clx)+1
  icely = int((pos_ip1i(2)-1)/cly)+1
  icelz = int((pos_ip1i(3)-1)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      x = pos_ip1i(1) - pos(k,1)
      y = pos_ip1i(2) - pos(k,2)
      z = pos_ip1i(3) - pos(k,3)
      x = periodic_x(x)
      y = periodic_y(y)
      z = Periodic_z(z)
      if ((x*x+y*y+z*z)<rcc2) then
        EE2=EE2+pos(k,4)*real_ij(x,y,z)
      end if
      j = cell_list_r(j)
    end do
  end do
  EE2 = EE2 * qq2
  DeltaE = DeltaE + EE1 + EE2

  !
  !interaction of the added two particles
  x = pos_ip1i(1) - pos_ip1(1)
  y = pos_ip1i(2) - pos_ip1(2)
  z = pos_ip1i(3) - pos_ip1(3)
  x = Periodic_x(x)
  y = Periodic_y(y)
  z = Periodic_z(z)
  if (x*x+y*y+z*z<rcc2) then
    DeltaE = DeltaE + qq1*qq2*real_ij(x,y,z)
  end if
  !
  !modified term of slab geometry
  dMz = pos_ip1i(4)*pos_ip1i(3)/sigma_unit + pos_ip1(4)*pos_ip1(3)/sigma_unit
  DeltaE = DeltaE + Mz_coef * (2*Mz*dMz + dMz*dMz)
  !
  !External field energy
  DeltaE=DeltaE-EF*pos_ip1i(4)*pos_ip1i(3)/sigma_unit & 
         -EF*pos_ip1(4)*pos_ip1(3)/sigma_unit
  !
  !
  DeltaE = DeltaE - sqrt(alpha2/pi)*(qq1**2+qq2**2)*lb/beta

end subroutine Delta_real_Energy_add


subroutine Delta_real_Energy_delete(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,x1,y1,z1,t,qq1,qq2,icelx,icely,icelz,ncel
  real*8 :: EE1,EE2
  real*8, dimension(3) :: rij

  DeltaE = 0
  qq1 = pos_ip0(4)         !polymer
  qq2 = pos_ip0i(4)        !ions
  !
  ! Real Space
  EE1 = 0
  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        x = pos_ip1(1) - pos(k,1)
        y = pos_ip1(2) - pos(k,2)
        z = pos_ip1(3) - pos(k,3)
        x = periodic_x(x)
        y = periodic_y(y)
        z = Periodic_z(z)
        if ((x*x+y*y+z*z)<rcc2) then
          EE1=EE1+pos(k,4)*real_ij(x,y,z)
        end if
      end if
      j = cell_list_r(j)
    end do
  end do
  EE1 = -EE1 * qq1

  EE2 = 0
  icelx = int((pos_ip1i(1)-1)/clx)+1
  icely = int((pos_ip1i(2)-1)/cly)+1
  icelz = int((pos_ip1i(3)-1)/clz)+1 
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip1) then
        x = pos_ip1i(1) - pos(k,1)
        y = pos_ip1i(2) - pos(k,2)
        z = pos_ip1i(3) - pos(k,3)
        x = periodic_x(x)
        y = periodic_y(y)
        z = Periodic_z(z)
        if ((x*x+y*y+z*z)<rcc2) then
          EE2=EE2+pos(k,4)*real_ij(x,y,z)
        end if
      end if
      j = cell_list_r(j)
    end do
  end do
  EE2 = -EE2 * qq2
  DeltaE = DeltaE + EE1 + EE2
  !
  !interaction of the added two particles
  x = pos_ip1i(1) - pos_ip1(1)
  y = pos_ip1i(2) - pos_ip1(2)
  z = pos_ip1i(3) - pos_ip1(3)
  x = Periodic_x(x)
  y = Periodic_y(y)
  z = Periodic_z(z)
  if ((x*x+y*y+z*z)<rcc2) then
    DeltaE = DeltaE + qq1*qq2*real_ij(x,y,z)
  end if 
  !
  !modified term of slab geometry
  dMz = - pos_ip0i(4)*pos_ip0i(3)/sigma_unit - pos_ip0(4)*pos_ip0(3)/sigma_unit
  DeltaE = DeltaE + Mz_coef * (2*Mz*dMz + dMz*dMz)
  !
  !External field energy
  DeltaE=DeltaE+EF*pos_ip0i(4)*pos_ip0i(3)/sigma_unit &
         +EF*pos_ip0(4)*pos_ip0(3)/sigma_unit

  DeltaE = DeltaE + sqrt(alpha2/pi)*(qq1*qq1+qq2*qq2)*lb/beta

end subroutine Delta_real_Energy_delete


subroutine update_real_cell_list_Ewald
  use global_variables
  implicit none
  integer :: icelx,icely,icelz

  icelx = int((pos_ip0(1)-1)/clx)+1
  icely = int((pos_ip0(2)-1)/cly)+1
  icelz = int((pos_ip0(3)-1)/clz)+1
  call update_real_cell_list_delete_Ewald(ip,icelx,icely,icelz)

  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  call update_real_cell_list_add_Ewald(ip,icelx,icely,icelz)

  Mz = Mz + dMz

end subroutine update_real_cell_list_Ewald


subroutine update_cell_list_pH_add_Ewald
  use global_variables
  implicit none
  integer :: icelx,icely,icelz

  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  call update_real_cell_list_add_Ewald(ip,icelx,icely,icelz)
  call update_charge_cell_list_add_Ewald(ip)

  icelx = int((pos_ip1i(1)-1)/clx)+1
  icely = int((pos_ip1i(2)-1)/cly)+1
  icelz = int((pos_ip1i(3)-1)/clz)+1
  call update_real_cell_list_add_Ewald(ip1,icelx,icely,icelz)
  call update_charge_cell_list_add_Ewald(ip1)

  Nq_net = Nq_net + 2
  Nq_net_pe = Nq_net_pe + 1

  Mz = Mz + dMz

end subroutine update_cell_list_pH_add_Ewald


subroutine update_cell_list_pH_delete_Ewald
  use global_variables
  implicit none
  integer :: icelx,icely,icelz

  icelx = int((pos_ip0(1)-1)/clx)+1
  icely = int((pos_ip0(2)-1)/cly)+1
  icelz = int((pos_ip0(3)-1)/clz)+1
  call update_real_cell_list_delete_Ewald(ip,icelx,icely,icelz)
  call update_charge_cell_list_delete_Ewald(ip)

  icelx = int((pos_ip0i(1)-1)/clx)+1
  icely = int((pos_ip0i(2)-1)/cly)+1
  icelz = int((pos_ip0i(3)-1)/clz)+1 
  call update_real_cell_list_delete_Ewald(ip1,icelx,icely,icelz)
  call update_charge_cell_list_delete_Ewald(ip1)

  Nq_net = Nq_net - 2
  Nq_net_pe = Nq_net_pe - 1

  Mz = Mz + dMz

end subroutine update_cell_list_pH_delete_Ewald


subroutine update_real_cell_list_add_Ewald(iq,icelx,icely,icelz)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer, intent(in) :: icelx, icely, icelz
  integer :: nti,bfi,ii,j       ! next particle of ii, before particle of ii
  integer :: ed, st 

  ii = inv_charge(iq)   !ii belongs to [1,Nq]

  inv_cell_list_r(ii) = 0
  if ( inv_hoc_r(icelx,icely,icelz) /=0 ) then
    inv_cell_list_r( hoc_r(icelx,icely,icelz) ) = ii
  else
    inv_hoc_r(icelx,icely,icelz) = ii
  end if

  cell_list_r(ii) = hoc_r(icelx,icely,icelz)
  hoc_r(icelx,icely,icelz) = ii

end subroutine update_real_cell_list_add_Ewald


subroutine update_real_cell_list_delete_Ewald(iq,icelx,icely,icelz)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer, intent(in) :: icelx, icely, icelz
  integer :: nti,bfi,ii,j       ! next particle of ii, before particle of ii
  integer :: ed, st 

  ii = inv_charge(iq)   !ii belongs to [1,Nq]

  bfi = cell_list_r(ii)
  nti = inv_cell_list_r(ii)

  if ( bfi/=0 .and. nti/=0 ) then        !middle
    cell_list_r(nti) = bfi
    inv_cell_list_r(bfi) = nti
  elseif ( bfi==0 .and. nti/=0 ) then    !the first one
    cell_list_r(nti) = bfi
    inv_hoc_r(icelx,icely,icelz) = nti
  elseif ( bfi/=0 .and. nti==0 ) then    !the last one
    hoc_r(icelx,icely,icelz) = bfi
    inv_cell_list_r(bfi) = nti
  else                                   !only one
    hoc_r(icelx,icely,icelz) = nti
    inv_hoc_r(icelx,icely,icelz) = bfi
  end if

end subroutine update_real_cell_list_delete_Ewald


subroutine update_charge_cell_list_add_Ewald(iq)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer :: ii, j     

  ii = inv_charge(iq)         ! ii belongs to [1,Nq]

  inv_cell_list_q(ii) = 0
  if ( cell_list_q(Nq+1)/=0 ) then
    inv_cell_list_q(cell_list_q(Nq+1)) = ii
  else
    inv_cell_list_q(Nq+1) = ii
  end if

  cell_list_q(ii) = cell_list_q(Nq+1)
  cell_list_q(Nq+1) = ii

end subroutine update_charge_cell_list_add_Ewald


subroutine update_charge_cell_list_delete_Ewald(iq)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer :: nti,bfi,ii,j       ! next particle of ii, before particle of ii

  ii = inv_charge(iq)         !ii belongs to [1,Nq]
  
  bfi = cell_list_q(ii)
  nti = inv_cell_list_q(ii)

  if ( bfi/=0 .and. nti/=0 ) then        !middle
    cell_list_q(nti) = bfi
    inv_cell_list_q(bfi) = nti
  elseif ( bfi==0 .and. nti/=0 ) then    !the first one
    cell_list_q(nti) = bfi
    inv_cell_list_q(Nq+1) = nti
  elseif ( bfi/=0 .and. nti==0 ) then
    cell_list_q(Nq+1) = bfi
    inv_cell_list_q(bfi) = nti
  else
    cell_list_q(Nq+1) = nti
    inv_cell_list_q(Nq+1) = bfi
  end if

end subroutine update_charge_cell_list_delete_Ewald


subroutine error_analysis_ewald(EE1)
  !--------------------------------------!
  !
  !Reference:
  !1.The error is estimated by Kolafa1992,
  !  JIRI KOLAFA, JOHN W. PERRAM, 'CUTOFF ERRORS IN THE 
  !  EWALD SUMMATION FORMULAE FOR POINT CHARGE SYSTEMS',
  !  Molecular Simulation, Vol. 9(5), pp.351-368, (1992).
  !2.The true error e.g. RMS energy is defined at:
  !  Ulrich Essmann, Lalith Perera et al, 'A smooth particle
  !  mesh Ewald method', Journal of Chemical Physics, Vol. 103(9),
  !  (1995).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: EE0,tol1,st,fn,tau_rf1
  real*8, intent(out) :: EE1

  tol1 = tol
  tau_rf1 = tau_rf
  tol = 4
  tau_rf = 50                
  !
  !Initialize ewald parameters and array allocate.
  call Initialize_ewald_parameters
  Kmax1=ceiling(Kmax1/3.)
  Kmax2=ceiling(Kmax2/3.)
  Kmax3=ceiling(Kmax3/3.)
  !
  !Construct the array totk_vectk(K_total,3), and allocate
  !rho_k(K_total), delta_rhok(K_total).
  call build_totk_vectk
  !
  !Construct the coefficients vector in Fourier space
  call build_exp_ksqr
  !
  !
  call pre_calculate_real_space
  !
  !Initialize cell list of charge
  call Initialize_cell_list_q_Ewald
  !
  !Initialize real cell list
  call Initialize_real_cell_list_Ewald
  !
  !Construct the structure factor rho_k
  call build_rho_k
  !
  !
  call total_energy_ewald(EE0)
  deallocate(rho_k)
  deallocate(totk_vectk)
!   call write_energy_parameters_Ewald

  tol = tol1
  tau_rf = tau_rf1
  !
  !Initialize ewald parameters and array allocate.
  call Initialize_ewald_parameters
  !
  !
  call pre_calculate_real_space
  !
  !Initialize cell list of charge
  call Initialize_cell_list_q_Ewald
  !
  !Initialize real cell list
  call Initialize_real_cell_list_Ewald
  !
  !
  call cpu_time(st)
  call total_energy_spme(EE1)
  call cpu_time(fn)
!   call write_energy_parameters_Ewald

  rmse = abs(EE1-EE0)/abs(EE0)

end subroutine error_analysis_ewald


subroutine read_energy_parameters_Ewald
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  open(unit=100, file='energy_data.txt')
    read(100,*) lb
    read(100,*) Ef
    read(100,*) ordr(1)
    read(100,*) ordr(2)
    read(100,*) ordr(3)
    read(100,*) tol
    read(100,*) tau_rf
  close(100)

end subroutine read_energy_parameters_Ewald


subroutine write_energy_parameters_Ewald
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  write(*,*) '******************  Potential  *********************'
  write(*,*) 'Kmax1      :', Kmax1
  write(*,*) 'Kmax2      :', Kmax2
  write(*,*) 'Kmax3      :', Kmax3
  write(*,*) 'K_total    :', K_total
  write(*,*) 'alpha      :', alpha
  write(*,*) 'tol        :', tol
  write(*,*) 'tau_rf     :', tau_rf
  write(*,*) 'rcc        :', rcc
  write(*,*) 'nclx       :', nclx
  write(*,*) 'ncly       :', ncly
  write(*,*) 'nclz       :', nclz
  write(*,*) 'clx        :', clx
  write(*,*) 'cly        :', cly
  write(*,*) 'clz        :', clz
  write(*,*) '****************************************************'

end subroutine write_energy_parameters_Ewald


subroutine Initialize_ewald_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rho, v_verlet

  !
  ! alpha and rcc
  alpha    = ( tau_rf * pi**3 * Nq / (Lx*Ly*Lz*Z_empty)**2 ) ** (1.D0/6)
  alpha2   = alpha * alpha
  rcc  = tol / alpha * 2
  rcc2 = rcc * rcc
  if (rcc/2<min(Lx/3,Lz/3)) then
    !
    !use verlet list in real space
    Kmax1 = ceiling(tol*Lx*alpha/pi)
    Kmax2 = ceiling(tol*Ly*alpha/pi)
    Kmax3 = ceiling(tol*Lz*Z_empty*alpha/pi)
  else
    rcc = min(Lx/3,Lz/3)-0.1
    Kmax1    = ceiling(tol*tol/pi*Lx/rcc)
    Kmax2    = ceiling(tol*tol/pi*Ly/rcc)
    Kmax3    = ceiling(tol*tol/pi*Lz*Z_empty/rcc)
    alpha    = tol / rcc
    alpha2   = alpha * alpha
    rcc = rcc*2
    rcc2 = rcc * rcc
  end if
  !
  !Cell list parameters
  nclx = int(Lx2/(rcc+0.1))     !cell numbers in x direction
  ncly = int(Ly2/(rcc+0.1))
  nclz = int(Lz2/(rcc+0.1))
  clx = 1.D0*Lx2/nclx         !cell length    
  cly = 1.D0*Ly2/ncly
  clz = 1.D0*Lz2/nclz
  Mz_coef = 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta

  Kmax1 = ceiling(Kmax1*3.)
  Kmax2 = ceiling(Kmax2*3.)
  Kmax3 = ceiling(Kmax3*3.)
  ng    = (/Kmax1,Kmax2,Kmax3/)
  gdim  = (/Lx,Ly,Lz*Z_empty/)

end subroutine Initialize_ewald_parameters


subroutine pre_calculate_real_space
  use global_variables
  implicit none
  integer :: i, j, k 
  integer :: nnl, nn_half
  real*8 :: rr, x, y, z

  nnl = nint(rcc+1)    ! cut off length, lattice unit
  if (allocated(real_ij)) deallocate(real_ij)
  allocate( real_ij(0:nnl, 0:nnl, 0:nnl) )

  do i = 0, nnl
    do j = 0, nnl
      do k = 0, nnl
        x = i/sigma_unit                      !sigma unit
        y = j/sigma_unit
        z = k/sigma_unit
        rr = sqrt(x*x+y*y+z*z)
        real_ij(i,j,k) = erfc(alpha*rr)/rr
      end do
    end do
  end do

  real_ij = real_ij * lb / beta
  real_ij(0,0,0) = 0

end subroutine pre_calculate_real_space


subroutine compute_Nq_net
  use global_variables
  implicit none
  integer :: i

  Nq_net = 0
  Nq_net_pe = 0

  do i = 1, NN
    if(pos(i,4)/=0) then
      Nq_net = Nq_net + 1
      if(i<=Npe) then
        Nq_net_pe = Nq_net_pe + 1
      end if 
    end if
  end do

end subroutine compute_Nq_net


subroutine Build_Charge_Ewald
  !--------------------------------------!
  !Initialize and charge.
  !   
  !Input
  !   pos
  !Output
  !   charge, inv_charge, Mz
  !External Variables
  !   pos, charge, NN, Nq
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j

  if (allocated(charge)) deallocate(charge)
  allocate(charge(Nq))
  if (allocated(inv_charge)) deallocate(inv_charge)
  allocate(inv_charge(NN))

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
      j = j + 1
      charge(j) = i
    end if
  end do

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
      j = j + 1
      inv_charge(i) = j
    else
      inv_charge(i) = 0
    end if
  end do

end subroutine Build_Charge_Ewald


subroutine Initialize_cell_list_q_Ewald
  use global_variables
  implicit none
  integer :: i, j, k

  if (allocated(cell_list_q)) deallocate(cell_list_q)
  if (allocated(inv_cell_list_q)) deallocate(inv_cell_list_q)
  allocate(cell_list_q(Nq+1))     ! the last one is head of the list
  allocate(inv_cell_list_q(Nq+1)) ! the last one is the head of the list

  !assume initial state, all particles are charged.
  cell_list_q = 0
  do i = 1, Nq
    if (pos(charge(i),4)/=0) then
      cell_list_q(i) = cell_list_q(Nq+1)
      cell_list_q(Nq+1) = i
    end if
  end do

  inv_cell_list_q(Nq+1) = 0
  do i = Nq, 1, -1
    if (pos(charge(i),4)/=0) then
      inv_cell_list_q(i) = inv_cell_list_q(Nq+1)
      inv_cell_list_q(Nq+1) = i
    end if
  end do

  open(112,file='./data/cell_list_q.txt')
    do i = 1, Nq + 1
      write(112,*) cell_list_q(i), inv_cell_list_q(i)
    end do
  close(112)

end subroutine Initialize_cell_list_q_Ewald



subroutine Initialize_real_cell_list_Ewald
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n, p, q, r, x, y, z
  integer :: icelx, icely, icelz

  !
  ! maxium situation, (125,125,100), 6.2Mb RAM is needed.
  if (allocated(hoc_r)) deallocate(hoc_r)
  if (allocated(inv_hoc_r)) deallocate(inv_hoc_r)
  allocate(hoc_r(nclx,ncly,nclz))
  allocate(inv_hoc_r(nclx,ncly,nclz))
  hoc_r = 0
  inv_hoc_r = 0

  if (allocated(cell_list_r)) deallocate(cell_list_r)
  if (allocated(inv_cell_list_r)) deallocate(inv_cell_list_r)
  allocate(cell_list_r(Nq))
  allocate(inv_cell_list_r(Nq))
  cell_list_r = 0
  inv_cell_list_r = 0

  Mz = 0
  do i = 1, NN
    Mz = Mz + pos(i,4)*pos(i,3)/sigma_unit !sigma unit
  end do

  do i = 1, Nq
    j = charge(i)
    if (pos(j,4)/=0) then
      icelx = int((pos(j,1)-1)/clx) + 1
      icely = int((pos(j,2)-1)/cly) + 1
      icelz = int((pos(j,3)-1)/clz) + 1
      cell_list_r(i) = hoc_r(icelx,icely,icelz)
      hoc_r(icelx,icely,icelz) = i
    end if
  end do

  do i = Nq, 1, -1
    j = charge(i)
    if (pos(j,4)/=0) then
      icelx = int((pos(j,1)-1)/clx) + 1
      icely = int((pos(j,2)-1)/cly) + 1
      icelz = int((pos(j,3)-1)/clz) + 1
      inv_cell_list_r(i) = inv_hoc_r(icelx,icely,icelz)
      inv_hoc_r(icelx,icely,icelz) = i
    end if
  end do

  !
  ! maxium situation, (125*125*100,28,3), 500Mb RAM is needed.
  if(allocated(cell_near_list)) deallocate(cell_near_list)
  allocate(cell_near_list(nclx*ncly*nclz,28,3))
  cell_near_list = 0
  m = 0
  do i = 1, nclx
    do j = 1, ncly
      do k = 1, nclz
        m = m + 1
        n = 0
        do p = -1, 1
          do q = -1, 1
            do r = -1, 1
              x = i + p
              y = j + q
              z = k + r
              if (z>0 .and. z<=nclz) then
                n = n + 1
                if (x>nclx) then
                  x = x - nclx
                elseif (x<=0) then
                  x = x + nclx
                end if
                if (y>ncly) then
                  y = y - ncly
                elseif (y<=0) then
                  y = y + ncly
                end if
                cell_near_list(m,n,1) = x
                cell_near_list(m,n,2) = y
                cell_near_list(m,n,3) = z
              end if
            end do
          end do
        end do
        cell_near_list(m,28,1) = n
      end do
    end do
  end do

  open(113,file='./data/cell_list_r.txt')
    do i = 1, Nq
      write(113,*) i, cell_list_r(i), inv_cell_list_r(i)
    end do
  close(113)


!   open(100,file='./data/hoc_r.txt')
!     do i = 1, nclx
!       do j = 1, ncly
!         do k = 1, nclz
!          write(100,*) i,j,k,hoc_r(i,j,k),inv_hoc_r(i,j,k)
!         end do
!       end do
!     end do
!   close(100)

end subroutine Initialize_real_cell_list_Ewald


subroutine build_totk_vectk
  !--------------------------------------!
  !exp_ksqr, rho_k, delta_rhok are all vectors with size of
  !K_total. For i = 1 to K_total, we often need to know 
  !corresponding wave number kx,ky,kz. This progam build a 
  !array totk_vectk(1:K_total,3) to store kx,ky,kz.
  !What's more, rho_k and delta_rhok are allocated here.
  !   
  !Input
  !   
  !Output
  !   totk_vectk
  !   K_total
  !External Variables
  !   K_total, Kmax1, Kmax2, Kmax3
  !   totk_vectk
  !Routine Referenced:
  !   
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l
  real*8  :: ksqr, k1, k2, k3, factor, kcut

  K_total=0
  do k = 0, Kmax3
    do i = -Kmax1, Kmax1
      do j = -Kmax2, Kmax2
        kcut = (1.D0*i/Kmax1) * (1.D0*i/Kmax1) &
             + (1.D0*j/Kmax2) * (1.D0*j/Kmax2) &
             + (1.D0*k/Kmax3) * (1.D0*k/Kmax3)
        if ( kcut>1 .or. kcut==0 ) cycle
        K_total = K_total + 1
      end do
    end do
  end do

  if ( allocated(totk_vectk) ) deallocate(totk_vectk)
  if ( allocated(rho_k)      ) deallocate(rho_k)
  allocate( totk_vectk( K_total, 3 ) )
  allocate( rho_k( K_total )         )
  totk_vectk = 0
  rho_k      = 0

  l=0
  do k = 0, Kmax3
    do i = -Kmax1, Kmax1
      do j = -Kmax2, Kmax2
        kcut = (1.D0*i/Kmax1) * (1.D0*i/Kmax1) &
             + (1.D0*j/Kmax2) * (1.D0*j/Kmax2) &
             + (1.D0*k/Kmax3) * (1.D0*k/Kmax3)
        if ( kcut>1 .or. kcut==0 ) cycle
        l = l + 1
        totk_vectk( l, 1 ) = i
        totk_vectk( l, 2 ) = j
        totk_vectk( l, 3 ) = k
      end do
    end do
  end do

end subroutine build_totk_vectk


subroutine build_exp_ksqr
  !--------------------------------------!
  !Reciprocal energy is divided to three parts: 
  !1.structrure factor is referred to rho_k.
  !2.difference of structure factor between new and old
  !position is referred to delta_rhok.
  !3.the other which includes exp(k^2/4/alpha) is referred 
  !to exp_ksqr.
  !This program is used to bulid the third part.
  !
  !Input
  !   
  !Output
  !   exp_ksqr
  !External Variables
  !   K_total
  !   Kmax1, Kmax2, Kmax3
  !   alpha2, lb
  !   Lx, Ly, Lz, Z_empty
  !   Beta
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.300(12.1.25),
  !however his alpha is alpha^2 in this program.b 
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l, ord(3)
  real*8  :: ksqr, k1, k2, k3, factor

  if ( allocated(exp_ksqr) ) deallocate(exp_ksqr)
  allocate( exp_ksqr(K_total) )
  exp_ksqr = 0

  l = 0
  do i = 1, K_total
    ord = totk_vectk(i,:)
    if ( ord(3) == 0 ) then
      factor = 1
    else
      factor = 2
    end if
    k1   = 2*pi*ord(1) / Lx
    k2   = 2*pi*ord(2) / Ly
    k3   = 2*pi*ord(3) / Lz/Z_empty 
    ksqr = k1*k1 + k2*k2 + k3*k3 
    exp_ksqr(i) = factor * 4*pi / (Lx*Ly*Lz*Z_empty) *  &
                  exp(-ksqr/4/alpha2) / ksqr * lb / Beta     
  end do

end subroutine build_exp_ksqr


subroutine build_rho_k
  !--------------------------------------!
  !Calculate the structure factor array.
  !   
  !Input
  !   
  !Output
  !   rho_k
  !External Variables
  !   pos, charge
  !   Nq, Lx, Ly, Lz, Z_empty, K_total
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none
  complex(kind=8) :: eikx(1:Nq, -Kmax1:Kmax1)
  complex(kind=8) :: eiky(1:Nq, -Kmax2:Kmax2)
  complex(kind=8) :: eikz(1:Nq, 0:Kmax3)
  integer i,j,l,m,n,p,q,r,ord(3)
  real*8 :: c1, c2, c3
  real*8 :: zq(Nq)
  rho_k = 0
  zq = 0
  eikx = 0
  eiky = 0
  eikz = 0

  m = cell_list_q(Nq+1)
  do while( cell_list_q(m)/=0 )
    zq(m) = pos(charge(m),4)
    m = cell_list_q(m)
  end do

  c1 = 2*pi/Lx
  c2 = 2*pi/Ly
  c3 = 2*pi/Lz/Z_empty
  m = cell_list_q(Nq+1)
  do while(m /= 0)
    i = charge(m)
    eikx(m,0)  = (1,0)
    eiky(m,0)  = (1,0)
    eikz(m,0)  = (1,0)

    eikx(m,1)  = cmplx( cos(c1*pos(i,1)/sigma_unit), &
                 sin(c1*pos(i,1)/sigma_unit), 8 )
    eiky(m,1)  = cmplx( cos(c2*pos(i,2)/sigma_unit), &
                 sin(c2*pos(i,2)/sigma_unit), 8 )
    eikz(m,1)  = cmplx( cos(c3*pos(i,3)/sigma_unit), &
                 sin(c3*pos(i,3)/sigma_unit), 8 )
    eikx(m,-1) = conjg(eikx(m,1))
    eiky(m,-1) = conjg(eiky(m,1))
    m = cell_list_q(m)
  end do

  do p=2, Kmax1
    m = cell_list_q(Nq+1)
    do while(m/=0)
      eikx(m,p)=eikx(m,p-1)*eikx(m,1)
      eikx(m,-p)=conjg(eikx(m,p))
      m = cell_list_q(m)
    end do
  end do
  do q=2, Kmax2
    m = cell_list_q(Nq+1)
    do while(m/=0)
      eiky(m,q)=eiky(m,q-1)*eiky(m,1)
      eiky(m,-q)=conjg(eiky(m,q))
      m = cell_list_q(m)
    end do
  end do
  do r=2, Kmax3
    m = cell_list_q(Nq+1)
    do while(m/=0)
      eikz(m,r)=eikz(m,r-1)*eikz(m,1)
      m = cell_list_q(m)
    end do
  end do

  do i = 1, K_total
    ord = totk_vectk(i,:)
    m = cell_list_q(Nq+1)
    do while(m/=0)
      rho_k(i) = rho_k(i) + &
                 zq(m) * eikx(m,ord(1)) * eiky(m,ord(2)) * eikz(m,ord(3))
      m = cell_list_q(m)
    end do
  end do

end subroutine build_rho_k

!
!The following program are the SPME program.
!

subroutine SPME_Ewald(Enrg)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
  implicit none
  real*8, intent(out) :: Enrg
  real*8, dimension(Nq_net,ordr(1)) :: SPx
  real*8, dimension(Nq_net,ordr(1)) :: dSPx
  real*8, dimension(Nq_net,ordr(2)) :: SPy
  real*8, dimension(Nq_net,ordr(2)) :: dSPy
  real*8, dimension(Nq_net,ordr(3)) :: SPz
  real*8, dimension(Nq_net,ordr(3)) :: dSPz
  integer, dimension(Nq_net,ordr(1)) :: iqmap
  integer, dimension(Nq_net,ordr(2)) :: jqmap
  integer, dimension(Nq_net,ordr(3)) :: kqmap
  real*8 :: st,fn,tm
  real*8, dimension(Nq_net,3) :: acc_f
  integer :: i,m

  if (Nq_net/=0 .and. qq/=0) then
    if (allocated(posq)) deallocate(posq)
    allocate(posq(Nq_net,4))
    m = 0
    do i = 1, NN
      if (pos(i,4)/=0) then
        m = m + 1
        posq(m,1:3) = pos(i,1:3)/sigma_unit
        posq(m,4) = pos(i,4)*1.D0
      end if
    end do

    call cpu_time(st)
    
    call splcof(SPx, dSPx, iqmap, ordr(1), 1)
    call splcof(SPy, dSPy, jqmap, ordr(2), 2)
    call splcof(SPz, dSPz, kqmap, ordr(3), 3)

    call MapCharges(SPx, SPy, SPz, iqmap, jqmap, kqmap)

    call pmeOrthoConvBC(Enrg)

    call cpu_time(fn)

    Enrg = Enrg * lb / beta
  else
    Enrg = 0
  end if

end subroutine SPME_Ewald



subroutine splcof(SP, dSP, qmap, ord, xyz)
  !-----------------------------------------!
  !ord指定数组大小可以吗，回去查查, 可以的
  !不能在定义时候初始化
  !-----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: ord
  real*8, dimension(Nq_net,ord), intent(out) :: SP
  real*8, dimension(Nq_net,ord), intent(out) :: dSP
  integer, dimension(Nq_net,ord), intent(out) :: qmap
  integer, intent(in) :: xyz
  real*8, dimension(Nq_net) :: uu,ww
  integer, dimension(Nq_net) :: uuc
  integer :: i,j,k
  SP=0
  dSP=0
  qmap=0
  
  uu = ng(xyz)/gdim(xyz)*posq(:,xyz)            !化到标定后的分数坐标
  uuc = ceiling(uu)                           !ceiling can also be used for array
  
  if (ord==4) then
    ww = uu-uuc+1
    
    SP(:,3) = 0.5*ww*ww
    SP(:,1) = 0.5*(1.0-ww)*(1.0-ww)
    SP(:,2) = 1.0-SP(:,1)-SP(:,3)

    dSP(:,1) = -SP(:,1)
    dSP(:,2) = SP(:,1) - SP(:,2)
    dSP(:,3) = SP(:,2) - SP(:,3)
    dSP(:,4) = SP(:,3)
    
    SP(:,4) = ww*SP(:,3)/3
    SP(:,3) = ((ww+1.0)*SP(:,2) + (3.0-ww)*SP(:,3))/3
    SP(:,1) = (1.0-ww)*SP(:,1)/3
    SP(:,2) = 1.0-SP(:,1)-SP(:,3)-SP(:,4)
  else
    ww = uu-uuc+1
    do k=1,ord
      SP(:,ord-k+1)=bspln_cof(k,1,xyz)*ww
      dSP(:,ord-k+1)=dspln_cof(k,1,xyz)*ww  
      do j=2,ord-2
        SP(:,ord-k+1)=ww*(SP(:,ord-k+1)+bspln_cof(k,j,xyz))
        dSP(:,ord-k+1)=ww*(dSP(:,ord-k+1)+dspln_cof(k,j,xyz))
      end do
      SP(:,ord-k+1)=ww*(SP(:,ord-k+1)+bspln_cof(k,ord-1,xyz))
      SP(:,ord-k+1)=SP(:,ord-k+1)+bspln_cof(k,ord,xyz)
      dSP(:,ord-k+1)=dSP(:,ord-k+1)+dspln_cof(k,ord-1,xyz)
      ww=ww+1   
    end do
  end if

  do i = 0,ord-1
    j = i + 1
    qmap(:,j) = uuc(:) - ord + i
    qmap(:,j) = qmap(:,j) - ng(xyz)*floor(1.*qmap(:,j)/ng(xyz))
  end do
  
  qmap = qmap + 1
  
  if (xyz == 1) then
    do i = 1,ord
      do j = 1,Nq_net
        SP(j,i) = SP(j,i)*posq(j,4)
        dSP(j,i) = dSP(j,i)*posq(j,4)
      end do
    end do
  end if
  
end subroutine splcof


subroutine MapCharges(SPx, SPy, SPz, iqmap, jqmap, kqmap)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
  implicit none
  real*8,  dimension(Nq_net,ordr(1)), intent(in) :: SPx
  real*8,  dimension(Nq_net,ordr(2)), intent(in) :: SPy
  real*8,  dimension(Nq_net,ordr(3)), intent(in) :: SPz
  integer, dimension(Nq_net,ordr(1)), intent(in) :: iqmap
  integer, dimension(Nq_net,ordr(2)), intent(in) :: jqmap
  integer, dimension(Nq_net,ordr(3)), intent(in) :: kqmap
  integer :: h,i,j,k,ii,jj,kk,t,ng1o,ng2o,ng3o,ih,jh,kh
  real*8  :: dqx,dqxy,dd,qmij(ordr(1),ordr(2))
  real*8  :: SPxh(ordr(1),1),SPyh(1,ordr(2)),SPzh(ordr(3))
  integer :: imp(ordr(1)),jmp(ordr(2)),kmp(ordr(3))
  real*8, dimension(ordr(1),ordr(2),ordr(3)) :: qmijk

  Q_PME=0
  
  ng1o=ng(1)-ordr(1)+2
  ng2o=ng(2)-ordr(2)+2
  ng3o=ng(3)-ordr(3)+2
  qmijk=0
  do h=1,Nq_net
    imp=iqmap(h,:)
    jmp=jqmap(h,:)
    kmp=kqmap(h,:)
    ih=imp(1)
    jh=jmp(1)
    kh=kmp(1) 
    SPxh(:,1)=SPx(h,:)
    SPyh(1,:)=SPy(h,:)
    SPzh=SPz(h,:)
    qmij=matmul(SPxh,SPyh)
    do i=1,ordr(3)
      qmijk(:,:,i)=SPzh(i)*qmij
    end do
    Q_PME(imp,jmp,kmp)=Q_PME(imp,jmp,kmp)+qmijk
  end do

end subroutine Mapcharges


subroutine pmeOrthoConvBC(Enrg)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
  implicit none
  include "fftw3.f90"
  real*8, intent(out) :: Enrg
  complex ( kind = 8 ), dimension(ng(1),ng(2),ng(3)) :: fQ_PME
  integer i,j,k
  integer ( kind = 8 ) plan_backward
  integer ( kind = 8 ) plan_forward
  integer ( kind = 8 ) plan
  real*8 :: st,fn
  
  call dfftw_plan_dft_3d_ ( plan_forward, ng(1), ng(2), ng(3), Q_PME, fQ_PME, FFTW_FORWARD, FFTW_Estimate ) 

  call dfftw_execute_ ( plan_forward )
  
  call dfftw_destroy_plan_ ( plan_forward )
  
  fQ_PME(1,1,1)=0   

  Enrg=sum(fQ_PME*conjg(fQ_PME)*BC_PME)/(ng(1)*ng(2)*ng(3))/2

end subroutine pmeOrthoConvBC


subroutine bspln_coeffs
  !----------------------------------------!
  !Calculate the coefficients of the b-spline basis function.
  !2010, Gradimir V. Milovanovic, Applied Mathematics Letters
  !----------------------------------------!
  use global_variables
  implicit none
  integer :: h,i,j,k,m,g,ii
  real*8, allocatable, dimension(:,:,:) :: a
  if ( allocated(bspln_cof) ) deallocate(bspln_cof)
  if ( allocated(dspln_cof) ) deallocate(dspln_cof)
  allocate(bspln_cof(maxval(ordr),maxval(ordr),3))
  allocate(dspln_cof(maxval(ordr),maxval(ordr),3))
  bspln_cof=0
  dspln_cof=0
  
  do h=1,3
    allocate(a(ordr(h),ordr(h),ordr(h)+1))
    a=0
    a(1,1,2)=1
    do m=2,ordr(h)
      a(m-1+1,m,-1+2)=0
      a(m-2+1,m,-1+2)=0
      g=floor(m/2.D0)-1
      do k=0,g
        a(m-1+1,m,k+2) = 1.D0/(m-1)/(m-k)* ( 1.D0*m*a(m-2+1,m-1,k+2) - 1.D0*k*(m-1)*a(m-1+1,m,k-1+2) &
                         - 1.D0*a(m-2+1,m,k-1+2) )
        a(m-1+1,m,m-k-1+2) = (-1.D0)**(m-1) * a(m-1+1,m,k+2)
        do i=m-2,0,-1
          a(i+1,m,k+2) = 1.D0*m/(i+1-m) * ( 1.D0*(i+1)*a(i+1+1,m,k+2)-1.D0*a(i+1,m-1,k+2) )
          do ii=0,m-i-1
            a(i+1,m,m-k-1+2) = 1.D0*a(i+1,m,m-k-1+2) + factorial(ii+i)/factorial(ii)  &
                               * a(i+ii+1,m,k+2)*m**ii
          end do
          a(i+1,m,m-k-1+2) = 1.D0*a(i+1,m,m-k-1+2) * (-1)**i/factorial(i)
        end do
      end do
      if (mod(m,2)/=0) then
        a(m-1+1,m,g+1+2) = 1.D0/(m-1)/(m-(g+1)) * ( 1.D0*m*a(m-2+1,m-1,g+1+2) - 1.D0*(g+1)*(m-1) &
                           * a(m-1+1,m,(g+1)-1+2) - 1.D0*a(m-2+1,m,(g+1)-1+2) )
        do i=m-2,0,-1
          a(i+1,m,g+1+2) = 1.D0*m/(i+1-m) * ( 1.D0*(i+1)*a(i+1+1,m,g+1+2)-1.D0*a(i+1,m-1,g+1+2) )
        end do
      end if
    end do

    do i=0,ordr(h)-1
      do k=0,ordr(h)-1
        bspln_cof(k+1,i+1,h)=a(ordr(h)-i,ordr(h),k+2) 
      end do
    end do
    do i=1,ordr(h)
      do j=1,ordr(h)
        dspln_cof(i,j,h)=bspln_cof(i,j,h)*(ordr(h)-j)
      end do
    end do
    deallocate(a)
  end do

end subroutine bspln_coeffs


subroutine PME_BC
  !----------------------------------------!
  !
  !----------------------------------------!
  use global_variables
  implicit none
  complex (kind=8) :: bx(ng(1)), by(ng(2)), bz(ng(3))
  real*8 :: mx(ng(1)), my(ng(2)), mz(ng(3))
  real*8 :: ex(ng(1)), ey(ng(2)), ez(ng(3))
  real*8 dx, dxy, dxyz
  integer :: i,j,k
  
  if ( allocated(BC_PME) ) deallocate(BC_PME)
  if ( allocated(Q_PME)  ) deallocate(Q_PME)

  allocate( BC_PME(Kmax1, Kmax2, Kmax3) )
  allocate( Q_PME(Kmax1, Kmax2, Kmax3)  )

  BC_PME = 0
  Q_PME  = 0
  
  call pmeOrthoTabBC(bx, mx, ex, 1)
  call pmeOrthoTabBC(by, my, ey, 2)
  call pmeOrthoTabBC(bz, mz, ez, 3)
  
  do i=1,ng(1)
    dx = mx(i)*mx(i)
    do j=1,ng(2)
      dxy = dx + my(j)*my(j)
      do k = 1,ng(3)
        dxyz = dxy + mz(k)*mz(k)
        if (i+j+k>3) then                         !避免除以0
          BC_PME(i,j,k) = bx(i)*conjg(bx(i)) * by(j)*conjg(by(j)) *        &
                          bz(k)*conjg(bz(k)) * ex(i)*ey(j)*ez(k) / dxyz
        end if
      end do
    end do
  end do
end subroutine PME_BC


subroutine pmeOrthoTabBC(bv, mv, ev, xyz)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: xyz
  complex (kind=8), dimension(ng(xyz)), intent(out) :: bv
  real*8, dimension(ng(xyz)), intent(out) :: mv
  real*8, dimension(ng(xyz)), intent(out) :: ev
  real*8 :: tpS,pivol
  integer :: i,j,ord
  complex (kind=8) :: bsp
  ord=ordr(xyz)
  bv=0
  mv=0
  ev=0
  
  if (xyz == 1) then
    pivol = ng(1)*ng(2)*ng(3)/(pi*gdim(1)*gdim(2)*gdim(3))
  end if
  tpS=pi/alpha
  do i = 1,ng(xyz)
    if (i-1 <= ng(xyz)/2) then            !这个m的取值范围有些讲究
      mv(i) = (i-1)/gdim(xyz)
    else
      mv(i) = (i-1-ng(xyz))/gdim(xyz)
    end if
    if (xyz == 1) then
      ev(i) = pivol*exp(-tpS*tpS*mv(i)*mv(i))
    else
      ev(i) = exp(-tpS*tpS*mv(i)*mv(i))
    end if
  end do
  
  do i=0, ng(xyz)-1
    bsp=0
    do j=0, ord-2
      bsp = bsp + bspln(j+1.D0,ord) * &
            cmplx( cos(2*pi*i*j/ng(xyz)), sin(2*pi*i*j/ng(xyz)), 8 )
    end do
    if (abs(bsp)<1e-10) then
      bsp=1e15
    end if
    bv(i+1)=cmplx( cos(2*pi*i*(ord-1)/ng(xyz)), sin(2*pi*i*(ord-1)/ng(xyz)), &
                   8 ) / bsp !complex number can be divided
  end do
  
end subroutine pmeOrthoTabBC


Function factorial(n)
  integer, intent(in) :: n
  real*8 :: factorial
  integer :: i
  
  factorial=1.D0
  if (n/=0) then
    do i=n,1,-1
      factorial=factorial*i
    end do
  end if
  
end function factorial


recursive real*8 Function bspln(x,ord) result(y)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  implicit none
  real*8, intent(in)  :: x
  integer, intent(in) :: ord

  if (ord > 2) then
    y = (x/(ord-1))*bspln(x,ord-1) + ((ord-x)/(ord-1))*bspln(x-1.0,ord-1)
  elseif (ord == 2) then
    if (x >= 0.0 .and. x <= 2.0) then
      y = 1.0 - abs(x-1.0)
    else
      y = 0.0
    end if
  else
    write(*,*) 'bspln error!  Order = ', ord
  end if

end Function bspln


end module compute_energy_ewald


















