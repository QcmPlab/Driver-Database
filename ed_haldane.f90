program ed_haldane
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS

  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso,Nlat
  logical                                       :: converged

  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: halHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk
  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8)                                       :: a,a0,bklen

  real(8),dimension(2)                          :: d1,d2,d3
  real(8),dimension(2)                          :: a1,a2,a3
  real(8),dimension(2)                          :: bk1,bk2,pointK,pointKp

  !variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: ts,tsp,phi,delta,Mh,wmixing
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym
  !
  real(8),dimension(2)                          :: Eout
  real(8),allocatable,dimension(:)              :: dens
  !
  real(8),dimension(:,:),allocatable            :: Zmats
  complex(8),dimension(:,:,:),allocatable       :: Zfoo
  complex(8),allocatable,dimension(:,:,:,:,:)   :: S0


  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputHALDANE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputHALDANE.conf",default=1d0)
  call parse_input_variable(tsp,"TSP","inputHALDANE.conf",default=1.d0/3/sqrt(3d0))
  call parse_input_variable(mh,"MH","inputHALDANE.conf",default=0d0)
  call parse_input_variable(phi,"PHI","inputHALDANE.conf",default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  phi=phi*pi

  if(Nspin/=1.OR.Norb/=1)stop "Wrong setup from input file: Nspin=Norb=1"
  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso

  !see (ed_graphene for backupd version)
  !LATTICE BASIS:
  ! nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]
  !
  !
  !next nearest-neighbor displacements: A-->A, B-->B, cell basis
  a1 = d2-d3                    !a*sqrt(3)[sqrt(3)/2,-1/2]
  a2 = d3-d1                    !a*sqrt(3)[-sqrt(3)/2,-1/2]
  a3 = d1-d2                    !a*sqrt(3)[0, 1]
  !
  !
  !RECIPROCAL LATTICE VECTORS:
  bklen=4d0*pi/sqrt(3d0)
  bk1=bklen*[ sqrt(3d0)/2d0 ,  1d0/2d0 ]
  bk2=bklen*[ sqrt(3d0)/2d0 , -1d0/2d0 ]


  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]



  !Allocate Weiss Field:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  allocate(S0(Nlat,Nspin,Nspin,Norb,Norb));S0=zero
  allocate(Zmats(Nlso,Nlso));Zmats=eye(Nlso)
  allocate(Zfoo(Nlat,Nso,Nso));Zfoo=0d0

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(halHloc,Nlat,Nspin,Norb)

  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver(Bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath,Hloc,iprint=1)
     call ed_get_sigma_matsubara(Smats,Nlat)
     call ed_get_sigma_real(Sreal,Nlat)
     S0 = dreal(Smats(:,:,:,:,:,1))
     call ChernNumber()
     call build_EigenBands()

     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=4)
     call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal,iprint=4)

     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc,iprint=4)
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc,iprint=4)
     endif

     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf(Bath,Weiss,Hloc,ispin=1)

     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo



contains




  !--------------------------------------------------------------------!
  !Haldane HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_haldane_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                       :: h0,hx,hy,hz
    real(8)                       :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !
    h0 = -2*tsp*cos(phi)*sum( cos(kdota(:)) )
    hx =-ts*sum( cos(kdotd(:)) )
    hy =-ts*sum( sin(kdotd(:)) )
    hz = -2*tsp*sin(phi)*sum( sin(kdota(:)) ) + Mh 
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
    hk = hk  +  nnn2lso_reshape(S0,Nlat,Nspin,Norb)
  end function hk_haldane_model






  !---------------------------------------------------------------------
  !PURPOSE: Get Haldane Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional             :: file
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky    
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    real(8)                               :: foo,n0(Nlso)
    integer                               :: unit
    complex(8),dimension(Nlso,Nlso,Lmats) :: Gmats,fooSmats
    complex(8),dimension(Nlso,Nlso,Lreal) :: Greal,fooSreal
    real(8),dimension(:),allocatable      :: kxgrid,kygrid
    real(8),dimension(2)                  :: kvec
    real(8),dimension(:,:),allocatable    :: kpath
    write(LOGfile,*)"Build H(k) for Haldane model (now also Nobel Laureate):"
    Lk=Nk*Nk
    write(*,*)"# of k-points     :",Lk
    write(*,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    !
    allocate(Hk(Nlso,Nlso,Lk))
    allocate(Wtk(Lk))
    print*,"Build Hk(Nlso,Nlso) for the Graphene model"
    ik=0
    do iy=1,Nk
       ky = dble(iy-1)/Nk
       do ix=1,Nk
          ik=ik+1
          kx=dble(ix-1)/Nk
          kvec = kx*bk1 + ky*bk2
          Hk(:,:,ik) = hk_haldane_model(kvec,Nlso)
       enddo
    enddo
    Wtk = 1d0/Lk
    !
    allocate(halHloc(Nlso,Nlso))
    halHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(halHloc))<1.d-9)halHloc=0d0
    call write_Hloc(halHloc)
    !
    !
    allocate(Kpath(4,2))
    KPath(1,:)=[0,0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointKp
    KPath(4,:)=[0d0,0d0]
    call TB_Solve_path(hk_haldane_model,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.nint")
    !
    !
    call ChernNumber()
    !
    !Build the local GF:
    fooSmats=zero
    call dmft_gloc_matsubara(Hk,Wtk,Gmats,fooSmats,iprint=4)
    fooSreal=zero
    call dmft_gloc_realaxis(Hk,Wtk,Greal,fooSreal,iprint=4)
    do iorb=1,Nlso
       n0(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
    enddo
    open(10,file="density.dat")
    write(10,"(10F20.12)")(n0(iorb),iorb=1,Nlso),sum(n0)
    close(10)
    !
  end subroutine build_hk




  subroutine build_EigenBands()
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky    
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    real(8)                               :: foo,n0(Nlso)
    integer                               :: unit
    complex(8),dimension(Nlso,Nlso,Lmats) :: Gmats,fooSmats
    complex(8),dimension(Nlso,Nlso,Lreal) :: Greal,fooSreal
    real(8),dimension(:),allocatable      :: kxgrid,kygrid
    real(8),dimension(2)                  :: kvec
    real(8),dimension(:,:),allocatable    :: kpath
    !
    allocate(Kpath(4,2))
    KPath(1,:)=[0,0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointKp
    KPath(4,:)=[0d0,0d0]
    call TB_Solve_path(hk_haldane_model,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.dat")
    !
    !
  end subroutine build_EigenBands



  subroutine ChernNumber
    complex(8),dimension(:,:,:),allocatable :: BlochStates
    real(8),dimension(:,:),allocatable      :: Berry_curvature
    integer                                 :: ik,ix,iy
    complex(8)                              :: Eigvec(Nlso,Nlso)
    real(8)                                 :: Chern,BZ_area,eigval(Nlso)
    real(8),dimension(:),allocatable        :: kxgrid
    !CHERN NUMBERS:
    allocate(BlochStates(Nlso,Nk,Nk))
    allocate(Berry_curvature(Nk,Nk))
    chern=0d0
    ik=0
    do ix=1,Nk
       do iy=1,Nk
          ik=ik+1
          Eigvec = Hk(:,:,ik) +  nnn2lso_reshape(S0,Nlat,Nspin,Norb)
          call eigh(Eigvec,Eigval)
          BlochStates(:,ix,iy) = Eigvec(:,1)
       enddo
    enddo
    call get_Chern_number(BlochStates,Chern,Berry_curvature,Nk/pi2*Nk/pi2)
    write(LOGfile,*)"Chern =",chern
    allocate(kxgrid(Nk))
    kxgrid=linspace(0d0,pi2,Nk)
    call splot3d("Berry_Curvature.nint",kxgrid,kxgrid,Berry_Curvature)
    open(10,file="ChernNumber.dat")
    write(10,"(I3,F16.12)")nint(chern),chern
    close(10)
  end subroutine ChernNumber


  subroutine Get_Chern_number(State,Chern_number,Berry_curvatures,one_over_area)
    complex(8),intent(in),dimension(:,:,:)                     :: state !(nhilb, nk1, nk2)
    real(8),intent(out)                                        :: chern_number
    real(8),intent(out),dimension(size(state,2),size(state,3)) :: berry_curvatures
    real(8),intent(in)                                         :: one_over_area
    integer                                                    :: nhilb,nk1,nk2
    integer                                                    :: i1,i2,i3,ix,i1p,i1m,i2m,i2p,it
    complex(8)                                                 :: path_scalar_products(4)
    real(8)                                                    :: berry_phase
    !
    Nhilb= size(state,1)
    Nk1  = size(state,2)
    Nk2  = size(state,3)
    chern_number = zero
    do i1= 1, nk1
       i1p = modulo(i1,nk1) + 1
       i1m = modulo(i1-2,nk1) + 1
       do i2= 1, nk2           !faccio l'integrale sulla bz
          i2p = modulo(i2,nk2) + 1
          i2m = modulo(i2-2,nk2) + 1
          path_scalar_products(1) = dot_product(state(:,i1,i2),state(:,i1, i2p))
          path_scalar_products(2) = dot_product(state(:,i1,i2p),state(:,i1p, i2p))
          path_scalar_products(3) = dot_product(state(:,i1p,i2p),state(:,i1p, i2))
          path_scalar_products(4) = dot_product(state(:,i1p,i2),state(:,i1,i2))
          berry_phase = -dimag(zlog( product(path_scalar_products)  ))
          berry_curvatures(i1,i2) = berry_phase*one_over_area
          chern_number = chern_number + berry_phase
       enddo
    enddo
    chern_number = chern_number/pi2
  end subroutine Get_Chern_number



end program ed_haldane



