program ed_triang_pxpy
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE SF_MPI
  implicit none
  !
  integer                                   :: Nkx,Nky,Nk,Nkpath,Npts,L,Nlso
  integer                                   :: iloop
  logical                                   :: converged
  !
  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)   :: Hk     !the Hamiltonian H(k)
  complex(8),allocatable,dimension(:,:)     :: pxpyHloc !loc Hamiltonian [Nlso][Nlso]
  complex(8),allocatable,dimension(:,:,:,:) :: Hloc   !loc Hamiltonian [Nspin][Nspin][Norb][Norb]
  real(8),allocatable,dimension(:)          :: Wtk    !weight of the k-points
  complex(8),allocatable                    :: sigmaPXPY(:,:),Zmats(:,:)
  !
  !Bath:
  integer                                   :: Nb
  real(8),allocatable                       :: Bath(:),Bath_(:)
  !
  !The local hybridization function:
  complex(8),allocatable                    :: Delta(:,:,:,:,:)
  complex(8),allocatable                    :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable                    :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !
  !parameters for the model:
  real(8)                                   :: Vsigma,Vpi,LamISB,LamSOC,wmixing
  !
  character(len=16)                         :: finput
  character(len=32)                         :: hkfile
  logical                                   :: spinsym,usez
  real(8),allocatable                       :: dens(:)
  !
  complex(8),dimension(4,4)                 :: Gamma0,GammaX,GammaY,GammaZ,GammaS

  !MPI Vars:
  integer                                   :: irank,comm,rank,size2
  logical                                   :: master

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  size2 = get_Size_MPI(comm)
  master = get_Master_MPI(comm)

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputPXPY.conf')
  !
  call parse_input_variable(Nkx,"NKX",finput,default=10,comment="# of k points along X")
  call parse_input_variable(Nky,"NKY",finput,default=10,comment="# of k points along Y")
  call parse_input_variable(nkpath,"NKPATH",finput,default=100,comment="# of points along BZ path segment")
  call parse_input_variable(L,"L",finput,default=2048,comment="# of frequencies")
  !
  call parse_input_variable(Vsigma,"Vsigma",finput,default=1d0)
  call parse_input_variable(Vpi,"Vpi",finput,default=-1d0)
  call parse_input_variable(LamISB,"LAMISB",finput,default=0.1d0)
  call parse_input_variable(LamSOC,"LAMSOC",finput,default=0d0)
  !
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  !
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.false.)
  call parse_input_variable(usez,"USEZ",finput,default=.false.)
  !
  call ed_read_input(trim(finput),MPI_COMM_WORLD)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nlso=Nspin*Norb


  Gamma0 = kron_pauli( pauli_tau_0, pauli_sigma_0)
  GammaX = kron_pauli( pauli_tau_x, pauli_sigma_0)
  GammaY = kron_pauli( pauli_tau_y, pauli_sigma_0)
  GammaZ = kron_pauli( pauli_tau_z, pauli_sigma_0)
  GammaS = kron_pauli( pauli_tau_y, pauli_sigma_z)

  !Allocate Weiss Field:
  allocate(Delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(dens(Nlso))
  allocate(SigmaPXPY(Nlso,Nlso)) !extrapolation of Sigma(iw_1) 
  allocate(Zmats(Nlso,Nlso))     !corresponding Z matrix

  !Buil the Hamiltonian on a grid or on  path
  !<TO BE UPDATED:
  call build_hk(trim(hkfile))



  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(comm,bath,Hloc=j2so(pxpyHloc)) ![Nlso,Nlso]-->[Nspin,Nspin,Norb,Norb]

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     call ed_get_dens(dens)

     call dmft_gloc_matsubara(comm,Hk,Wtk,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     if(master)then
        call dmft_self_consistency(Gmats,Smats,Delta,j2so(pxpyHloc),cg_scheme)
        call dmft_print_gf_matsubara(Delta,"Weiss",iprint=1)
     endif
     call Bcast_MPI(comm,Delta)

     !Fit the new bath, starting from the old bath + the supplied delta

     call ed_chi2_fitgf(comm,delta,bath,ispin=1)
     if(.not.spinsym)then
        call ed_chi2_fitgf(comm,delta,bath,ispin=2)
     else
        call spin_symmetrize_bath(bath,save=.true.)
     endif
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath
     !
     if(master)then
        converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
        if(nread/=0.d0)call search_chemical_potential(xmu,sum(dens),converged)
     endif
     call Bcast_MPI(comm,bath)
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)

     if(master)call end_loop
  enddo

  call dmft_gloc_realaxis(comm,Hk,Wtk,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)


  call dmft_kinetic_energy(comm,Hk,Wtk,Smats)

  ! call solve_hk_topological(so2j(Smats(:,:,:,:,1),Nlso))

  call finalize_MPI()

contains



  !---------------------------------------------------------------------
  !PURPOSE: GET BHZ HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j
    integer                             :: iorb,jorb
    integer                             :: isporb,jsporb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    !
    !STEP 1: set up the direct- and reciprocal-lattice vectors basis (a=1)
    !> setup real space basis
    call TB_set_ei([1d0,0d0],[-0.5d0,sqrt(3d0)/2d0])
    !> get reciprocal lattice basis
    call TB_build_bk(verbose=.true.)
    !
    if(master)write(LOGfile,*)"Build H(k) for PxPy triangular model:"
    Nk = Nkx*Nky
    if(master)write(*,*)"# of k-points     :",Nk
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nlso,Nlso,Nk))
    allocate(Wtk(Nk))
    call TB_build_model(Hk,hk_triang_pxpy,Nlso,[Nkx,Nky],wdos=.false.)
    Wtk = 1d0/Nk

    if(master.AND.present(file))then
       call TB_write_hk(Hk,trim(file),&
            Nlat=1,&
            Nspin=Nspin,&
            Norb=Norb,&
            Nkvec=[Nkx,Nky])
    endif
    !
    allocate(pxpyHloc(Nlso,Nlso))
    pxpyHloc = sum(Hk(:,:,:),dim=3)/Nk
    where(abs(dreal(pxpyHloc))<1.d-9)pxpyHloc=0d0
    if(master)  call TB_write_Hloc(pxpyHloc)
  end subroutine build_hk




  !---------------------------------------------------------------------
  !PURPOSE: GET THE BHZ HAMILTONIAN ALONG THE Gamma-X-M-Gamma path
  !---------------------------------------------------------------------
  subroutine build_hk_GXMG(kpath_)
    integer                            :: i,j
    integer                            :: Npts
    real(8),dimension(:,:),optional        :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    character(len=64)                      :: file
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is constructed along this path.
    if(present(kpath_))then
       if(master)write(LOGfile,*)"Build H(k) BHZ along a given path:"
       Npts = size(kpath_,1)
       Nk=(Npts-1)*Nkpath
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="Eig_path.nint"
    else
       if(master)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
       Npts = 5
       allocate(Kpath(Npts,2))
       KPath(1,:)=[0d0,0d0]
       KPath(2,:)=[0.5d0,0d0]
       Kpath(3,:)=[1d0/3d0,1d0/3d0]
       KPath(4,:)=[-1d0/3d0,2d0/3d0]
       Kpath(5,:)=[0d0,0d0]
       Kpath=Kpath*pi2
       file="PxPyBands"
    endif
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nlso,Nlso,Nk))
    allocate(wtk(Nk))
    ! call set_sigmaPxPy()
    call TB_build_model(Hk,hk_triang_pxpy,Nlso,kpath,Nkpath)
    Wtk = 1d0/Nk
    if(master)   call TB_Solve_model(hk_triang_pxpy,Nlso,KPath,Nkpath,&
         colors_name=[red,blue,orange,green],&
         points_name=[character(len=10) :: "{/Symbol G}","M","K","K`","{/Symbol G}"],&
         file=reg(file),iproject=.false.)

  end subroutine build_hk_GXMG


  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  !
  subroutine set_SigmaPxPy(sigma)
    complex(8),dimension(Nlso,Nlso),optional :: sigma(Nlso,Nlso)
    sigmaPxPy = zero;if(present(sigma))sigmaPxPy=sigma
  end subroutine set_SigmaPxPy



  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  !Fortran FUNCTION with two arguments mandatory: Kvec (dble, dim=2), N (integer) 
  function hk_triang_pxpy(kvec,N) result(Hk)
    !input
    real(8),dimension(:)      :: kvec
    integer                   :: N
    !output
    complex(8),dimension(N,N) :: Hk !the matrix H(k) for a given k=kvec
    !aux.
    real(8) :: kx,ky
    real(8) :: cx,cy,cxy
    real(8) :: sx,sy,sxy
    !
    if(N/=4) stop "hk_triang_pxpy error: N != 4"
    !
    kx  = kvec(1)
    ky  = kvec(2)
    cx  = cos(kx)
    cy  = cos(ky)
    cxy = cos(kx+ky)
    sx  = sin(kx)
    sy  = sin(ky)
    sxy = sin(kx+ky)
    !
    Hk = (Vsigma+Vpi)*(cx + cy + cxy)*Gamma0 + &
         sqrt(3d0)/2d0*(Vsigma-Vpi)*(cxy-cy)*GammaX + &
         LamISB*(sx + sy - sxy)*GammaY + &
         0.5d0*(Vsigma-Vpi)*(2*cx - cy - cxy)*GammaZ + &
         LamSOC*GammaS
    !
  end function hk_triang_pxpy




  ! !--------------------------------------------------------------------!
  ! !PURPOSE: Solve the topological Hamiltonian
  ! !--------------------------------------------------------------------!
  ! subroutine solve_hk_topological(sigma)
  !   integer                                :: i,j
  !   integer                                :: Npts
  !   complex(8),dimension(Nlso,Nlso)          :: sigma(Nlso,Nlso)
  !   real(8),dimension(:,:),allocatable     :: kpath
  !   !
  !   !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
  !   write(LOGfile,*)"Build H_TOP(k) BHZ along path:"
  !   !
  !   call set_sigmaBHZ()
  !   !
  !   Npts = 8
  !   Lk=(Npts-1)*Nkpath
  !   allocate(kpath(Npts,3))
  !   kpath(1,:)=kpoint_m1
  !   kpath(2,:)=kpoint_x2
  !   kpath(3,:)=kpoint_gamma
  !   kpath(4,:)=kpoint_x1
  !   kpath(5,:)=kpoint_m2
  !   kpath(6,:)=kpoint_r
  !   kpath(7,:)=kpoint_x3
  !   kpath(8,:)=kpoint_gamma
  !   call set_sigmaBHZ(sigma)
  !   call TB_solve_model(hk_bhz,Nlso,kpath,Nkpath,&
  !        colors_name=[red1,blue1,red1,blue1],&
  !        points_name=[character(len=20) :: "M","X","G","X1","A","R","Z","G"],&
  !        file="Eig_Htop.ed")
  !   if (usez) then
  !      write(*,*) "Z11=",Zmats(1,1)
  !      write(*,*) "Z22=",Zmats(2,2)
  !      write(*,*) "Z33=",Zmats(3,3)
  !      write(*,*) "Z44=",Zmats(4,4)
  !   endif
  ! end subroutine solve_hk_topological

















  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg,Nlso) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nlso,Nlso)               :: g
    integer                                     :: Nlso,i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nlso,Nlso)               :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so


end program



