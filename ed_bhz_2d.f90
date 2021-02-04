program ed_bhz
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE SF_MPI
  implicit none
  integer                                     :: iloop,Lk,Nso
  logical                                     :: converged
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,print_mode
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:),Weiss_(:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable                      :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  complex(8),allocatable,dimension(:)         :: Gtest
  !hamiltonian input:
  complex(8),allocatable                      :: Hk(:,:,:),bhzHloc(:,:),sigmaBHZ(:,:),Zmats(:,:)
  real(8),allocatable                         :: Wtk(:)
  real(8),allocatable                         :: kxgrid(:),kygrid(:)
  integer,allocatable                         :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: mh,lambda,wmixing,akrange,rh
  character(len=16)                           :: finput
  character(len=32)                           :: hkfile
  logical                                     :: spinsym,usez,mixG0
  !
  real(8),dimension(2)                        :: Eout
  real(8),allocatable                         :: dens(:)
  complex(8),dimension(4,4)                   :: Gamma1,Gamma2,Gamma5,GammaS
  complex(8),dimension(4,4)                   :: GammaN,GammaTz,GammaSz,GammaRz
  complex(8),dimension(4,4)                   :: GammaE0,GammaEx,GammaEy,GammaEz
  real(8),dimension(:),allocatable            :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable :: Hsym_basis
  !MPI Vars:
  integer                                     :: irank,comm,rank,size2,ierr
  logical                                     :: master

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  size2 = get_Size_MPI(comm)
  master = get_Master_MPI(comm)


  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(usez,"USEZ",finput,default=.false.)
  !
  call ed_read_input(trim(finput),comm)
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
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(dens(Norb))
  allocate(SigmaBHZ(Nso,Nso))
  allocate(Zmats(Nso,Nso))

  gamma1=kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2=kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaS=kron_pauli( pauli_sigma_z, pauli_tau_0)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z )
  !
  gammaN=kron_pauli( pauli_sigma_0, pauli_tau_0 )
  gammaTz=kron_pauli( pauli_sigma_0, pauli_tau_z )
  gammaSz=kron_pauli( pauli_sigma_z, pauli_tau_0 )
  gammaRz=kron_pauli( pauli_sigma_z, pauli_tau_z )
  !
  gammaE0=kron_pauli( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron_pauli( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron_pauli( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron_pauli( pauli_sigma_z, pauli_tau_x )

  !Buil the Hamiltonian on a grid or on  path
  call set_sigmaBHZ()
  call build_hk(trim(hkfile))

  print_mode=1
  if(ed_mode=="nonsu2")print_mode=4

  !Setup solver
  if(bath_type=="replica")then
     !Setup HLOC symmetries
     allocate(lambdasym_vector(4))
     allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,4))
     !
     lambdasym_vector(1)=Mh
     lambdasym_vector(2:)=sb_field
     Hsym_basis(:,:,:,:,1)=j2so(Gamma5)
     Hsym_basis(:,:,:,:,2)=j2so(GammaE0)
     Hsym_basis(:,:,:,:,3)=j2so(GammaEx)
     Hsym_basis(:,:,:,:,4)=j2so(GammaEz)
     ! Hsym_basis(:,:,:,:,5)=j2so(GammaEy)



     call ed_set_Hloc(Hsym_basis,lambdasym_vector)
     Nb=ed_get_bath_dimension(Hsym_basis)
     allocate(Bath(Nb))
     allocate(Bath_(Nb))
     call ed_init_solver(comm,bath)    
  else     
     Nb=ed_get_bath_dimension()
     allocate(Bath(Nb))
     allocate(Bath_(Nb))
     call ed_init_solver(comm,bath,Hloc=j2so(bhzHloc))
  endif


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     call ed_get_dens(dens)


     !Get GLOC:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=print_mode)


     !Update WeissField:
     call dmft_self_consistency(Gmats,Smats,Weiss,j2so(bhzHloc),cg_scheme)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=print_mode)


     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif

     !Fit the new bath, starting from the old bath + the supplied Weiss
     select case(ed_mode)
     case default
        stop "ed_mode!=Normal/Nonsu2"
     case("normal")
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
        if(.not.spinsym)then
           call ed_chi2_fitgf(comm,Weiss,bath,ispin=2)
        else
           call ed_spin_symmetrize_bath(bath,save=.true.)
        endif
     case("nonsu2")
        call ed_chi2_fitgf(comm,Weiss,bath)
     end select

     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     !

     !Check convergence (if required change chemical potential)
     ! Gtest=zero
     ! do ispin=1,Nspin
     !    do iorb=1,Norb
     !       Gtest=Gtest+Weiss(ispin,ispin,iorb,iorb,:)/Norb/Nspin
     !    enddo
     ! enddo
     Gtest=Weiss(1,1,1,1,:)
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo




  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=print_mode)

  call dmft_kinetic_energy(Hk,Smats)

  call solve_hk_topological(so2j(Smats(:,:,:,:,1)))


  call finalize_MPI()



contains


  !---------------------------------------------------------------------
  !PURPOSE: GET BHZ HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky    
    integer                             :: iorb,jorb
    integer                             :: isporb,jsporb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Nso,Nso,Lmats) :: Gmats
    complex(8),dimension(Nso,Nso,Lreal) :: Greal
    real(8)                             :: wm(Lmats),wr(Lreal),dw
    !
    call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
    !
    if(master)write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk**2
    if(master)write(*,*)"# of k-points     :",Lk
    if(master)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
    allocate(wtk(Lk))
    !
    call TB_build_model(Hk,hk_bhz,Nso,[Nk,Nk])
    wtk = 1d0/Lk
    if(master.AND.present(file))then
       call TB_write_hk(Hk,trim(file),&
            Nlat=1,&
            Nspin=1,&
            Norb=Norb,&
            Nkvec=[Nk,Nk])
    endif
    allocate(bhzHloc(Nso,Nso))
    bhzHloc = zero
    bhzHloc = sum(Hk,dim=3)/Lk
    where(abs(dreal(bhzHloc))<1d-6)bhzHloc=zero
    if(master)  call TB_write_Hloc(bhzHloc)
  end subroutine build_hk






  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  subroutine set_SigmaBHZ(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    sigmaBHZ = zero;if(present(sigma))sigmaBHZ=sigma
  end subroutine set_SigmaBHZ


  !--------------------------------------------------------------------!
  !PURPOSE: Solve the topological Hamiltonian
  !--------------------------------------------------------------------!
  subroutine solve_hk_topological(sigma)
    integer                                :: i,j
    integer                                :: Npts
    complex(8),dimension(Nso,Nso)          :: sigma(Nso,Nso)
    real(8),dimension(:,:),allocatable     :: kpath
    !
    if(master)then
       !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
       write(LOGfile,*)"Build H_TOP(k) BHZ along path:"
       !
       call set_sigmaBHZ()
       !
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,2))
       kpath(1,:)=kpoint_gamma
       kpath(2,:)=kpoint_x1
       kpath(3,:)=kpoint_m1
       kpath(4,:)=kpoint_gamma
       call set_sigmaBHZ(sigma)
       call TB_solve_model(hk_bhz,Nso,kpath,Nkpath,&
            colors_name=[red,blue,red,blue],&
            points_name=[character(len=20) :: "{\Symbol G}","X","M","{\Symbol G}"],&
            file="Eig_Htop.ed")
    endif
  end subroutine solve_hk_topological




  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_bhz(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: ek,kx,ky
    integer                   :: ii
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    kx=kvec(1)
    ky=kvec(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
    ! !add the SigmaBHZ term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(SigmaBHZ)
    if (usez) then
       Zmats=zero
       do ii=1,Nso
          Zmats(ii,ii)  = 1.d0/abs( 1.d0 +  abs(dimag(sigmaBHZ(ii,ii))/(pi/beta)) )
       end do
       Hk = matmul(Zmats,Hk)
    endif
  end function hk_bhz











  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine read_sigma(sigma)
    complex(8)        :: sigma(:,:,:,:,:)
    integer           :: iorb,ispin,i,L,unit
    real(8)           :: reS(Nspin),imS(Nspin),ww
    character(len=20) :: suffix
    if(size(sigma,1)/=Nspin)stop "read_sigma: error in dim 1. Nspin"
    if(size(sigma,3)/=Norb)stop "read_sigma: error in dim 3. Norb"
    L=size(sigma,5);print*,L
    if(L/=Lmats.AND.L/=Lreal)stop "read_sigma: error in dim 5. Lmats/Lreal"
    do iorb=1,Norb
       unit=free_unit()
       if(L==Lreal)then
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_realw.ed"
       elseif(L==Lmats)then
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed"
       endif
       write(*,*)"read from file=","impSigma"//reg(suffix)
       open(unit,file="impSigma"//reg(suffix),status='old')
       do i=1,L
          read(unit,"(F26.15,6(F26.15))")ww,(imS(ispin),reS(ispin),ispin=1,Nspin)
          forall(ispin=1:Nspin)sigma(ispin,ispin,iorb,iorb,i)=dcmplx(reS(ispin),imS(ispin))
       enddo
       close(unit)
    enddo
  end subroutine read_sigma


  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
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
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
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


end program ed_bhz



