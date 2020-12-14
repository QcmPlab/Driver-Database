!  1UP----2DW
program ed_bhz_afm2_2d
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                       :: ip,iloop,ilat,ineq,Lk,Nso,Nlso,ispin,iorb
  logical                                       :: converged
  integer                                       :: Nineq,Nlat
  !Bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath_ineq(:,:),Bath_prev(:,:)
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_ineq
  !Hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk ![Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nk]
  complex(8),allocatable,dimension(:,:)         :: bhzHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc_ineq
  real(8),allocatable,dimension(:)              :: Wtk
  !variables for the model:
  integer                                       :: Nktot,Nkx,Nkpath,unit
  real(8)                                       :: mh,lambda,wmixing
  character(len=16)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: waverage,spinsym,fullsym
  !Dirac matrices:
  complex(8),dimension(4,4)                     :: Gamma1,Gamma2,Gamma3,Gamma4,Gamma5
  integer                                       :: comm,rank
  logical                                       :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.conf')
  call parse_input_variable(nkx,"NKX",finput,default=25)
  call parse_input_variable(Nkpath,"Nkpath",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(waverage,"WAVERAGE",finput,default=.false.)
  call parse_input_variable(spinsym,"spinsym",finput,default=.false.)
  call parse_input_variable(fullsym,"fullsym",finput,default=.true.)
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call ed_read_input(trim(finput),comm)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  Nlat=2                      !number of independent sites, 4 for AFM ordering
  Nineq=Nlat
  Nso=Nspin*Norb
  Nlso=Nlat*Nso
  if(Norb/=2)stop "Norb != 2"
  if(Nspin/=2)stop "Nspin != 2"
  if(Nso/=4)stop "Nso != 4"
  if(Nlso/=8)stop "Nlso != 8"

  if(fullsym)then
     Nineq=1
     write(*,*)"Using Nineq sites=",Nineq
     open(free_unit(unit),file="symmetries.used")
     write(unit,*)"Symmetries used are:"
     write(unit,*)"(site=2,l,s)=(site=1,l,-s)"
     close(unit)
  endif

  if(spinsym)sb_field=0.d0

  Gamma1 = kron_pauli(pauli_z,pauli_x)
  Gamma2 =-kron_pauli(pauli_0,pauli_y)
  Gamma3 = kron_pauli(pauli_x,pauli_x)
  Gamma4 = kron_pauli(pauli_y,pauli_x)
  Gamma5 = kron_pauli(pauli_0,pauli_z)

  !Allocate Weiss Field:
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))


  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(bhzHloc,Nlat,Nspin,Norb)
  do ip=1,Nineq
     Hloc_ineq(ip,:,:,:,:) = Hloc(ip,:,:,:,:)
  enddo


  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  call ed_init_solver(Bath_ineq)
  do ip=1,Nineq
     call break_symmetry_bath(Bath_ineq(ip,:),sb_field,(-1d0)**(ip+1))
  enddo


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master) call start_loop(iloop,nloop,"DMFT-loop")
     !
     !
     call ed_solve(Comm,Bath_ineq,Hloc_ineq,iprint=1)
     call ed_get_sigma_matsubara(Smats_ineq,Nineq)
     call ed_get_sigma_real(Sreal_ineq,Nineq)
     do ip=1,Nineq
        Smats(ip,:,:,:,:,:) = Smats_ineq(ip,:,:,:,:,:)
        Sreal(ip,:,:,:,:,:) = Sreal_ineq(ip,:,:,:,:,:)
     enddo
     if(fullsym)then
        do ispin=1,2
           Smats(2,ispin,ispin,:,:,:)=Smats(1,3-ispin,3-ispin,:,:,:)
           Sreal(2,ispin,ispin,:,:,:)=Sreal(1,3-ispin,3-ispin,:,:,:)
        enddo
     endif
     !
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats,iprint=4) !tridiag option off
     do ip=1,Nineq
        Gmats_ineq(ip,:,:,:,:,:) = Gmats(ip,:,:,:,:,:)
     enddo
     !
     !
     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=4)
     else
        call dmft_delta(Comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=4)
     endif
     !
     !
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     if(spinsym)then
        call spin_symmetrize_bath(Bath_ineq,save=.true.)
     else
        call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     endif
     !
     !
     ! Mixing:
     if(iloop>1)Bath_ineq = wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq
     !
     ! Convergence
     if(master)converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     call Bcast_MPI(Comm,converged)
     !
     if(master)call end_loop
  enddo


  call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal,iprint=4)


  call finalize_MPI()


  print*,"Bravo"

contains




  !--------------------------------------------------------------------!
  !PURPOSE: BUILD THE H(k) FOR THE BHZ-AFM MODEL.
  !--------------------------------------------------------------------!
  subroutine build_hk(file)
    character(len=*)                        :: file
    integer                                 :: Npts
    integer                                 :: i,j,k,ik,iorb,jorb
    integer                                 :: ix,iy,iz
    real(8)                                 :: kx,ky,kz
    real(8),dimension(:),allocatable        :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable      :: kpath
    real(8),dimension(2)                    :: bk1,bk2,kvec
    real(8)                                 :: n(Nlso)
    complex(8)                              :: w
    complex(8)                              :: Gmats(Nlso,Nlso,Lmats),Greal(Nlso,Nlso,Lreal)
    complex(8)                              :: Smats(Nlso,Nlso,Lmats),Sreal(Nlso,Nlso,Lreal)
    !
    Nktot=Nkx*Nkx
    allocate(Hk(Nlso,Nlso,Nktot))
    allocate(kxgrid(Nkx),kygrid(Nkx))
    write(LOGfile,*)"Build H(k) AFM2-BHZ 2d:"
    write(LOGfile,*)"Using Nk_total="//txtfy(Nktot)
    !
    ! kxgrid = kgrid(Nkx)
    ! Hk = build_hk_model(hk_model,Nlso,kxgrid,kxgrid,[0d0])
    ! call write_hk_w90(trim(file),Nlso,&
    !      Nd=Nso,&
    !      Np=0,   &
    !      Nineq=2,&
    !      hk=Hk,  &
    !      kxgrid=kxgrid,&
    !      kygrid=kygrid,&
    !      kzgrid=[0d0])
    bk1 = pi*[1,-1]
    bk2 = 2*pi*[0,1]
    ik=0
    do iy=1,Nkx
       ky = dble(iy-1)/Nkx
       do ix=1,Nkx
          kx=dble(ix-1)/Nkx
          ik=ik+1
          kvec = kx*bk1 + ky*bk2
          Hk(:,:,ik) = hk_model(kvec,Nlso)
       enddo
    enddo
    !
    !
    allocate(bhzHloc(Nlso,Nlso))
    bhzHloc = sum(Hk(:,:,:),dim=3)/Nktot
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0.d0
    !
    allocate(Wtk(Nktot))
    Wtk=1.d0/dble(Nktot)
    !
    !
    ! Gmats=zero
    ! Greal=zero
    ! Smats=zero
    ! Sreal=zero
    ! call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=1)
    ! call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal,iprint=1)
    ! do iorb=1,Nlso
    !    n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
    ! enddo
    ! !
    ! open(10,file="observables.nint")
    ! write(10,"(20F20.12)")(n(iorb),iorb=1,Nlso),sum(n)
    ! close(10)
    ! write(*,"(A,20F14.9)")"Occupations =",(n(iorb),iorb=1,Nlso),sum(n)/size(n)
    !
    !
    !solve along the standard path in the 2D BZ.
    Npts=4
    allocate(kpath(Npts,3))
    kpath(1,:)=kpoint_Gamma
    kpath(2,:)=kpoint_X1
    kpath(3,:)=kpoint_M1
    kpath(4,:)=kpoint_Gamma
    call solve_Hk_along_BZpath(Hk_model,Nlso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1, red1,blue1,red1,blue1],&
         points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
         file="Eigenbands_afm2.nint")
  end subroutine build_hk





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  !+-----------------------------------------------------------------------------+!
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    complex(8),dimension(Nso,Nso) :: M
    complex(8),dimension(Nso,Nso) :: tx,ty,thx,thy
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = kpoint(1)
    ky = kpoint(2)
    !
    !
    M  = Mh*Gamma5
    tx = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    thx= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    !
    ty = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    thy= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    !
    ! H2 =  | m1                       tx + tx^+.e^i.2.kx + ty^+.e^i.(kx+ky) + ty^+.e^i.(kx-ky) |
    !       | tx^+ + tx.e^-i.2.kx + ty.e^-i.(kx+ky)+ ty^+.e^-i.(kx-ky)          m2              |
    !
    hk(1:4,1:4)    = M
    hk(1:4,5:8)    = tx  + thx*exp(xi*2*kx) + thy*exp(xi*(kx+ky)) + ty*exp(xi*(kx-ky))
    !
    hk(5:8,1:4)    = thx + tx*exp(-xi*2*kx) + ty*exp(-xi*(kx+ky)) + thy*exp(-xi*(kx-ky))
    hk(5:8,5:8)    = M
    !
  end function hk_model






end program ed_bhz_afm2_2d



