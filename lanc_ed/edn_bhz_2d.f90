program edn_bhz
  USE LANC_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE SF_MPI
  implicit none
  integer                :: iloop,Lk,Nso
  logical                :: converged
  !Bath:
  integer                :: Nb
  real(8),allocatable    :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:),sigmaBHZ(:,:),Zmats(:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:)
  integer,allocatable    :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                :: Nk,Nkpath
  real(8)                :: mh,lambda,wmixing,akrange,rh
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym,getak,getdeltaw,getpoles,usez
  type(finter_type)      :: finter_func
  !
  real(8),dimension(2)   :: Eout
  real(8),allocatable    :: dens(:)
  !MPI Vars:
  integer                :: irank,comm,rank,size2
  logical                :: master

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
  call parse_input_variable(rh,"RH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
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
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(dens(Norb))
  allocate(SigmaBHZ(Nso,Nso))
  allocate(Zmats(Nso,Nso))

  !Buil the Hamiltonian on a grid or on  path
  call set_sigmaBHZ()
  call build_hk(trim(hkfile))

  !Setup solver
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(comm,bath,Hloc=j2so(bhzHloc))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath)
     !
     call ed_get_sigma_matsubara(Smats)
     call dmft_gloc_matsubara(comm,Hk,Wtk,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     if(master)then
        call dmft_self_consistency(Gmats,Smats,Delta,j2so(bhzHloc),cg_scheme)
        call dmft_print_gf_matsubara(Delta,"Weiss",iprint=1)
     endif
     call Bcast_MPI(comm,Delta)

     !Fit the new bath, starting from the old bath + the supplied delta
     call ed_chi2_fitgf(comm,delta,bath,ispin=1)
     if(.not.spinsym)then
        call ed_chi2_fitgf(comm,delta,bath,ispin=2)
     else
        call ed_spin_symmetrize_bath(bath,save=.true.)
     endif
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath
     !
     if(master)then
        converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
        call ed_get_dens(dens)     
        if(nread/=0.d0)call ed_search_variable(xmu,sum(dens),converged)
     endif
     call Bcast_MPI(comm,bath)
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)

     if(master)call end_loop
  enddo

  call ed_get_sigma_realaxis(Sreal)
  call dmft_gloc_realaxis(comm,Hk,Wtk,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)


  call dmft_kinetic_energy(comm,Hk,Wtk,Smats)

  call solve_hk_topological(so2j(Smats(:,:,:,:,1),Nso))

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
    if(master)write(LOGfile,*)"# of k-points     :",Lk
    if(master)write(LOGfile,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
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
    bhzHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0d0
    if(master)  call TB_write_Hloc(bhzHloc)
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
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="Eig_path.nint"
    else
       if(master)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_Gamma
       kpath(2,:)=kpoint_M1
       kpath(3,:)=kpoint_X1
       kpath(4,:)=kpoint_Gamma
       file="Eigenbands.nint"
    endif
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))
    call set_sigmaBHZ()
    call TB_build_model(Hk,hk_bhz,Nso,kpath,Nkpath)
    wtk = 1d0/Lk
    if(master)  call TB_Solve_model(hk_bhz,Nso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
  end subroutine build_hk_GXMG


  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  !
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
    !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
    write(LOGfile,*)"Build H_TOP(k) BHZ along path:"
    !
    call set_sigmaBHZ()
    !
    Npts = 8
    Lk=(Npts-1)*Nkpath
    allocate(kpath(Npts,3))
    kpath(1,:)=kpoint_m1
    kpath(2,:)=kpoint_x2
    kpath(3,:)=kpoint_gamma
    kpath(4,:)=kpoint_x1
    kpath(5,:)=kpoint_m2
    kpath(6,:)=kpoint_r
    kpath(7,:)=kpoint_x3
    kpath(8,:)=kpoint_gamma
    call set_sigmaBHZ(sigma)
    call TB_solve_model(hk_bhz,Nso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=20) :: "M","X","G","X1","A","R","Z","G"],&
         file="Eig_Htop.ed")
    if (usez) then
       write(LOGfile,*) "Z11=",Zmats(1,1)
       write(LOGfile,*) "Z22=",Zmats(2,2)
       write(LOGfile,*) "Z33=",Zmats(3,3)
       write(LOGfile,*) "Z44=",Zmats(4,4)
    endif
  end subroutine solve_hk_topological




  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_bhz(kvec,N) result(hk)
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky
    integer                   :: N,ii
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    kx=kvec(1)
    ky=kvec(2)
    Hk          = zero
    Hk(1:2,1:2) = hk_bhz2x2(kx,ky)
    Hk(3:4,3:4) = conjg(hk_bhz2x2(-kx,-ky))
    ! Hk(1,4) = -delta ; Hk(4,1)=-delta
    ! Hk(2,3) =  delta ; Hk(3,2)= delta
    Hk(1,3) = xi*rh*(sin(kx)-xi*sin(ky))
    Hk(3,1) =-xi*rh*(sin(kx)+xi*sin(ky))
    !add the SigmaBHZ term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(SigmaBHZ)
    !
    if (usez) then
       Zmats=zero
       do ii=1,Nso
          Zmats(ii,ii)  = 1.d0/abs( 1.d0 +  abs(dimag(sigmaBHZ(ii,ii))/(pi/beta)) )
       end do
       Hk = matmul(Zmats,Hk)
    endif
  end function hk_bhz

  function hk_bhz2x2(kx,ky) result(hk)
    real(8)                   :: kx,ky,epsik
    complex(8),dimension(2,2) :: hk
    epsik   = cos(kx)+cos(ky)
    hk(1,1) = mh - epsik
    hk(2,2) =-mh + epsik
    hk(1,2) = lambda*(sin(kx)-xi*sin(ky))
    hk(2,1) = lambda*(sin(kx)+xi*sin(ky))
  end function hk_bhz2x2











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
       write(LOGfile,*)"read from file=","impSigma"//reg(suffix)
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


  function so2j(fg,Nso) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nso,Nso)               :: g
    integer                                     :: Nso,i,j,iorb,jorb,ispin,jspin
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
    complex(8),dimension(Nso,Nso)               :: fg
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


end program edn_bhz



