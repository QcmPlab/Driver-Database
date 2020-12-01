program ed_bhz_3d
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                :: iloop,Lk,Nso
  logical                :: converged
  !Bath:
  integer                :: Nb
  real(8),allocatable    :: Bath(:),Bath_Prev(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:),sigmaBHZ(:,:)
  real(8),allocatable    :: Wtk(:)
  integer,allocatable    :: ik2ix(:),ik2iy(:),ik2iz(:)
  !variables for the model:
  integer                :: Nk,Nkpath
  real(8)                :: ez,mh,lambda
  real(8)                :: wmixing
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym,getpoles
  !
  type(finter_type)      :: finter_func
  !
  real(8),dimension(2)   :: Eout
  !
#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  !
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  !
  call parse_input_variable(ez,"ez",finput,default=1d0)
  call parse_input_variable(mh,"MH",finput,default=0.25d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  !
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(getpoles,"GETPOLES",finput,default=.false.)
  !
  call ed_read_input(trim(finput))
  !
  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb
  !
  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(SigmaBHZ(Nso,Nso))
  call set_sigmaBHZ()           !this set sigma_BHZ(0) to zero
  !
  !
  !
  if(getpoles)then
     call read_sigma(Sreal)
     call read_sigma(Smats)
     !Get 3d bands for U=0 Hamiltonian
     call build_Hk_path()
     !Get 3d Bands from Top. Hamiltonian
     call solve_hk_topological( so2j(Smats(:,:,:,:,1),Nso) )
     !Get Poles
     call get_poles

     stop
  endif
  !
  !
  !Buil the Hamiltonian on a grid or on path
  call build_hk(trim(hkfile))!
  !
  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nb))
  allocate(Bath_Prev(Nb))
  call ed_init_solver(bath)
  call set_hloc(j2so(bhzHloc))
  !
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     call build_z2_indices( so2j(Smats(:,:,:,:,1),Nso) ) !so2j transforms a Nspin:Nspin:Norb:Norb into a Nso:Nso matrix

     !Get the Local GF and the Weiss field/Delta function to be fitted
     call ed_get_gloc(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=1)
     call ed_get_weiss(Gmats,Smats,Delta,Hloc=j2so(bhzHloc),iprint=1)

     !Fit the new bath, starting from the old bath + the supplied delta
     call ed_chi2_fitgf(delta(1,1,:,:,:),bath,ispin=1)
     if(.not.spinsym)then
        call ed_chi2_fitgf(delta(2,2,:,:,:),bath,ispin=2)
     else
        call spin_symmetrize_bath(bath,save=.true.)
     endif

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
     Bath_Prev=Bath

     if(ED_MPI_ID==0)converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     if(ED_MPI_ID==0)call end_loop
  enddo
  !
  !Get kinetic energy:
  Eout = ed_kinetic_energy(Hk,Wtk,Smats)
  !Get 3d Bands from Top. Hamiltonian
  call solve_hk_topological( so2j(Smats(:,:,:,:,1),Nso) )
  !
#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif


contains



  !---------------------------------------------------------------------
  !PURPOSE: GET BHZ HAMILTONIAN IN THE FULL BZ
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
    real(8)                             :: kxgrid(Nk),kygrid(Nk),kzgrid(Nk)
    real(8)                             :: wm(Lmats),wr(Lreal),dw,n0(Nso)
    !
    !get H(k) and solve the non-interacting problem along a path in 3d:
    call build_hk_path()
    !
    !Get H(k) in the BZ:    
    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk**3
    if(ED_MPI_ID==0)write(*,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))
    kxgrid = kgrid(Nk);kygrid = kgrid(Nk);kzgrid = kgrid(Nk)
    Hk     = build_hk_model(hk_bhz,Nso,kxgrid,kygrid,kzgrid)
    wtk    = 1d0/Lk
    if(ED_MPI_ID==0.AND.present(file))then
       call write_hk_w90(file,Nso,&
            Nd=Norb,&
            Np=1,   &
            Nineq=1,&
            hk=Hk,  &
            kxgrid=kxgrid,&
            kygrid=kygrid,&
            kzgrid=kzgrid)
    endif
    !
    !Get the local part of H(k)
    if(allocated(bhzHloc))deallocate(bhzHloc)
    allocate(bhzHloc(Nso,Nso))
    bhzHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0d0
    if(ED_MPI_ID==0)call write_Hloc(bhzHloc)
    !
    !Build the local GF:
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i)+xmu,Hk(:,:,ik))/Lk
       enddo
       do i=1,Lreal
          Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu,Hk(:,:,ik))/Lk
       enddo
    enddo
    do iorb=1,Nso
       call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_iw.ed",wm,Gmats(iorb,iorb,:))
       call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_realw.ed",wr,&
            -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
    enddo
    !
  end subroutine build_hk




  !---------------------------------------------------------------------
  !PURPOSE: GET THE BHZ HAMILTONIAN ALONG A 3D PATH IN THE BZ
  !---------------------------------------------------------------------
  subroutine build_hk_path(kpath_)
    integer                                :: i,j
    integer                                :: Npts
    real(8),dimension(:,:),optional        :: kpath_
    real(8),dimension(:,:),allocatable     :: kpath
    character(len=64)                      :: file
    !
    !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
    !
    sigmaBHZ=zero
    if(present(kpath_))then
       if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) BHZ along a given path:"
       Npts = size(kpath_,1)
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,size(kpath_,2)))      
       kpath=kpath_
       file="Eig_path.nint"
    else
       if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) BHZ along the X-G-M-R-Z-A-G-Z path:"
       Npts = 8
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=[0,1,0]!X
       kpath(2,:)=[0,0,0]!G
       kpath(3,:)=[1,1,0]!M
       kpath(4,:)=[1,1,1]!R
       kpath(5,:)=[0,0,1]!Z
       kpath(6,:)=[1,0,1]!A
       kpath(7,:)=[0,0,0]!G
       kpath(8,:)=[0,0,1]!Z
       kpath=kpath*pi
       file="Eigenbands.nint"
    endif
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))
    Hk     = build_hk_model(hk_bhz,Nso,kpath,Nkpath)
    wtk    = 1d0/Lk
    if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_bhz,Nso,kpath,Lk,&
         colors_name=[character(len=20) :: 'red','blue','red','blue'],&
         points_name=[character(len=20) :: "X","G","M","R","Z","A","G","Z"],&
         file=reg(file))
  end subroutine build_hk_path


  subroutine solve_hk_topological(sigma)
    integer                                :: i,j
    integer                                :: Npts
    complex(8),dimension(Nso,Nso)          :: sigma(Nso,Nso)
    real(8),dimension(:,:),allocatable     :: kpath
    !
    !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
    if(ED_MPI_ID==0)then
       if(ED_MPI_ID==0)write(LOGfile,*)"Build H_TOP(k) BHZ along the X-G-M-R-Z-A-G-Z path:"
       !
       Npts = 8
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=[0,1,0]!X
       kpath(2,:)=[0,0,0]!G
       kpath(3,:)=[1,1,0]!M
       kpath(4,:)=[1,1,1]!R
       kpath(5,:)=[0,0,1]!Z
       kpath(6,:)=[1,0,1]!A
       kpath(7,:)=[0,0,0]!G
       kpath(8,:)=[0,0,1]!Z
       kpath=kpath*pi
       call set_sigmaBHZ(sigma)
       if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_bhz,Nso,kpath,Lk,&
            colors_name=[character(len=20) :: 'red','blue','red','blue'],&
            points_name=[character(len=20) :: "X","G","M","R","Z","A","G","Z"],&
            file="Eig_Htop.ed")
    endif
  end subroutine solve_hk_topological




  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  subroutine set_SigmaBHZ(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    sigmaBHZ = zero;if(present(sigma))sigmaBHZ=sigma
  end subroutine set_SigmaBHZ

  function hk_bhz(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !
    Hk          = zero
    Hk(1:2,1:2) = &
         (Mh - cos(kx) - cos(ky) - ez*cos(kz))*pauli_tau_z +&
         lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y
    Hk(3:4,3:4) = conjg( &
         (Mh-cos(-kx) - cos(-ky) - ez*cos(-kz))*pauli_tau_z +&
         lambda*sin(-kx)*pauli_tau_x + lambda*sin(-ky)*pauli_tau_y)
    Hk(1:2,3:4) = lambda*sin(kz)*pauli_tau_x
    Hk(3:4,1:2) = lambda*sin(kz)*pauli_tau_x
    !add the SigmaBHZ term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(SigmaBHZ)
    !
  end function hk_bhz




  !---------------------------------------------------------------------
  !PURPOSE: GET POLES ON THE REAL AXIS
  !---------------------------------------------------------------------
  subroutine get_poles
    USE IOFILE
    integer                                     :: i,j,ik,ix,iy
    integer                                     :: iorb,jorb
    integer                                     :: isporb,jsporb
    integer                                     :: ispin,jspin
    integer                                     :: iso,unit
    real(8),dimension(Nso)                      :: dzeta
    complex(8),dimension(Nso,Nso)               :: zeta,fg,gdelta,fgk
    complex(8),dimension(:,:,:,:,:),allocatable :: gloc,Sreal,Smats
    complex(8),dimension(:,:,:,:),allocatable   :: gk,gfoo,ReSmat
    complex(8)                                  :: iw
    complex(8),dimension(:,:),allocatable       :: detGiw
    real(8)                                     :: wr(Lreal),wm(Lmats)
    real(8),dimension(Lreal)                    :: Den
    real(8),dimension(:),allocatable            :: Ipoles,Xcsign,Iweight
    real(8),dimension(:,:),allocatable          :: Ipoles3d
    real(8),dimension(:,:),allocatable          :: Mpoles,Mweight
    real(8),dimension(:,:,:),allocatable        :: Mpoles3d
    integer                                     :: Linterval
    integer                                     :: count,Ninterval,maxNinterval,int
    real(8)                                     :: sign,sign_old
    real(8),dimension(:,:),allocatable :: kpath
    wr = linspace(wini,wfin,Lreal)
    wm = pi/beta*(2*arange(1,Lmats)-1)
    allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
    call read_sigma(Sreal)
    call read_sigma(Smats)
    !
    allocate(kpath(6,3))
    kpath(1,:)=[-0.125,-0.125,-0.125]!G-e<-R
    kpath(2,:)=[0,0,0]!G
    kpath(3,:)=[0.125,0.125,0.125]!G+e->R
    kpath(4,:)=[1-0.125,1-0.125,1-0.125]!R-e<-G
    kpath(5,:)=[1,1,1]!R
    kpath(6,:)=[1+0.125,1+0.125,1+0.125]!R+e->G
    kpath=kpath*pi
    call build_hk_path(kpath)
    !
    Linterval = 150 !Maximum number of allowed intervals to look for zeros&poles
    !
    allocate(Xcsign(0:Linterval))
    allocate(Ipoles(Lk),Iweight(Lk))
    allocate(Mpoles(Lk,Linterval),Mweight(Lk,Linterval))
    !
    !FINDING THE POLES:
    !assume \eps=0.d0 ==> the ImSigma(poles)=0 this condition should be automatically
    !verified at the pole from definition of the pole (the ImSigma at the pole is just
    !an artificial broadening of an otherwise delta function, whose position should be 
    !determined by ReSigma only.
    Ipoles=0.d0   
    Mpoles=0.d0
    write(LOGfile,*)"Solving for the poles..."
    maxNinterval=-1
    do ik=1,Lk
       do i=1,Lreal
          zeta(:,:) = (wr(i)+xmu)*eye(Nso) - dreal(so2j(Sreal(:,:,:,:,i),Nso))
          Den(i) = dreal((zeta(1,1) - Hk(1,1,ik))*(zeta(2,2) - Hk(2,2,ik))) - Hk(1,2,ik)*Hk(2,1,ik)
       enddo
       Xcsign(0)=0.d0
       count=0
       sign_old=sgn(Den(Lreal/2+1))
       do i=Lreal/2+1,Lreal
          sign=sgn(Den(i))
          if(sign*sign_old<1)then
             count=count+1
             if(count>Linterval)stop "Allocate Xcsign to a larger array."
             Xcsign(count)=wr(i)
          endif
          sign_old=sign
       enddo
       Ninterval=count
       if(count>maxNinterval)maxNinterval=count
       call init_finter(finter_func,wr,Den,3)
       do int=1,Ninterval
          Mpoles(ik,int) = fzero_brentq(det_poles,Xcsign(int-1),Xcsign(int))
          Mweight(ik,int)= get_weight(hk(:,:,ik)-so2j(Smats(:,:,:,:,1),Nso))
       enddo
       ipoles(ik) = fzero_brentq(det_poles,0.d0,wr(Lreal))
       iweight(ik)= get_weight(hk(:,:,ik)-so2j(Smats(:,:,:,:,1),Nso))
       call delete_finter(finter_func)
    enddo
    call splot("BHZpoles.ed",(/(ik-1,ik=1,Lk)/),ipoles(:),iweight(:))
    do int=1,maxNinterval
       unit=free_unit()
       open(unit,file="BHZpoles_int"//reg(txtfy(int))//".ed")
       if(any((Mpoles(:,int)/=0.d0)))then
          do ik=1,Lk
             if(Mpoles(ik,int)/=0.d0)write(unit,*)ik-1,Mpoles(ik,int),Mweight(ik,int)
          enddo
       endif
       close(unit)
    enddo
    !
  end subroutine get_poles

  function det_poles(w) result(det)
    real(8),intent(in) :: w
    real(8)            :: det
    det = finter(finter_func,w)
  end function det_poles

  function get_weight(hk) result(wt)
    complex(8),dimension(4,4) :: hk,foo
    real(8),dimension(4)      :: eigv
    real(8) :: wt
    foo = hk
    call matrix_diagonalize(foo,eigv)
    wt = sum(foo(:,1))
  end function Get_Weight







  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine build_z2_indices(sigma0)
    integer                                :: unit
    complex(8),dimension(Nso,Nso),optional :: sigma0(Nso,Nso)
    complex(8),dimension(Nso,Nso)          :: sigma0_
    !
    integer,dimension(4)                   :: z2
    !
    sigma0_=zero;if(present(sigma0))sigma0_=sigma0
    !Evaluate the Z2 index:
    !STRONG TI
    z2(1) = z2_number(reshape( [ [0,0,0] , [1,0,0] , [1,1,0] , [0,1,0] , [0,1,1] , [0,0,1] , [1,0,1] , [1,1,1] ] , [3,8] )*pi,sigma0_)
    !WEAK TI
    !K=1: n_1=1, n_2,3=0,1
    z2(2) = z2_number(reshape( [ [1,0,0] , [1,1,0] , [1,1,1] , [1,0,1] ] , [3,4])*pi,sigma0_)
    !K=2: n_2=1, n_1,2=0,1
    z2(3) = z2_number(reshape( [ [0,1,0] , [0,1,1] , [1,1,1] , [1,1,0] ] , [3,4])*pi,sigma0_)
    !k=3: n_3=1, n_1,2=0,1
    z2(4) = z2_number(reshape( [ [0,0,1] , [0,1,1] , [1,1,1] , [1,0,1] ] , [3,4])*pi,sigma0_)
    unit=free_unit()
    open(unit,file="z2_invariant.ed")
    write(unit,*)z2
    close(unit)
  end subroutine build_z2_indices

  function z2_number(ktrims,sigma0) result(z2)
    real(8),dimension(:,:),intent(in)            :: ktrims
    complex(8),dimension(Nso,Nso)                :: sigma0
    complex(8),dimension(Nso,Nso,size(Ktrims,2)) :: Htrims
    complex(8),dimension(size(Ktrims,2))         :: Delta
    integer                                      :: z2
    integer                                      :: i,j,Ntrim,itrim,Nocc
    !
    Ntrim=size(Ktrims,2)
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_bhz(Ktrims(:,itrim),Nso) + sigma0(:,:)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    !
    z2=product(Delta(:))
    if(z2>0)then
       z2=0
    else
       z2=1
    end if
    !
  end function z2_number

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

  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k = (iw + xmu)*zeye(4)-hk
    call inv(g0k)
  end function inverse_g0k



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


end program ed_bhz_3d



