program lancED
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Nb,Le,Nso,iorb
  logical                                     :: converged
  real(8),dimension(5)                        :: Wbethe,Dbethe
  !Bath:
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !  
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: H0
  real(8),dimension(:),allocatable            :: de,dens
  !
  real(8),dimension(:),allocatable            :: Wband
  !
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal,Weiss_
  complex(8),allocatable,dimension(:) :: Gtest
  character(len=16)                           :: finput
  real(8)                                     :: wmixing
  !
  integer                                     :: comm,rank
  logical                                     :: master
  logical :: betheSC,wGimp,mixG0,symOrbs

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wbethe,"WBETHE",finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Dbethe,"DBETHE",finput,default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(betheSC,"BETHESC",finput,default=.false.)
  call parse_input_variable(wGimp,"wGimp",finput,default=.false.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.)
  !
  call ed_read_input(trim(finput),comm)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  if(Nspin/=1.OR.Norb>5)stop "Wrong setup from input file: Nspin=1 OR Norb>5"
  Nso=Nspin*Norb


  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(Wband(Nso))
  allocate(H0(Nso))
  allocate(de(Nso))
  !
  Wband = Wbethe(:Norb)
  H0    = Dbethe(:Norb)
  do iorb=1,Norb
     Ebands(iorb,:) = linspace(-Wband(iorb),Wband(iorb),Le,mesh=de(iorb))
     Dbands(iorb,:) = dens_bethe(Ebands(iorb,:),Wband(iorb))*de(iorb)
  enddo
  if(master)call TB_write_Hloc(one*diag(H0))


  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  Hloc(1,1,:,:)=diag(H0)


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bath_(Nb))
  call ed_init_solver(comm,bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(comm,Ebands,Dbands,H0,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     !
     call dmft_gloc_realaxis(Comm,Ebands,Dbands,H0,Greal,Sreal)
     if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
     !
     !Get the Weiss field/Delta function to be fitted
     if(.not.betheSC)then
        call dmft_self_consistency(comm,Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)
     else
        if(wGimp)call ed_get_gimp_matsubara(Gmats)
        call dmft_self_consistency(comm,Gmats,Weiss,Hloc,SCtype=cg_scheme,wbands=Wband)
     endif
     if(master)call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
     !
     !
     !
     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(symOrbs)then
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1,iorb=1)
        call orb_equality_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
     endif
     !
     !MIXING:
     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     !
     !Check convergence (if required change chemical potential)
     if(master)then
        Gtest=zero
        do iorb=1,Norb
           Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
        enddo
        converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
        if(nread/=0d0)then
           call ed_get_dens(dens)
           call ed_search_variable(xmu,sum(dens),converged)
        endif
     endif
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo


  call dmft_gloc_realaxis(Comm,Ebands,Dbands,H0,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
  call dmft_kinetic_energy(comm,Ebands,Dbands,H0,Smats(1,1,:,:,:))

  call finalize_MPI()


end program lancED

  






  ! !+----------------------------------------+
  ! subroutine get_delta_bethe
  !   integer                     :: i,j,iorb
  !   complex(8)                  :: iw,zita,g0loc
  !   complex(8),dimension(Lmats) :: gloc,sigma,Tiw
  !   complex(8),dimension(Lreal) :: grloc
  !   real(8)                     :: wm(Lmats),wr(Lreal),tau(0:Lmats),C0,C1,n0
  !   real(8),dimension(0:Lmats)  :: sigt,gtau,Ttau
  !   real(8),dimension(3)  :: Scoeff
  !   wm = pi/beta*(2*arange(1,Lmats)-1)
  !   wr = linspace(wini,wfin,Lreal)
  !      do i=1,Lmats
  !         iw = xi*wm(i)
  !         zita    = iw + xmu - impSmats(1,1,1,1,i)
  !         gloc(i) = gfbethe(wm(i),zita,Wband)
  !         if(cg_scheme=='weiss')then
  !            delta(i)= one/(one/gloc(i) + impSmats(1,1,1,1,i))
  !         else
  !            delta(i)= iw + xmu - impSmats(1,1,1,1,i) - one/gloc(i)
  !         endif
  !      enddo
  !      do i=1,Lreal
  !         iw=cmplx(wr(i),eps)
  !         zita     = iw + xmu - impSreal(1,1,1,1,i)
  !         grloc(i) = gfbether(wr(i),zita,Wband)
  !      enddo
  !      if(ED_MPI_ID==0)then
  !         call splot("Gloc_iw.ed",wm,gloc)
  !         call splot("Delta_iw.ed",wm,delta)
  !         call splot("Gloc_realw.ed",wr,-dimag(grloc)/pi,dreal(grloc))
  !      endif
  !      ! tau(0:) = linspace(0.d0,beta,Lmats+1)
  !      ! C0=Uloc(1)*(ed_dens_up(1)-0.5d0)
  !      ! C1=Uloc(1)**2*ed_dens_up(1)*(1.d0-ed_dens_dw(1))
  !      ! Tiw=dcmplx(C0,-C1/wm)
  !      ! call splot("Tail_iw.ed",wm,Tiw)
  !      ! Ttau = -C1/2.d0
  !      ! Sigma = impSmats(1,1,1,1,:)  - Tiw
  !      ! call fftgf_iw2tau(Sigma,Sigt(0:),beta,notail=.true.)
  !      ! Sigt=Sigt + Ttau
  !      ! call splot("Sigma_tau.ed",tau,sigt)
  !      ! Sigt=Sigt !+ Ttau
  !      ! call fftgf_tau2iw(sigt(0:),sigma,beta)
  !      ! Sigma=Sigma !+ Tiw
  !      ! call splot("Sigma_iw.ed",wm,sigma)
  !      ! call fftgf_iw2tau(gloc,gtau(0:),beta)
  !      ! call splot("Gloc_tau.ed",tau(0:),gtau(0:))
  ! end subroutine get_delta_bethe
  ! !+----------------------------------------+






  ! !<DEBUG:
  ! ! subroutine get_ed_energy(Lk) 
  ! !   integer               :: Lk
  ! !   real(8),dimension(Lk) :: ek
  ! !   real(8)               :: de
  ! !   real(8),dimension(Lk) :: Wtk
  ! !   ek  = linspace(-Wband,Wband,Lk,mesh=de)
  ! !   Wtk = dens_bethe(ek,wband)*de
  ! !   call ed_kinetic_energy(impSmats(1,1,1,1,:),ek,wtk)
  ! ! end subroutine get_ed_energy


  ! function get_energy(Lk) result(H0)
  !   integer                     :: Lk
  !   complex(8),dimension(Lk)    :: Hk
  !   complex(8),dimension(Lmats) :: Sigma
  !   real(8)                     :: H0
  !   real(8),dimension(Lk)       :: Wtk
  !   real(8)                     :: Tail0,Tail1
  !   real(8)                     :: Sigma_HF,Ak,Bk
  !   complex(8)                  :: Ck,Dk,Zk
  !   complex(8)                  :: Zeta,Gk,Tk
  !   integer                     :: i,ik,iorb
  !   real(8),dimension(Lmats)    :: wm
  !   real(8)                     :: de
  !   !
  !   wm = pi/beta*dble(2*arange(1,Lmats)-1)
  !   !
  !   Hk  = one*linspace(-Wband,Wband,Lk,mesh=de)
  !   Wtk = dens_bethe(dreal(Hk(:)),wband)*de
  !   Sigma = impSmats(1,1,1,1,:) 
  !   Sigma_HF = dreal(Sigma(Lmats))
  !   !
  !   H0=0.d0
  !   do ik=1,Lk
  !      Ak = Hk(ik)
  !      Bk =-Hk(ik) - Sigma_hf
  !      do i=1,Lmats
  !         Gk = one/(xi*wm(i) + xmu - Hk(ik) - Sigma(i) )
  !         Tk = one/(xi*wm(i)) - Bk/(xi*wm(i))**2
  !         Ck = Ak*(Gk - Tk)
  !         H0 = H0 + Ck*Wtk(ik)
  !      enddo
  !   enddo
  !   H0=H0/beta*4d0
  !   !
  !   Tail0=zero
  !   Tail1=zero
  !   do ik=1,Lk
  !      Ak= Hk(ik)
  !      Bk =-Hk(ik) - Sigma_hf
  !      Ck= Ak*Bk
  !      Tail0 = Tail0 + 0.5d0*Ak*Wtk(ik)
  !      Tail1 = Tail1 + 0.25d0*Ck*Wtk(ik)*beta
  !   enddo
  !   Tail0=2d0*Tail0
  !   Tail1=2d0*Tail1
  !   H0 = H0 + Tail0 + Tail1
  ! end function get_energy
  ! !>DEBUG




