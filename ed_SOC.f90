program ed_SOC
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  !
  !#########   VARIABLEs DECLARATION   #########
  !
  integer                                        :: iloop,i,j
  integer                                        :: Nlat,ilat
  integer                                        :: Nso,io,jo
  integer                                        :: iorb,jorb,ispin,jspin
  logical                                        :: converged
  real(8)                                        :: wmixing
  character(len=60)                              :: finput
  character(len=32)                              :: hkfile
  !Mpi:
  integer                                        :: comm,rank,ier
  logical                                        :: master
  !Bath:
  integer                                        :: Nb,unit
  real(8),allocatable                            :: Bath(:)
  real(8),allocatable                            :: Bath_old(:)
  !Local functions:
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Sreal
  complex(8),allocatable,dimension(:)            :: Saux_wr,Saux_wm
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Greal
  !Weiss&Hybridization functions
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Weiss_old
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Delta
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Delta_old
  !Hmiltonian input:
  integer                                        :: Nk
  integer                                        :: Nkpath
  complex(8),allocatable,dimension(:,:,:)        :: Hk
  real(8),allocatable,dimension(:)               :: Wtk
  complex(8),allocatable,dimension(:,:)          :: d_t2g_Hloc_so
  complex(8),allocatable,dimension(:,:,:,:)      :: d_t2g_Hloc_nn
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Sigma_correction
  !Variables for the model:
  real(8)                                        :: soc,ivb
  !custom variables for rotations:
  logical                                        :: surface
  logical                                        :: Hk_test
  logical                                        :: rotateG0loc
  complex(8),allocatable,dimension(:,:,:,:,:)    :: impG
  !custom variables for convergence test:
  complex(8),allocatable,dimension(:)            :: conv_funct
  !custom variables for chempot search:
  character(len=32)                              :: ed_file_suffix
  logical                                        :: converged_n,upprshft
  integer                                        :: conv_n_loop=0
  integer                                        :: shift_n_loop=0
  real(8)                                        :: zJ1_2=0.d0,zJ3_2=0.d0
  real(8)                                        :: bottom,top
  real(8)                                        :: Alvl=0.d0
  real(8)                                        :: shift
  real(8)                                        :: dw,sumdens,xmu_old
  real(8),allocatable,dimension(:)               :: w
  real(8),allocatable,dimension(:)               :: orb_dens
  logical                                        :: look4n=.true.
  !custom variables for density matrix:
  real(8),allocatable,dimension(:)               :: dm_eig
  complex(8),allocatable,dimension(:,:)          :: dm,dm_rot
  complex(8),allocatable,dimension(:,:)          :: dm_custom_rot
  !custom variables for SOC expectations:
  complex(8),allocatable,dimension(:,:,:)        :: Stot
  complex(8),allocatable,dimension(:,:,:)        :: Ltot
  complex(8),allocatable,dimension(:)            :: jz
  !non interacting analysis:
  real(8)                                        :: mu
  logical                                        :: nonint_mu_shift!=.false.
  !
  !#########   MPI INITIALIZATION   #########
  !
#ifdef _MPI
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  !rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#else
  master=.true.
#endif
  !
  !#########    VARIABLE PARSING    #########
  !
  call parse_cmd_variable(finput,           "FINPUT",              default='inputED_SOC.in')
  call parse_input_variable(hkfile,         "HKFILE",finput,       default="hkfile.in")
  call parse_input_variable(nk,             "NK",finput,           default=10)
  call parse_input_variable(nkpath,         "NKPATH",finput,       default=500)
  call parse_input_variable(wmixing,        "WMIXING",finput,      default=0.5d0)
  call parse_input_variable(soc,            "SOC",finput,          default=0.0d0)
  call parse_input_variable(ivb,            "IVB",finput,          default=0.0d0)
  call parse_input_variable(surface,        "SURFACE",finput,      default=.false.)
  call parse_input_variable(Hk_test,        "HK_TEST",finput,      default=.true.)
  call parse_input_variable(upprshft,       "upprshft",finput,     default=.false.)
  call parse_input_variable(rotateG0loc,    "ROTATEG0loc",finput,  default=.false.)
  call parse_input_variable(nonint_mu_shift,"NONINTMUSHIFT",finput,default=.false.)
  !
#ifdef _MPI
  call ed_read_input(trim(finput),comm)
#else
  call ed_read_input(trim(finput))
#endif
  !
  Nso=Nspin*Norb
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(Lfit,"Lfit")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(ed_para,"ED_PARA")
  call add_ctrl_var(Jz_basis,"JZ_BASIS")
  call add_ctrl_var(ed_file_suffix,"ed_file_suffix")
  !
  !#########       ALLOCATION       #########
  !
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));            Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));            Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));            Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));            Greal=zero
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));            Weiss=zero
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats));            delta=zero
  allocate(Saux_wm(Lmats)); Saux_wm=zero
  allocate(Saux_wr(Lreal)); Saux_wr=zero
  !
  allocate(weiss_old(Nspin,Nspin,Norb,Norb,Lmats));        weiss_old=zero
  allocate(delta_old(Nspin,Nspin,Norb,Norb,Lmats));        delta_old=zero
  allocate(Sigma_correction(Nspin,Nspin,Norb,Norb,Lmats)); Sigma_correction=zero
  !
  allocate(conv_funct(Lmats));                             conv_funct=zero
  !
  allocate(dm(Nspin*Norb,Nspin*Norb));                     dm=zero
  allocate(dm_eig(Nspin*Norb));                            dm_eig=zero
  allocate(dm_rot(Nspin*Norb,Nspin*Norb));                 dm_rot=zero
  allocate(dm_custom_rot(Nspin*Norb,Nspin*Norb));          dm_custom_rot=zero
  !
  allocate(Stot(3,Norb,Norb));                             Stot=zero
  allocate(Ltot(3,Nspin,Nspin));                           Ltot=zero
  allocate(jz(3));                                         jz=zero
  !
  allocate(d_t2g_Hloc_so(Nspin*Norb,Nspin*Norb));          d_t2g_Hloc_so=zero
  allocate(d_t2g_Hloc_nn(Nspin,Nspin,Norb,Norb));          d_t2g_Hloc_nn=zero
  !
  allocate(w(Lreal));                                      w=0.0d0
  w = linspace(wini,wfin,Lreal,mesh=dw)
  !
  Lfit=min(int((Uloc(1)+1.2+(3*SOC/2))*(beta/pi))+100,Lmats)
  if(master)write(LOGfile,'(a12,I6,2(a12,F10.3))')"Lfit:",Lfit,"iwmax:",(pi/beta)*(2*Lfit-1),"U+2D+Dsoc:",Uloc(1)+1.2+(3*SOC/2)
  !
  !#########        BUILD Hk        #########
  !
  call build_hk(trim(hkfile),trim(hlocfile))               !;stop
  if(surface)then
      allocate(Wtk(Nk*Nk));Wtk=1.d0/(Nk*Nk)
  else
      allocate(Wtk(Nk*Nk*Nk));Wtk=1.d0/(Nk*Nk*Nk)
  endif
  !stop
  if(nonint_mu_shift)stop
  !
  !#########          BATH          #########
  !
  if (bath_type/="replica") then
     Nb=get_bath_dimension()
  else
     Nb=get_bath_dimension(d_t2g_Hloc_nn)
  endif
  if(master)write(LOGfile,*)"Bath_size:",Nb
  allocate(Bath(Nb));     Bath=0.0d0
  allocate(Bath_old(Nb)); Bath_old=0.0d0
  !
  !#########      INIT SOLVER       #########
  !
#ifdef _MPI
  call ed_init_solver(Comm,Bath,d_t2g_Hloc_nn)
#else
  call ed_init_solver(Bath,d_t2g_Hloc_nn)
#endif
  !
  !#########          DMFT          #########
  !
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !solve impurity
#ifdef _MPI
     call ed_solve(comm,Bath)
#else
     call ed_solve(Bath)
#endif
     !
     !get sigmas
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     if(ed_mode=="normal".and.SOC==0.0d0)then
        Saux_wm=zero
        Saux_wr=zero
        do ispin=1,Nspin
           do iorb=1,Norb
              Saux_wm=Saux_wm+Smats(ispin,ispin,iorb,iorb,:)/6.d0
              Saux_wr=Saux_wm+Sreal(ispin,ispin,iorb,iorb,:)/6.d0
           enddo
        enddo
        do ispin=1,Nspin
           do iorb=1,Norb
              Smats(ispin,ispin,iorb,iorb,:)=Saux_wm
              Sreal(ispin,ispin,iorb,iorb,:)=Saux_wr
           enddo
        enddo
     endif
     !
     !get local Gf's
#ifdef _MPI
     call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats);              call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=3)
     call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal);               call dmft_print_gf_realaxis(Greal,"Gloc",iprint=3)
#else
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats);                   call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=3)
     call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal) ;                   call dmft_print_gf_realaxis(Greal,"Gloc",iprint=3)
#endif
     !
     !operations on Weiss/Delta
     if(cg_scheme=='weiss')then
        !get Weiss
        call dmft_weiss(Gmats,Smats,Weiss,d_t2g_Hloc_nn) ;           call dmft_print_gf_matsubara(Weiss,"WeissG0",iprint=3)
        !mix Weiss
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_old
        !old Weiss
        Weiss_old=Weiss
        !fit Weiss
        if (ed_mode=="normal") then
#ifdef _MPI
           call ed_chi2_fitgf(Comm,Weiss,bath,ispin=1)
#else
           call ed_chi2_fitgf(Weiss,bath,ispin=1)
#endif
           call spin_symmetrize_bath(bath,save=.true.)
           !if(SOC==0.0d0)call  orb_symmetrize_bath(bath,save=.true.)
        else
#ifdef _MPI
           call ed_chi2_fitgf(Comm,Weiss,bath)
#else
           call ed_chi2_fitgf(Weiss,bath)
#endif
        endif
     else
        !get Delta
        call dmft_delta(Gmats,Smats,Delta,d_t2g_Hloc_nn);            call dmft_print_gf_matsubara(Delta,"Delta",iprint=3)
        !mix Delta
        if(iloop>1)Delta = wmixing*Delta + (1.d0-wmixing)*Delta_old
        !old Delta
        Delta_old=Delta
        !fit Delta
        if (ed_mode=="normal") then
#ifdef _MPI
           call ed_chi2_fitgf(Comm,Delta,bath,ispin=1)
#else
           call ed_chi2_fitgf(Delta,bath,ispin=1)
#endif
           call spin_symmetrize_bath(bath,save=.true.)
           !if(SOC==0.0d0)call  orb_symmetrize_bath(bath,save=.true.)
        else
#ifdef _MPI
           call ed_chi2_fitgf(Comm,Delta,bath)
#else
           call ed_chi2_fitgf(Delta,bath)
#endif
        endif
     endif
     !
     !each loop operations
     if(master)then
        !
        !
        !+ get rotation:
        !  imp_basis --> Jdiag_basis
        call build_rotation(dm_custom_rot)
        !
        !+ get rho in the impurity basis 
        !  and rotate it in the diagonal J basis
        call ed_get_density_matrix(dm,dm_custom_rot,dm_eig,dm_rot)
        !
        !+ print operators Simp, Limp, Jimp in the {a,s} basis
        call ed_get_quantum_SOC_operators_single()
        !
        !
        write(LOGfile,*)
        write(LOGfile,*) "   -------------------- rotations ---------------------"
        if(bath_type=="replica")then
           !
           ! rotation of impSmats in the {J} basis
           call Jz_rotate(Smats,"impS","wm")
           !
           ! rotation of impGmats in the {J} basis
           if(allocated(impG))deallocate(impG);allocate(impG(Nspin,Nspin,Norb,Norb,Lmats));impG=zero
           call ed_get_gimp_matsubara(impG)
           call Jz_rotate(impG,"impG","wm")
           deallocate(impG)
           !
           ! rotation of Greal in the {J} basis
           Alvl=0.2d0;bottom=0.d0;top=0.d0
           call Jz_rotate(Greal,"Gloc","wr",bottom,top,Alvl)
           !
        elseif(bath_type=="normal")then
           call compute_spectral_moments_nn(Greal,w,-1.d0/pi,dw)
        endif
        write(LOGfile,*) "   ----------------------------------------------------"
        write(LOGfile,*)
        !
        !e - chemical potential find
        converged_n=.true.
        xmu_old=xmu
        allocate(orb_dens(Norb));orb_dens=0.d0
        call ed_get_dens(orb_dens);sumdens=sum(orb_dens)
        deallocate(orb_dens)
        if(nread/=0.d0.and.look4n)then
           converged_n=.false.
           if(iloop>=2)call search_chempot(xmu,sumdens,converged_n)
        endif
        if(converged_n)then
           conv_n_loop=conv_n_loop+1
        else
           conv_n_loop=0
        endif
        !
        !f - convergence
        write(LOGfile,*)
        write(LOGfile,*) "   ------------------- convergence --------------------"
        if(cg_scheme=='weiss')then
           do i=1,Lmats
              conv_funct(i)=sum(nn2so_reshape(weiss(:,:,:,:,i),Nspin,Norb))
           enddo
        else
           do i=1,Lmats
              conv_funct(i)=sum(nn2so_reshape(delta(:,:,:,:,i),Nspin,Norb))
           enddo
        endif
        if(converged_n)converged = check_convergence(conv_funct,dmft_error,nsuccess,nloop)
        write(LOGfile,'(a35,L3)') "sigma converged",converged
        write(LOGfile,'(a35,L3)') "dens converged",converged_n
        converged = converged .and. converged_n
        write(LOGfile,'(a35,L3)') "total converged",converged
        write(LOGfile,'(a35,I3)') "global iloop",iloop
        write(LOGfile,'(a35,I3)') "times dens is ok",conv_n_loop
        write(LOGfile,'(a35,I3)') "times rigid shift",shift_n_loop
        write(LOGfile,*) "   ----------------------------------------------------"
        write(LOGfile,*)
        !
        !g - final mu shift
        !if(converged_n.and.upprshft.and.((nread==5.d0.and.zJ1_2<=0.01).or.(nread==2.d0.and.zJ3_2<=0.01)))then
        if(converged_n.and.upprshft) then !.and.(abs(nread-2.d0)<=nerr).or.(abs(nread-5.d0)<=nerr)))then
           write(LOGfile,*)
           write(LOGfile,*) "   -------------------- uppershift --------------------"
           !shift_n_loop=shift_n_loop+1
           !
           write(LOGfile,'(2(a10,F10.5))')"top:",top,"bottom:",bottom
           shift      = bottom + ( top - bottom ) / 2.d0
           xmu_old    = xmu
           if(abs(shift)>2*dw)then
              shift_n_loop=shift_n_loop+1
              xmu        = xmu_old + shift
              converged  = .false.
              look4n     = .false.
              write(LOGfile,'(6(a10,F10.5))')"shift:",shift,"xmu_old:",xmu_old,"xmu_new:",xmu
              unit=free_unit()
              open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
              write(unit,'(3F25.12,a10,1I5)')xmu,sumdens,shift,"shift",shift_n_loop
              close(unit)
           else
              write(LOGfile,'(6(a10,F10.5))')"NO shift:",shift,"2dw:",2*dw,"xmu_old:",xmu_old,"xmu_new:",xmu
           endif
           !
           !unit=free_unit()
           !open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
           !write(unit,*)xmu,sumdens,shift,"shift"
           !close(unit)
           !
        write(LOGfile,*) "   ----------------------------------------------------"
           write(LOGfile,*)
        endif
     endif
     !
#ifdef _MPI
     call Bcast_MPI(Comm,top)
     call Bcast_MPI(Comm,bottom)
     call Bcast_MPI(Comm,xmu)
     call Bcast_MPI(Comm,converged)
     call MPI_Barrier(Comm,ier)
#endif
     !
     if(master)call end_loop
     !
  enddo
  !
  !#########    BUILD Hk ON PATH    #########
  !
  call build_eigenbands()
  !
  !
#ifdef _MPI
  call finalize_MPI()
#endif
  !
  !
contains



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: build the Non interacting Hamiltonian with SOC and IVSB
  !+------------------------------------------------------------------------------------------+!
  subroutine build_hk(file1,file2)
    implicit none
    character(len=*),optional                    :: file1
    character(len=*),optional                    :: file2
    real(8),dimension(3)                         :: bk_x,bk_y,bk_z
    integer                                      :: ik,Lk
    !integer                                      :: i_mu,max_mu=500
    !real(8)                                      :: mu_edge=0.5
    integer                                      :: i_mu,max_mu=100
    real(8)                                      :: mu_edge=2.0d0
    complex(8),dimension(Nso,Nso,Lmats)          :: Gmats
    complex(8),dimension(Nso,Nso,Lreal)          :: Greal
    complex(8),allocatable                       :: Gso(:,:,:,:,:)
    real(8)                                      :: wm(Lmats),wr(Lreal),dw,mu
    !
    if(master)then
       write(LOGfile,*)"Build H(Nso,Nso,k)"
       write(LOGfile,*)"# of k-points per direction :",Nk
       write(LOGfile,*)"# of SO-bands               :",Nso
    endif
    if(allocated(Bath))stop" H(K) must be build before bath allocation, errors shall come otherwise"
    !
    bk_x = [1.d0,0.d0,0.d0]*2*pi
    bk_y = [0.d0,1.d0,0.d0]*2*pi
    bk_z = [0.d0,0.d0,1.d0]*2*pi
    !
    if(allocated(Hk))deallocate(Hk)
    !
    if(surface) then
       call TB_set_bk(bk_x,bk_y)
       Lk=Nk*Nk
       if(master)write(LOGfile,*)"surface tot k-points:",Lk
       allocate(Hk(Nso,Nso,Lk));Hk=zero
       Sigma_correction=zero
       call TB_build_model(Hk,hk_Ti3dt2g,Nso,[Nk,Nk])
       if(master.AND.present(file1)) call TB_write_hk(Hk,file1,Nso,Norb,1,1,[Nk,Nk])
    else
       call TB_set_bk(bk_x,bk_y,bk_z)
       Lk=Nk*Nk*Nk
       if(master)write(LOGfile,*)"bulk tot k-points:",Lk
       allocate(Hk(Nso,Nso,Lk));Hk=zero
       Sigma_correction=zero
       call TB_build_model(Hk,hk_Ti3dt2g,Nso,[Nk,Nk,Nk])
       if(master.AND.present(file1)) call TB_write_hk(Hk,file1,Nso,Norb,1,1,[Nk,Nk,Nk])
    endif
    !
    d_t2g_Hloc_so = sum(Hk(:,:,:),dim=3)/Lk
    where(abs((d_t2g_Hloc_so))<1.d-9)d_t2g_Hloc_so=0d0
    d_t2g_Hloc_nn=so2nn_reshape(d_t2g_Hloc_so,Nspin,Norb)
    call TB_write_hloc(d_t2g_Hloc_so,file2)
    !
    !-----  Build the local GF in the spin-orbital Basis   -----
    !
    !matsu freq
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    Gmats=zero
    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i) , Hk(:,:,ik) )/Lk
       enddo
    enddo
    !
    allocate(Gso(Nspin,Nspin,Norb,Norb,Lmats));Gso=zero
    do i=1,Lmats
       Gso(:,:,:,:,i)=so2nn_reshape(Gmats(:,:,i),Nspin,Norb)
    enddo
    !
    if(master)      call dmft_print_gf_matsubara(Gso,"G0loc",iprint=3)
    if(rotateG0loc) call Jz_rotate(Gso,"G0lc","wm")
    deallocate(Gso)
    !
    !real freq
    if(.not.nonint_mu_shift)then
      max_mu=1
      mu_edge=0.d0
    endif
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    !
    do i_mu=1,max_mu
       if(nonint_mu_shift)then
          !mu=-mu_edge+(2*abs(mu_edge)/max_mu)*float(i_mu)
          mu=0.25d0+(2*abs(mu_edge)/max_mu)*float(i_mu)
       else
          if(SOC==0.3d0)then
             mu=0.3585 !(L=0.3)
          elseif(SOC==0.1d0)then
             mu=0.2819 !(L=0.1)
          else
             mu=xmu
          endif
       endif
       Greal=zero
       do ik=1,Lk
          do i=1,Lreal
             Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps),Hk(:,:,ik),mu)/Lk
          enddo
       enddo
       !
       allocate(Gso(Nspin,Nspin,Norb,Norb,Lreal));Gso=zero
       do i=1,Lreal
          Gso(:,:,:,:,i)=so2nn_reshape(Greal(:,:,i),Nspin,Norb)
       enddo
       !
       if(master)      call dmft_print_gf_realaxis(Gso,"G0loc",iprint=3)
       if(rotateG0loc) call Jz_rotate(Gso,"G0lc","wr")
       deallocate(Gso)
       !
    enddo
    !
  end subroutine build_hk



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: function that produces the full non interacting Hamiltonian
  !+------------------------------------------------------------------------------------------+!
  function hk_Ti3dt2g(kvec,Nso_) result(hk)
    real(8),dimension(:)                         :: kvec
    real(8)                                      :: kx,ky,kz
    integer                                      :: Nso_,ndx
    complex(8),dimension(Nso_,Nso_)              :: hk
    real(8),allocatable                          :: HoppingMatrix(:,:)
    complex(8),allocatable                       :: U(:,:),Udag(:,:)
    !
    if(surface)then
       kx=kvec(1);ky=kvec(2)
    else
       kx=kvec(1);ky=kvec(2);kz=kvec(3)
    endif
    !
    allocate(HoppingMatrix(Norb,0:6));HoppingMatrix=0.0d0
    call get_hopping(HoppingMatrix)
    !
    Hk=zero
    do i=1,Norb
       ndx=2*i-1
       Hk(ndx:ndx+1,ndx:ndx+1) = diagonal_orbital_dispersion(kx,ky,kz,HoppingMatrix(i,:))
    enddo
    !
    if((SOC/=zero.or.IVB/=zero).and.bath_type=="replica")then
       Hk = Hk + SOC*H_LS() + IVB*H_IVB(kx,ky,kz)
    endif
    !
    !correction with Sigma(iw=0)
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = ispin + (iorb-1)*Nspin
                jo = jspin + (jorb-1)*Nspin
                Hk(io,jo) = Hk(io,jo) + real(Sigma_correction(ispin,jspin,iorb,jorb,1))
             enddo
          enddo
       enddo
    enddo
    !
    !shape:  [Norb*Norb]*Nspin
    !
    !basis:  {a,Sz}
    Hk = so2os_reshape(Hk,Nspin,Norb)
    !basis:  {Lz,Sz}
    if(Jz_basis)then
       allocate(U(Nspin*Norb,Nspin*Norb));U=zero
       allocate(Udag(Nspin*Norb,Nspin*Norb));Udag=zero
       U=orbital_Lz_rotation_NorbNspin()
       Udag=transpose(conjg(orbital_Lz_rotation_NorbNspin()))
       Hk=matmul(Udag,matmul(Hk,U))
    endif
    !
  end function hk_Ti3dt2g



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: (DIAGONAL) build local SOC contribution in the Z formulation
  !+------------------------------------------------------------------------------------------+!
  function diagonal_orbital_dispersion(kx,ky,kz,t) result(hk)
    real(8),intent(in)                           :: kx,ky,kz
    real(8),intent(in),dimension(0:6)            :: t
    complex(8),dimension(2,2)                    :: hk
    real(8)                                      :: t0
    !
    if(Hk_test)then
       t0=0.1
    else
       t0=1.d0
    endif
    !
    !perovskite dispersion
    hk = zero
    if (surface) then
       if(Hk_test)then
          !surface model dispersion cosine on x, y
          !up
          hk(1,1) = t(0)+(                       & !onsite_orbX
                -2.*t(1)*cos(kx)                 & !t_100_orbX
                -2.*t(2)*cos(ky)                 & !t_010_orbX
          !      -1.*t(3))*t0                       !t_001_orbX
                )*t0                       !t_001_orbX

          !dw
          hk(2,2) = hk(1,1)
       else
          !surface realistic dispersion on x, y
          !up
          hk(1,1) = t(0)+(                       & !onsite_orbX
                -2.*t(1)*cos(kx)                 & !t_100_orbX
                -2.*t(2)*cos(ky)                 & !t_010_orbX
                -1.*t(3)                         & !t_001_orbX
                -2.*t(4)*cos(ky)                 & !t_011_orbX
                -2.*t(5)*cos(kx)                 & !t_101_orbX
                -4.*t(6)*cos(kx)*cos(ky))*t0       !t_110_orbX
          !dw
          hk(2,2) = hk(1,1)
       endif
    else
       if(Hk_test)then
          !bulk model dispersion cosine on x, y, z
          !up
          hk(1,1) = t(0)+(                       & !onsite_orbX
                -2.*t(1)*cos(kx)                 & !t_100_orbX
                -2.*t(2)*cos(ky)                 & !t_010_orbX
                -2.*t(3)*cos(kz))*t0               !t_001_orbX
          !dw
          hk(2,2) = hk(1,1)
       else
          !bulk realistic dispersion on x, y, z
          !up
          hk(1,1) = t(0)+(                       & !onsite_orbX
                -2.*t(1)*cos(kx)                 & !t_100_orbX
                -2.*t(2)*cos(ky)                 & !t_010_orbX
                -2.*t(3)*cos(kz)                 & !t_001_orbX
                -4.*t(4)*cos(ky)*cos(kz)         & !t_011_orbX
                -4.*t(5)*cos(kx)*cos(kz)         & !t_101_orbX
                -4.*t(6)*cos(kx)*cos(ky))*t0       !t_110_orbX
          !dw
          hk(2,2) = hk(1,1)
       endif
    endif
  end function diagonal_orbital_dispersion



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: (OFF DIAGONAL) build local SOC contribution in the Z formulation
  !+------------------------------------------------------------------------------------------+!
  function H_LS() result(hls)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: hls
    hls=zero
    hls(1:2,3:4) = +Xi * pauli_z / 2.
    hls(1:2,5:6) = -Xi * pauli_y / 2.
    hls(3:4,5:6) = +Xi * pauli_x / 2.
    !hermiticity
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          hls(j,i)=conjg(hls(i,j))
       enddo
    enddo
  end function H_LS



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: (OFF DIAGONAL) build local IVSB band mixing contribution in the Z formulation
  !+------------------------------------------------------------------------------------------+!
  function H_IVB(kx,ky,kz) result(hivb)
    real(8),intent(in)                           :: kx,ky,kz
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: hivb
    hivb=zero
    hivb(1:2,3:4) = zero
    hivb(1:2,5:6) = 2*xi*sin(kx)*eye(2) 
    hivb(3:4,5:6) = 2*xi*sin(ky)*eye(2)
    !hermiticity
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          hivb(j,i)=conjg(hivb(i,j))
       enddo
    enddo
  end function H_IVB



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: solve H(k) along path in the BZ.
  !+------------------------------------------------------------------------------------------+!
  subroutine build_eigenbands()
    real(8),dimension(:,:),allocatable           :: kpath
    integer                                      :: Npts,Lk
    type(rgb_color),dimension(:),allocatable     :: colors
    !
    if (master)then
       allocate(colors(Nso))
       colors=[red1,green1,blue1,red1,green1,blue1]
       !
       if(surface)then
          write(LOGfile,*)"Build surface H(k) along the path M-X-G-X"
          Npts = 4
          Lk=(Npts-1)*Nkpath
          allocate(kpath(Npts,3))
          kpath(1,:)=kpoint_M1
          kpath(2,:)=kpoint_X1
          kpath(3,:)=kpoint_Gamma
          kpath(4,:)=kpoint_X1
          !
          Sigma_correction=zero
          call TB_solve_model(hk_Ti3dt2g,Nso,kpath,Nkpath,colors_name=colors,points_name=[character(len=20) :: 'M', 'X', 'G', 'X'],file="Eigenband_surf.nint")
          !call compute_Akw_along_BZpath(hk_Ti3dt2g,Nso,kpath,Nkpath,Sreal,points_name=[character(len=20) :: 'M', 'X', 'G', 'X'])
          !
          Sigma_correction=Smats
          call TB_solve_model(hk_Ti3dt2g,Nso,kpath,Lk,colors_name=colors,points_name=[character(len=20) :: 'M', 'X', 'G', 'X'],file="Eigenband_surf.sigma")
          !
          Sigma_correction=zero
          !
       else
          write(LOGfile,*)"Build bulk H(k) along the path M-R-G-M-X-G-X"
          Npts = 7
          Lk=(Npts-1)*Nkpath
          allocate(kpath(Npts,3))
          kpath(1,:)=kpoint_M1
          kpath(2,:)=kpoint_R
          kpath(3,:)=kpoint_Gamma
          kpath(4,:)=kpoint_M1
          kpath(5,:)=kpoint_X1
          kpath(6,:)=kpoint_Gamma
          kpath(7,:)=kpoint_X1
          !
          Sigma_correction=zero
          call TB_solve_model(hk_Ti3dt2g,Nso,kpath,Nkpath,colors_name=colors,points_name=[character(len=20) :: 'M', 'R', 'G', 'M', 'X', 'G', 'X'],file="Eigenband_bulk.nint")
          call compute_Akw_along_BZpath(hk_Ti3dt2g,kpath,Nkpath,Sreal)
          !
          Sigma_correction=Smats
          call TB_solve_model(hk_Ti3dt2g,Nso,kpath,Nkpath,colors_name=colors,points_name=[character(len=20) :: 'M', 'R', 'G', 'M', 'X', 'G', 'X'],file="Eigenband_bulk.sigma")
          !
          Sigma_correction=zero
          !
       endif
       !
       write(LOGfile,*)"Im done on the path"
       stop
    endif
    !
  end subroutine build_eigenbands



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Compute A(k,w) along path in the BZ.
  !+------------------------------------------------------------------------------------------+!
  subroutine compute_Akw_along_BZpath(hk_model,kpath,Nk,Self)
    implicit none
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)                    :: kpoint
         integer                                 :: N
         complex(8),dimension(N,N)               :: hk_model
       end function hk_model
    end interface
    real(8),dimension(:,:),intent(in)            :: kpath
    integer               ,intent(in)            :: Nk
    complex(8)            ,intent(in)            :: Self(Nspin,Nspin,Norb,Norb,Lreal)
    !internal
    character(len=3)                             :: file_
    integer                                      :: Npts,Lfreq
    real(8),dimension(size(kpath,2))             :: kstart,kstop,kpoint,kdiff
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: h,U,Udag
    complex(8)                                   :: Self_basis(Nspin*Norb,Nspin*Norb,Lreal)
    integer                                      :: ipts,ik,ic,unit,u1,u2,iorb,ifreq
    real(8)                                      :: wr(Lreal),dw,mu
    real(8),allocatable,dimension(:,:,:,:)       :: Akw,Akw_rot
    !
    file_="Akw"
    Npts=size(kpath,1)
    Lfreq=size(Self,dim=5)
    if(allocated(Akw))    deallocate(Akw)    ;allocate(Akw(Nspin*Norb,Nspin*Norb,(Npts-1)*Nk,Lfreq))    ;Akw=zero
    if(allocated(Akw_rot))deallocate(Akw_rot);allocate(Akw_rot(Nspin*Norb,Nspin*Norb,(Npts-1)*Nk,Lfreq));Akw_rot=zero
    !
    U=zero;Udag=zero
    if(ed_para)then
       if(Jz_basis)then
          U=matmul(transpose(conjg(orbital_Lz_rotation_NorbNspin())),atomic_SOC_rotation())
       else
          U=atomic_SOC_rotation()
       endif
    else
       U=eye(Nspin*Norb)
    endif
    Udag=transpose(conjg(U))
    !
    do ifreq=1,Lfreq
       Self_basis(:,:,ifreq)=nn2so_reshape(Self(:,:,:,:,ifreq),Nspin,Norb)
    enddo
    !
    wr = linspace(wini,wfin,Lfreq,mesh=dw)
    ic = 0
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/dble(Nk)
       do ik=1,Nk
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hk_model(kpoint,Nspin*Norb)
          do ifreq=1,Lfreq
             Akw(:,:,ic,ifreq) = -aimag(inverse_g0k(dcmplx(wr(ifreq),eps),h+Self_basis(:,:,ifreq),xmu))/pi
             Akw_rot(:,:,ic,ifreq) = matmul(Udag,matmul(Akw(:,:,ic,ifreq),U))
          enddo
       enddo
    enddo
    !Tr[Akw] and Tr[Akw_rot] give the same result the trace is invariant !! Wht changes is the relative weight of the orbital components
    open(unit=106,file='Akw.dat',status='unknown',action='write',position='rewind')
    do ifreq=1,Lfreq
       write(106,'(9000F18.12)')wr(ifreq),(trace(Akw(:,:,ic,ifreq)*2.)+5.*ic,ic=1,(Npts-1)*Nk)
    enddo
    close(106)
    !test
    open(unit=106,file='Akwbis.dat',status='unknown',action='write',position='rewind')
    do ifreq=1,Lfreq
       write(106,'(9000F18.12)')wr(ifreq),(trace(Akw_rot(:,:,ic,ifreq)*2.)+5.*ic,ic=1,(Npts-1)*Nk)
    enddo
    close(106)
    !band 1
    open(unit=106,file='Akw_b1.dat',status='unknown',action='write',position='rewind')
    do ifreq=1,Lfreq
       write(106,'(9000F18.12)')wr(ifreq),((Akw_rot(1,1,ic,ifreq)*12.)+5.*ic,ic=1,(Npts-1)*Nk)
    enddo
    close(106)
    !band 2
    open(unit=106,file='Akw_b2.dat',status='unknown',action='write',position='rewind')
    do ifreq=1,Lfreq
       write(106,'(9000F18.12)')wr(ifreq),((Akw_rot(3,3,ic,ifreq)*12.)+5.*ic,ic=1,(Npts-1)*Nk)
    enddo
    close(106)

  end subroutine compute_Akw_along_BZpath



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE:   Build the hopping integrals matrix for realistic bandstructure
  !STRUCTURE: T[orbital(yz,zx,xy) nn direction(100,010,001),nnn direction(011,101,110)]
  !+------------------------------------------------------------------------------------------+!
  subroutine get_hopping(T)
    real(8),dimension(Norb,0:6),intent(out)      ::  T
    real(8),dimension(3,0:6)                     ::  T_bulk,T_LAOSTO
    real(8)                                      ::  Eo,t1,t2,t3
    real(8)                                      ::  t_010_yz,t_001_yz
    real(8)                                      ::  t_100_zx,t_001_zx
    real(8)                                      ::  t_100_xy,t_010_xy,t_001_xy
    !
    T=0.d0
    !
    !pristine lattice
    Eo = 3.31
    t1 = 0.276536
    t2 = 0.031329
    t3 = 0.076842
    !
    !lattice distortion
    t_010_yz = 0.232 !se c'è solo l'abbassamento del Ti questo dovrebbe essere uguale a t1, magari c'è anche altro dovuto ad LAO
    t_001_yz = 0.475
    !
    t_100_zx = 0.232
    t_001_zx = 0.475
    !
    t_100_xy = 0.286
    t_010_xy = 0.286
    t_001_xy = 0.03
    !
    !####  BULK STO  ####
    !orbital_1 = YZ
    T_bulk(1,0) = Eo
    T_bulk(1,1) = t2
    T_bulk(1,2) = t1
    T_bulk(1,3) = t1
    T_bulk(1,4) = t3
    T_bulk(1,5) = 0.d0
    T_bulk(1,6) = 0.d0
    !orbital_2 = ZX
    T_bulk(2,0) = Eo
    T_bulk(2,1) = t1
    T_bulk(2,2) = t2
    T_bulk(2,3) = t1
    T_bulk(2,4) = 0.d0
    T_bulk(2,5) = t3
    T_bulk(2,6) = 0.d0
    !orbital_3 = XY
    T_bulk(3,0) = Eo
    T_bulk(3,1) = t1
    T_bulk(3,2) = t1
    T_bulk(3,3) = t2
    T_bulk(3,4) = 0.d0
    T_bulk(3,5) = 0.d0
    T_bulk(3,6) = t3
    !
    !####  LAO/STO  ####
    !orbital_1 = YZ
    T_LAOSTO(1,0) = 1.087
    T_LAOSTO(1,1) = t2
    T_LAOSTO(1,2) = t_010_yz
    T_LAOSTO(1,3) = t_001_yz
    T_LAOSTO(1,4) = t3
    T_LAOSTO(1,5) = 0.d0
    T_LAOSTO(1,6) = 0.d0
    !orbital_2 = ZX
    T_LAOSTO(2,0) = 1.087
    T_LAOSTO(2,1) = t_100_zx
    T_LAOSTO(2,2) = t2
    T_LAOSTO(2,3) = t_001_zx
    T_LAOSTO(2,4) = 0.d0
    T_LAOSTO(2,5) = t3
    T_LAOSTO(2,6) = 0.d0
    !orbital_3 = XY
    T_LAOSTO(3,0) = 1.035
    T_LAOSTO(3,1) = t_100_xy
    T_LAOSTO(3,2) = t_010_xy
    T_LAOSTO(3,3) = t_001_xy
    T_LAOSTO(3,4) = 0.d0
    T_LAOSTO(3,5) = 0.d0
    T_LAOSTO(3,6) = t3
    !
    !
    if(Hk_test)then
       if(bath_type=="replica")then
          T=1.0d0
          T(1,0) = 0.0d0
          T(2,0) = 0.0d0
          T(3,0) = 0.0d0
       else
          T=1.0d0
          if(surface)then
             T(1,0) = +SOC
             T(2,0) = 0.0d0
             T(3,0) = 0.0d0
          else
             T(1,0) = +SOC
             T(2,0) = -SOC/2.d0
             T(3,0) = -SOC/2.d0
          endif
       endif
    else
       if(surface)then
          T=T_LAOSTO(1:Norb,:)
       else
          T=T_bulk(1:Norb,:)
       endif
    endif
    !
    !
  end subroutine get_hopping


  !____________________________________________________________________________________________!
  !                         Operators & Operations related to SOC
  !____________________________________________________________________________________________!
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Build the rotations
  !+------------------------------------------------------------------------------------------+!
  subroutine build_rotation(theta_C_,impHloc_rot_)
    implicit none
    complex(8),dimension(6,6),intent(out)          ::   theta_C_
    complex(8),dimension(6,6),intent(out),optional ::   impHloc_rot_
    real(8),dimension(6)                           ::   impHloc_eig
    !
    theta_C_ = zero
    if(Jz_basis)then
       !rotation {Lz,Sz}-->{a,Sz}-->{J}
       theta_C_ = matmul(transpose(conjg(orbital_Lz_rotation_NorbNspin())),atomic_SOC_rotation())
    else
       !rotation {a,Sz}-->{J}
       theta_C_ = atomic_SOC_rotation()
    endif
    !
    if(present(impHloc_rot_))then
       impHloc_rot_=zero
       impHloc_rot_=d_t2g_Hloc_so
       call eigh(impHloc_rot_,impHloc_eig,'V','U')
    endif
    !
  end subroutine build_rotation



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: rotations on DMFT functions G0loc(w/iw), Gloc(w), impG(iw), Sigma(iw)
  !         with theta_C or theta_rho
  !NOTE:    Compulsory input variables are:
  !         type_funct = G0lc/Gloc/impG/impS
  !         type_freq  = wm/wr
  !+------------------------------------------------------------------------------------------+!
  subroutine Jz_rotate(Fso,type_funct,type_freq,bottom_,top_,lvl_)
    implicit none
    !input variables
    complex(8),allocatable,intent(in)                ::   Fso(:,:,:,:,:)
    character(len=4)      ,intent(in)                ::   type_funct
    character(len=2)      ,intent(in)                ::   type_freq
    real(8)               ,intent(out),optional      ::   bottom_,top_
    real(8)               ,intent(in) ,optional      ::   lvl_
    !aux rotated fnct in so/nso
    complex(8),allocatable                           ::   f_in(:,:,:),f_out(:,:,:),Gimp(:,:,:,:,:)
    !output fnct in nn/nnn
    complex(8),allocatable                           ::   Fso_out(:,:,:,:,:)
    !array & indexes
    real(8)   ,allocatable,dimension(:)              ::   w
    integer                                          ::   io,jo,ndx,Lfreq,ik
    integer                                          ::   ispin,jspin,iorb,jorb
    !flags & names
    integer                                          ::   isetup=0
    logical                                          ::   level
    character(len=10)                                ::   file_rotation
    character(len=91)                                ::   header
    !rotations & dm
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      ::   theta_C
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      ::   rho_ab
    real(8)   ,dimension(Nspin*Norb,Nspin*Norb)      ::   dens_rot
    !top bottom
    integer                                          ::   posupper,poslower
    integer                                          ::   icount,max_count
    integer   ,dimension(200)                        ::   posmax
    real(8)                                          ::   second_derivative
    !observables
    real(8)                                          ::   norm,fact,lvl,dw
    real(8)   ,dimension(Nspin*Norb)                 ::   z_rot
    complex(8),dimension(Nspin*Norb)                 ::   Luttinger
    real(8)                                          ::   LS_0,jz_0,jz_0_sq!,Ek0
    !non interacting Ek
!    real(8)                                      ::   Ek0bis
!    complex(8),allocatable                       ::   Hkj(:,:,:)
    !
    if(type_funct=="G0lc") isetup=1
    if(type_funct=="Gloc") isetup=2
    if(type_funct=="impS") isetup=3
    if(type_funct=="impG") isetup=4
    lvl=1.0d0;if(present(lvl_))lvl=lvl_
    !
    !if(Jz_basis)       theta_C: {Lz,Sz}-->{a,Sz}-->{J}
    !if(.not.Jz_basis)  theta_C: {a,Sz}------------>{J}
    call build_rotation(theta_C)
    header="#    {1/2,-1/2}     {1/2,+1/2}     {3/2,-3/2}     {3/2,+3/2}     {3/2,-1/2}     {3/2,+1/2}"
    !
    !function allocation
    Lfreq=size(Fso,dim=5)
    if(allocated( f_in))  deallocate( f_in);  allocate(   f_in(Nspin*Norb,Nspin*Norb,Lfreq));f_in=zero
    if(allocated(f_out))  deallocate(f_out);  allocate(  f_out(Nspin*Norb,Nspin*Norb,Lfreq));f_out=zero
    if(allocated(Fso_out))deallocate(Fso_out);allocate(Fso_out(Nspin,Nspin,Norb,Norb,Lfreq));Fso_out=zero
    !
    !observables allocation
    if(isetup==1) then
       write(LOGfile,*) "  G0loc rotation"
    elseif(isetup==2) then
       write(LOGfile,*) "  Gloc rotation"
    elseif(isetup==3) then
       write(LOGfile,*) "  impS rotation"
    elseif(isetup==4) then
       write(LOGfile,*) "  impG rotation"
       if(allocated(Gimp))deallocate(Gimp);allocate(Gimp(Nspin,Nspin,Norb,Norb,Lfreq));Gimp=zero
       call ed_get_gimp_matsubara(Gimp)
    endif
    !
    !meshes
    if(type_freq=="wr")then
       if(allocated(w))deallocate(w)
       allocate(w(Lreal));w=0.d0
       w = linspace(wini,wfin,Lreal,mesh=dw)
       norm=dw
       fact=-1.d0/pi
       write(LOGfile,'(A11,2F9.4)') "   real freq",norm,fact
    elseif(type_freq=="wm")then
       if(allocated(w))deallocate(w)
       allocate(w(Lmats));w=0.d0
       w = pi/beta*(2*arange(1,Lmats)-1)
       norm=1.d0/beta
       fact=1.d0
       write(LOGfile,'(A11,2F9.4)') "   imag freq",norm,fact
    endif
    !
    !function intake
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                f_in(io,jo,:)=Fso(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    !
    !save the integral before rotation
    if(isetup==1) then
       open(unit=106,file='sum_'//type_freq//'_G0loc.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==2) then
       open(unit=106,file='sum_'//type_freq//'_Gloc.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==3) then
       open(unit=106,file='sum_'//type_freq//'_impS.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==4) then
       open(unit=106,file='sum_'//type_freq//'_impG.dat',status='unknown',action='write',position='rewind')
    endif
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,"(I3,I3,A3,4I3,F18.12)")io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(f_in(io,jo,:)))*norm
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    !
    !
    !###############################################################
    !#                                                             #
    !#                           ROTATIONS                         #
    !#                                                             #
    !###############################################################
    !
    write(LOGfile,*) "  Rotation with analytic LS"
    !
    !1)rotation
    f_out=zero
    do i=1,Lfreq
       f_out(:,:,i)=matmul(transpose(conjg(theta_C)),matmul(f_in(:,:,i),theta_C))
    enddo
    !
    !2)Zqp save in the case f_in = Smats(iw)
    if(isetup==3 .and. type_freq=="wm")then
       z_rot=0d0
       do io=1,Nspin*Norb
          z_rot(io)   = 1.d0/( 1.d0 + abs( dimag(f_out(io,io,1))/(pi/beta) ))
       enddo
       zJ1_2=z_rot(1)
       zJ3_2=z_rot(4)
       open(unit=106,file='Zqp_rot.dat',status='unknown',action='write',position='rewind')
       write(106,'(A)') header
       write(106,'(90F15.9,1X)')(z_rot(io),io=1,Nspin*Norb)
       close(106)
    endif
    !
    !3)Luttinger save in the case f_in = impG(iw)
    if(isetup==4 .and. type_freq=="wm")then
       Luttinger=zero
       do io=1,Nspin*Norb
          Luttinger(Nspin*Norb)=f_out(io,io,1)  
       enddo
       open(unit=106,file='Luttinger.dat',status='unknown',action='write',position='rewind')
       write(106,'(A)') header
       write(106,'(6F15.9,1X,6F15.9)') (real(Luttinger(io)),io=1,Nspin*Norb),(aimag(luttinger(io)),io=1,Nspin*Norb)
       close(106)
    endif
    !
    !4)top-bottom find of the half-filled band in the case f_in = Gloc(w) and N=2,5
    if(isetup==2 .and. type_freq=="wr" )then ! .and. upprshft )then
       if( (abs(nread-2.d0)<=2*nerr).or.(abs(nread-5.d0)<=2*nerr) )then
          if(abs(nread-5.d0)<=2*nerr)ndx=1
          if(abs(nread-2.d0)<=2*nerr)ndx=3
          if(present(top_).and.present(bottom_))then
             top_=0.d0;bottom_=0.d0
             posupper=10*Lfreq;poslower=-posupper
             max_count=0;posmax=0
             !
             freqloop:do i=Lfreq,1,-1
                if(( abs(real(f_out(ndx,ndx,i))).lt.lvl  ).and.( real(f_out(ndx,ndx,i))>0.d0)) then
                   max_count=max_count+1
                   posmax(max_count)=i
                   if(w(i)<-wfin) exit freqloop
                endif
             enddo freqloop
             !
             maxloop:do icount=1,max_count
                level=.false.
                if( -aimag(f_out(ndx,ndx,posmax(icount)))/pi>0.85 )level=.true.
                second_derivative = (real(f_out(ndx,ndx,posmax(icount)+1))-real(f_out(ndx,ndx,posmax(icount)-1)))/(w(posmax(icount)+1)-w(posmax(icount)-1))
                if(second_derivative>0.d0)then
                   if((posmax(icount)<posupper).and.(w(posmax(icount))>0.d0).and.level)then
                      posupper=posmax(icount)
                      top_=w(posupper)
                   elseif((w(posmax(icount))<0.d0).and.level)then
                      poslower=posmax(icount)
                      bottom_=w(poslower)
                      exit maxloop
                   endif
                endif
             enddo maxloop
             !
          endif
       endif
    endif
    !
    !5)first moment of A(w) in the case f_in = Gloc(w) and N=2,5
    if(isetup==2 .and. type_freq=="wr" )then
       call compute_spectral_moments_so(f_out,w,fact,norm)
    endif
    !
    !6)save the integral after rotation
    if(isetup==1) then
       open(unit=106,file='sum_'//type_freq//'_G0loc_rot.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==2) then
       open(unit=106,file='sum_'//type_freq//'_Gloc_rot.dat' ,status='unknown',action='write',position='rewind')
    elseif(isetup==3) then
       open(unit=106,file='sum_'//type_freq//'_impS_rot.dat' ,status='unknown',action='write',position='rewind')
    elseif(isetup==4) then
       open(unit=106,file='sum_'//type_freq//'_impG_rot.dat' ,status='unknown',action='write',position='rewind')
    endif
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,"(I3,I3,A3,4I3,F18.12)")io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(f_out(io,jo,:)))*norm
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    !7)save the rotated function
    if(isetup==1) then
       file_rotation="G0lc_rot"
     elseif(isetup==2) then
       file_rotation="Gloc_rot"
    elseif(isetup==3) then
       file_rotation="impS_rot"
    elseif(isetup==4) then
       file_rotation="impG_rot"
    endif
    Fso_out=zero
    do i=1,Lfreq
       Fso_out(:,:,:,:,i)=so2nn_reshape(f_out(:,:,i),Nspin,Norb)
    enddo
    if(type_freq=="wr") call dmft_print_gf_realaxis( Fso_out,file_rotation,iprint=3)
    if(type_freq=="wm") call dmft_print_gf_matsubara(Fso_out,file_rotation,iprint=3)
    !
    !8)non interacting rho and observables in the case f_in = G0loc(w)
    if(isetup==1.and.type_freq=="wr")then
       dens_rot=0.d0
       do io=1,Nspin*Norb
          do i=1,Lreal
             dens_rot(io,io)=dens_rot(io,io)+fact*dimag(f_out(io,io,i))*norm
             if(abs(w(i))<dw) exit
          enddo
       enddo
       !
       rho_ab = matmul(theta_C,matmul(dens_rot,transpose(conjg(theta_C))))
       LS_0=trace(matmul(rho_ab,atomic_SOC()))
       jz_0=trace(matmul(rho_ab,atomic_j("z")))
       jz_0_sq=trace(matmul(rho_ab,matmul(atomic_j("z"),atomic_j("z"))))
       !
!               Ek0=zero
!               do ik=1,Nk*Nk*Nk
!                  Ek0 = Ek0 + trace(matmul(rho_ab,Hk(:,:,ik)))/(Nk*Nk*Nk)
!               enddo
!               !
!               if(allocated(HkJ))deallocate(HkJ);allocate(HkJ(Nspin*Norb,Nspin*Norb,Nk*Nk*Nk));HkJ=zero
!               Ek0bis=zero
!               do ik=1,Nk*Nk*Nk
!                  HkJ(:,:,ik)=matmul(transpose(conjg(theta_C)),matmul(Hk(:,:,ik),theta_C))-mu*eye(Nspin*Norb)
!                  do io=1,Nspin*Norb
!                     if(real(HkJ(io,io,ik))>0.d0)cycle
!                     Ek0bis = Ek0bis + HkJ(io,io,ik)/(Nk*Nk*Nk)
!                  enddo
!               enddo
       !
       write(LOGfile,'(1A12,A72,1A12)') "mu","non-interacting J basis densities","Ntot"
       write(LOGfile,'(20F12.4)') mu,(dens_rot(io,io),io=1,Nso),sum(dens_rot)
       write(LOGfile,'(20A12)')   " LS_0","jz_0","jz_0_sq"!,"Ek0"
       write(LOGfile,'(20F12.4)') LS_0,jz_0,jz_0_sq!,Ek0,Ek0bis
       !
       open(unit=106,file='nonint_dens_rot.dat',status='unknown',action='write',position='append')
       write(106,'(20F18.12)')   mu,(dens_rot(io,io),io=1,Nso),LS_0,jz_0,jz_0_sq!,Ek0
       close(106)
       !
       open(106,file='nonint_density_matrix.dat',action='write',position='append',status='unknown')
       write(106,*)
       write(106,"(A10,F22.12)")"# mu:",mu
       write(106,*)
       write(106,"(A10)")"# Re{rho_nonint}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(106,"(90(F15.9,1X))")  (real(rho_ab(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(106,"(A10)")
         write(106,"(A10)")"# Im{rho_nonint}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(106,"(90(F15.9,1X))") (aimag(rho_ab(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(106,"(A10)")
       write(106,"(4A22)")"#mu","rho_tilda","trace{rho_tilda}","LS_nonint"
       write(106,'(10F22.12)')  mu,(dens_rot(io,io),io=1,Nso),sum(dens_rot),LS_0
       write(106,"(A10)")
       close(106)
       !
    endif
    !
  end subroutine Jz_rotate


  !+------------------------------------------------------------------------------------------+!
  !PURPOSE:
  !NOTE:    
  !+------------------------------------------------------------------------------------------+!
  subroutine compute_spectral_moments_so(A_,w_,fact_,norm_)
    implicit none
    complex(8),intent(in)                         :: A_(Nspin*Norb,Nspin*Norb,Lreal)
    real(8),intent(in)                            :: w_(Lreal)
    real(8),intent(in)                            :: fact_,norm_
    real(8)                                       :: moment1_J12=0.d0,moment2_J12=0.d0
    real(8)                                       :: moment1_J32=0.d0,moment2_J32=0.d0
    integer                                       :: i
    !
    moment1_J12=0.d0;moment2_J12=0.d0
    moment1_J32=0.d0;moment2_J32=0.d0
    do i=1,Lreal
       moment1_J12=moment1_J12+fact_*aimag(A_(1,1,i))*w_(i)*norm_
       moment1_J32=moment1_J32+fact_*aimag(A_(3,3,i))*w_(i)*norm_
       moment2_J12=moment2_J12+fact_*aimag(A_(1,1,i))*w_(i)*w_(i)*norm_
       moment2_J32=moment2_J32+fact_*aimag(A_(3,3,i))*w_(i)*w_(i)*norm_
    enddo
    if(master)then
       open(unit=106,file='Spectral_moment.dat',status='unknown',action='write',position='rewind')
       write(106,'(90A15,1X)')"#A_J1/2(w)*w","A_J3/2(w)*w","A_J1/2(w)*w^2","A_J3/2(w)*w^2"
       write(106,'(90F15.9,1X)')moment1_J12,moment1_J32,moment2_J12,moment2_J32
       close(106)
    endif
    !
  end subroutine compute_spectral_moments_so
  !
  subroutine compute_spectral_moments_nn(A_,w_,fact_,norm_)
    implicit none
    complex(8),intent(in)                         :: A_(Nspin,Nspin,Norb,Norb,Lreal)
    real(8),intent(in)                            :: w_(Lreal)
    real(8),intent(in)                            :: fact_,norm_
    real(8)                                       :: moment1_J12=0.d0,moment2_J12=0.d0
    real(8)                                       :: moment1_J32=0.d0,moment2_J32=0.d0
    integer                                       :: i
    !
    moment1_J12=0.d0;moment2_J12=0.d0
    moment1_J32=0.d0;moment2_J32=0.d0
    do i=1,Lreal
       moment1_J12=moment1_J12+fact_*aimag(A_(1,1,1,1,i))*w_(i)*norm_
       moment1_J32=moment1_J32+fact_*aimag(A_(1,1,2,2,i))*w_(i)*norm_
       moment2_J12=moment2_J12+fact_*aimag(A_(1,1,1,1,i))*w_(i)*w_(i)*norm_
       moment2_J32=moment2_J32+fact_*aimag(A_(1,1,2,2,i))*w_(i)*w_(i)*norm_
    enddo
    if(master)then
       open(unit=106,file='Spectral_moment.dat',status='unknown',action='write',position='rewind')
       write(106,'(90A15,1X)')"#A_J1/2(w)*w","A_J3/2(w)*w","A_J1/2(w)*w^2","A_J3/2(w)*w^2"
       write(106,'(90F15.9,1X)')moment1_J12,moment1_J32,moment2_J12,moment2_J32
       close(106)
    endif
    !
  end subroutine compute_spectral_moments_nn



  !____________________________________________________________________________________________!
  !                                       Gfs
  !____________________________________________________________________________________________!
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: G0_loc functions
  !+------------------------------------------------------------------------------------------+!
  function inverse_g0k(iw,hk,mu_) result(g0k)
    implicit none
    complex(8)                                    :: iw
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: hk
    real(8),intent(in),optional                   :: mu_
    !real(8)                                       :: mu
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: g0k,g0k_tmp
    integer                                       :: i,ndx
    integer (kind=4), dimension(6)                :: ipiv
    integer (kind=1)                              :: ok
    integer (kind=4), parameter                   :: lwork=2000
    complex (kind=8), dimension(lwork)            :: work
    real    (kind=8), dimension(lwork)            :: rwork
    !
    mu=0.d0
    if(present(mu_))mu=mu_
    !
    g0k=zero;g0k_tmp=zero
    !
    g0k=(iw+mu)*eye(Nspin*Norb)-hk
    g0k_tmp=g0k
    !
    call inv(g0k)
    call inversion_test(g0k,g0k_tmp,1.e-9)
  end function inverse_g0k



  !____________________________________________________________________________________________!
  !                                     utilities
  !____________________________________________________________________________________________!
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Inversion test
  !+------------------------------------------------------------------------------------------+!
  subroutine inversion_test(A,B,tol)
    implicit none
    complex (kind=8), intent(in)   ::   A(Nspin*Norb,Nspin*Norb)
    complex (kind=8), intent(in)   ::   B(Nspin*Norb,Nspin*Norb)
    real    (kind=4), intent(in)   ::   tol
    integer (kind=2)               ::   dime

    if (size(A).ne.size(B)) then
       write(LOGfile,*) "Matrices not equal cannot perform inversion test"
       stop
    endif
    dime=maxval(shape(A))
    if (abs(float(dime)-real(sum(matmul(A,B)))).gt.tol) write(LOGfile,'(A30)') "inversion test fail"
  end subroutine inversion_test


end program ed_SOC
