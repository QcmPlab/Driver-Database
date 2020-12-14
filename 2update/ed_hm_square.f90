program ed_hm_square
  USE DMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                       :: iloop,Nb,Lk,Nx,Nso,ik
  logical                                       :: converged
  real(8)                                       :: wband,ts,wmixing
  !Bath:
  real(8),allocatable                           :: Bath(:),BathOld(:)
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gkmats
  !
  character(len=16)                             :: finput
  complex(8),allocatable                        :: Hk(:,:,:)
  real(8),allocatable                           :: Wt(:)

  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
  call parse_input_variable(Nx,"Nx",finput,default=10,comment="Number of kx point for 2d BZ integration")
  !
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  Nso=1

  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))



  !Build Hk
  call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
  Lk = Nx*Nx
  allocate(Hk(Nso,Nso,Lk),Wt(Lk),Hloc(Nspin,Nspin,Norb,Norb))
  call TB_build_model(Hk(:,:,:),hk_model,Nso,[Nx,Nx])
  Wt = 1d0/Lk
  Hloc   = zero
  call TB_write_hk(Hk(:,:,:),"Hk2d.dat",1,&
       Nd=1,Np=0,Nineq=1,&
       Nkvec=[Nx,Nx])

  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bathold(Nb))
  call ed_init_solver(bath,Hloc)



  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath,Hloc) 
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)

     !Compute the local gfs:
     call dmft_gloc_matsubara(Hk,Wt,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc)
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc)
     endif

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(Weiss,bath)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*BathOld
     BathOld=Bath

     !Check convergence (if required change chemical potential)
     converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)

     call end_loop
  enddo

  !Compute the local gfs:
  call dmft_gloc_realaxis(Hk,Wt,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)


  !Get kinetic energy:
  call dmft_kinetic_energy(Hk,Wt,Smats)


  allocate(Gkmats(Lk,Nspin,Nspin,Norb,Norb,Lmats))
  do ik=1,Lk
     call dmft_gk_matsubara(Hk(:,:,ik),Wt(ik),Gkmats(ik,:,:,:,:,:),Smats)
  enddo



contains




  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky
    complex(8)           :: hk(N,N)
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = -one*2d0*ts*(cos(kx)+cos(ky))
  end function hk_model





end program ed_hm_square



