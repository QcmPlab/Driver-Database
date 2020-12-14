program ed_ahm_2bands
  USE DMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                       :: iloop,Nb,Lk,Nx,Nso
  logical                                       :: converged
  real(8)                                       :: wband,ts,alpha,wmixing,Eout(2),dens
  !Bath:
  real(8),allocatable                           :: Bath(:),BathOld(:)
  !
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal,Smats,Sreal,Delta
  character(len=16)                             :: finput,fhloc
  logical                                       :: phsym,normal_bath
  real(8),allocatable                           :: wt(:)
  complex(8),allocatable                        :: Hk(:,:,:,:)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.5d0,comment="hopping parameter")
  call parse_input_variable(alpha,"ALPHA",finput,default=1d0,comment="bandwidth ratio t_2 = alpha*t_1")
  call parse_input_variable(Nx,"Nx",finput,default=10,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(phsym,"phsym",finput,default=.false.,comment="Flag to enforce p-h symmetry of the bath.")
  call parse_input_variable(normal_bath,"normal",finput,default=.false.,comment="Flag to enforce no symmetry braking in the bath.")
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


  Nso=Nspin*Norb
  if(Nspin/=1.OR.Norb/=2)stop "This code is intended as a driver for the Norb=2 and Nspin=1 problem"


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bathold(Nb))
  call ed_init_solver(bath)


  !Allocate Weiss Field:
  allocate(delta(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(2,Nspin,Nspin,Norb,Norb,Lmats),Greal(2,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(2,Nspin,Nspin,Norb,Norb,Lmats),Sreal(2,Nspin,Nspin,Norb,Norb,Lreal))




  !Build Hk
  call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
  Lk = Nx*Nx
  allocate(Hk(2,Nso,Nso,Lk))
  allocate(Wt(Lk))
  call TB_build_model(Hk(1,:,:,:),hk_model,Nso,[Nx,Nx])
  Wt = 1d0/Lk
  Hk(2,:,:,:) = -Hk(1,:,:,:)
  !
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc   = zero
  call TB_write_hk(Hk(1,:,:,:),"Hk2d.dat",Nso,&
       Nd=Norb,Np=0,Nineq=Norb,&
       Nkvec=[Nx,Nx])


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath) 
     call ed_get_sigma_matsubara(Smats(1,:,:,:,:,:))
     call ed_get_self_matsubara(Smats(2,:,:,:,:,:))
     call ed_get_sigma_real(Sreal(1,:,:,:,:,:))
     call ed_get_self_real(Sreal(2,:,:,:,:,:))


     !Compute the local gfs:
     call dmft_gloc_matsubara_superc(Hk,Wt,Gmats,Smats,iprint=1)

     if(cg_scheme=='weiss')then
        call dmft_weiss_superc(Gmats,Smats,Delta,Hloc,iprint=1)
     else
        call dmft_delta_superc(Gmats,Smats,Delta,Hloc,iprint=1)
     endif

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(delta,bath,ispin=1)
     if(phsym)call ph_symmetrize_bath(bath)
     if(normal_bath)call enforce_normal_bath(bath)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*BathOld
     BathOld=Bath

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0)then
        call ed_get_dens(dens,iorb=1)
        call search_chemical_potential(xmu,dens,converged)
     endif
     call end_loop
  enddo

  !Compute the local gfs:
  call dmft_gloc_realaxis_superc(Hk,Wt,Greal,Sreal,iprint=1)

  !Compute the Kinetic Energy:
  call dmft_kinetic_energy(Hk(1,:,:,:),Wt,Smats(1,:,:,:,:,:),Smats(2,:,:,:,:,:))



contains




  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky
    complex(8)           :: hk(N,N)
    if(N/=2)stop "hk_model error: N!=2"
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = zero
    Hk(1,1) = -one*2d0*ts*(cos(kx)+cos(ky))
    Hk(2,2) = -one*2d0*alpha*ts*(cos(kx)+cos(ky))
  end function hk_model





end program ed_ahm_2bands



