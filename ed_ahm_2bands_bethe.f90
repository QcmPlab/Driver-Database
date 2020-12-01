program ed_ahm_2bands_bethe
  USE DMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                       :: iloop,Nb,Le,Nso
  logical                                       :: converged
  real(8)                                       :: alpha,wmixing,Eout(2),dens
  !Bath:
  real(8),allocatable                           :: Bath(:),BathOld(:)
  !
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:,:,:)
  real(8),dimension(2)                          :: Wband,de
  real(8),dimension(:,:),allocatable            :: Dbands
  real(8),dimension(:,:),allocatable            :: Ebands ![Nso][Le]
  real(8),dimension(:),allocatable              :: H0     ![Nso]
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal,Smats,Sreal,Delta
  character(len=16)                             :: finput
  logical                                       :: phsym,normal_bath


  call parse_cmd_variable(finput,"FINPUT",default='inputAHM.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(alpha,"ALPHA",finput,default=1d0,comment="bandwidth ratio W_2 = alpha*W_1=alpha*1.0")
  call parse_input_variable(Le,"LE",finput,default=500)
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




  !Allocate Weiss Field:
  allocate(delta(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(2,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(2,Nspin,Nspin,Norb,Norb,Lreal))


  !Build Hamiltonian structure:
  Wband(1)=1d0
  Wband(2)=alpha*Wband(1)
  !
  allocate(Ebands(Nso,Le))
  Ebands(1,:) = linspace(-Wband(1),Wband(1),Le,mesh=de(1))
  Ebands(2,:) = linspace(-Wband(2),Wband(2),Le,mesh=de(2))
  !
  allocate(Dbands(Nso,Le))    
  Dbands(1,:) = dens_bethe(Ebands(1,:),Wband(1))*de(1)
  Dbands(2,:) = dens_bethe(Ebands(2,:),Wband(2))*de(2)
  !
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc   = zero
  !
  allocate(H0(Nso))
  H0=zero


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

     !Retrieve impurity self-energies (normal, anomalous)
     call ed_get_Smats(Smats)

     !Compute the local gfs:
     call dmft_gloc_matsubara(Ebands,Dbands,H0,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats(1,:,:,:,:,:),"Gloc",iprint=1)
     call dmft_print_gf_matsubara(Gmats(2,:,:,:,:,:),"Floc",iprint=1)


     call dmft_self_consistency(&
          Gmats(1,:,:,:,:,:),Gmats(2,:,:,:,:,:),&
          Smats(1,:,:,:,:,:),Smats(2,:,:,:,:,:),&
          Delta,Hloc,trim(cg_scheme))


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
  call ed_get_Sreal(Sreal)
  call dmft_gloc_realaxis(Ebands,Dbands,H0,Greal,Sreal)

  call dmft_print_gf_matsubara(Gmats(1,:,:,:,:,:),"Gloc",iprint=1)
  call dmft_print_gf_matsubara(Gmats(2,:,:,:,:,:),"Floc",iprint=1)
  call dmft_print_gf_realaxis(Greal(1,:,:,:,:,:),"Gloc",iprint=1)
  call dmft_print_gf_realaxis(Greal(2,:,:,:,:,:),"Floc",iprint=1)

  !Compute the Kinetic Energy:
  call dmft_kinetic_energy(Ebands,Dbands,H0,Smats(1,:,:,:,:,:),Smats(2,:,:,:,:,:))



end program ed_ahm_2bands_bethe



