program ed_ahm_bethe
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                       :: iloop
  integer                                       :: Nb
  logical                                       :: converged
  !Bath:
  real(8),allocatable,dimension(:)              :: Bath,Bath_Prev
  !The local functions:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  !
  character(len=16)                             :: finput
  integer                                       :: Le
  real(8)                                       :: wband,wmixing
  logical                                       :: phsym,normal_bath
  real(8)                                       :: de,dens
  real(8),dimension(:,:),allocatable            :: Ebethe,Dbethe
  real(8),dimension(:),allocatable              :: H0
  complex(8),dimension(:,:,:,:),allocatable     :: Hloc
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputAHM.conf')
  call parse_input_variable(wband,"wband",finput,default=1.d0,comment="Bethe Lattice bandwidth")
  call parse_input_variable(Le,"Le",finput,default=500,comment="Number of energy levels for Bethe DOS integration")
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(phsym,"phsym",finput,default=.false.,comment="Flag to enforce p-h symmetry of the bath.")
  call parse_input_variable(normal_bath,"normal",finput,default=.false.,comment="Flag to enforce no symmetry braking in the bath.")
  !
  call ed_read_input(trim(finput))
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  !
  !
  !Allocate local functions (Nambu)
  allocate(Gmats(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss(2,Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(Greal(2,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Sreal(2,Nspin,Nspin,Norb,Norb,Lreal))
  !
  !
  !
  allocate(Ebethe(1,Le))  
  Ebethe(1,:) = linspace(-Wband,Wband,Le,mesh=de)
  !
  allocate(Dbethe(1,Le))
  Dbethe(1,:) = dens_bethe(Ebethe(1,:),wband)*de
  !
  allocate(Hloc(1,1,1,1))
  Hloc        = 0d0
  !
  allocate(H0(1))
  H0=0d0
  !

  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_prev(Nb))
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
     call dmft_gloc_matsubara(Ebethe,Dbethe,H0,Gmats,Smats)

     call dmft_self_consistency(&
          Gmats(1,:,:,:,:,:),Gmats(2,:,:,:,:,:),&
          Smats(1,:,:,:,:,:),Smats(2,:,:,:,:,:),&
          Weiss,Hloc,trim(cg_scheme))

     !Perform the self-consistency fitting the new bath
     call ed_chi2_fitgf(Weiss,bath,ispin=1)
     !if it holds apply symmetrizations 
     if(phsym)call ph_symmetrize_bath(bath,save=.true.)
     if(normal_bath)call enforce_normal_bath(bath,save=.true.)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     !Check convergence (if required change chemical potential)
     converged = check_convergence(Weiss(1,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)

     if(nread/=0.d0)then
        call ed_get_dens(dens,iorb=1)
        call search_chemical_potential(xmu,dens,converged)
     end if

     !Close this DMFT loop
     call end_loop
  enddo

  !Compute the local gfs:
  call ed_get_Sreal(Sreal)
  call dmft_gloc_realaxis(Ebethe,Dbethe,H0,Greal,Sreal)

  call dmft_print_gf_matsubara(Gmats(1,:,:,:,:,:),"Gloc",iprint=1)
  call dmft_print_gf_matsubara(Gmats(2,:,:,:,:,:),"Floc",iprint=1)
  call dmft_print_gf_realaxis(Greal(1,:,:,:,:,:),"Gloc",iprint=1)
  call dmft_print_gf_realaxis(Greal(2,:,:,:,:,:),"Floc",iprint=1)

  !Compute the Kinetic Energy:
  call dmft_kinetic_energy(Ebethe,Dbethe,H0,Smats(1,:,:,:,:,:),Smats(2,:,:,:,:,:))



end program ed_ahm_bethe



