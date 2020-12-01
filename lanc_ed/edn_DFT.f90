program ed_W90
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                                     :: Nlat,Nktot,Nkpath,Nkvec(3),Npts,Nlso
  integer                                     :: iloop
  integer                                     :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  integer                                     :: unit  
  real(8)                                     :: wmixing
  real(8),dimension(3)                        :: e1,e2,e3
  real(8),dimension(:,:),allocatable          :: kpath
  real(8)                                     :: ef,filling
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:),allocatable       :: Hloc
  real(8),dimension(:),allocatable            :: Wtk,dens
  character(len=60)                           :: w90file,InputFile,latfile,kpathfile,hkfile
  character(len=40),allocatable               :: points_name(:)
  logical                                     :: converged
  logical                                     :: EFflag
  !Bath:
  integer                                     :: Nb
  real(8)   ,allocatable,dimension(:)         :: Bath
  !local dmft Weisss:
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Weiss_prev
  !Mpi:
  integer                                     :: comm,rank,ier
  logical                                     :: master=.true.,bool


  !
#ifdef _MPI
  call init_MPI
  comm = MPI_COMM_WORLD
  master = get_master_MPI()
  rank = get_Rank_MPI()
#endif
  !
  !
  !#########    VARIABLE PARSING    #########
  !
  call parse_cmd_variable(InputFile,"INPUTFILE",default="inputDFT.conf")
  call parse_input_variable(NLAT,"NLAT",InputFile,         default=1)
  call parse_input_variable(Nkvec,"NKVEC",InputFile,default=[10,10,10])
  call parse_input_variable(nkpath,"NKPATH",InputFile,default=500)
  call parse_input_variable(wmixing,"WMIXING",InputFile,      default=0.5d0)
  call parse_input_variable(w90file,"w90file",InputFile,default="hij.conf")
  call parse_input_variable(hkfile,"hkfile",InputFile,default="hk.conf")
  call parse_input_variable(latfile,"latfile",InputFile,default="lat.conf")
  call parse_input_variable(kpathfile,"kpathfile",InputFile,default="kpath.conf")
  call parse_input_variable(EFflag,"EFflag",InputFile,default=.false.)
  call parse_input_variable(filling,"filling",InputFile,default=1d0)
  call ed_read_input(reg(InputFile),comm)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  !>DEVEL
  if(Nlat/=1)stop "This driver is now only for Nlat=1"
  !<DEVEL
  Nlso = Nlat*Nspin*Norb

  !METHOD 1 (setup W90 --> use internal W90 model)
  inquire(file=reg(latfile),exist=bool)
  if(bool)then
     open(free_unit(unit),file=reg(latfile))
     read(unit,*)e1
     read(unit,*)e2
     read(unit,*)e3
     close(unit)
  else
     e1 = [1d0,0d0,0d0]
     e2 = [0d0,1d0,0d0]
     e3 = [0d0,0d0,1d0]
  endif
  call TB_set_ei(e1,e2,e3)
  call TB_build_bk(verbose=.true.)
  call start_timer
  call TB_w90_setup(reg(w90file),nlat=Nlat,nspin=Nspin,norb=Norb,verbose=.true.)
  call stop_timer("TB_w90_setup")
  !
  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  inquire(file=reg(hkfile),exist=bool)
  if(bool)then
     call TB_read_hk(Hk,reg(hkfile),Nkvec)
     if(size(Hk,1)/=Nlat*Nspin*Norb)stop "ed_DFT error: wrong size in Hk as read from file"
     Nktot = size(Hk,3)
  else
     if(EFflag)then
        write(*,"(A)")"get Fermi Level"
        call start_timer  
        call TB_w90_FermiLevel(Nkvec,filling,Ef)
        call stop_timer("TB_w90_FermiLevel")
     endif
     !
     Nktot=product(Nkvec)
     allocate(Hk(Nlso,Nlso,Nktot))
     call TB_set_dos_lreal(256)
     call start_timer
     call TB_build_model(Hk,Nlso,Nkvec,wdos=.true.)
     call TB_write_hk(reg(hkfile),Nkvec)
     call stop_timer("TB_build_model")
  endif

  allocate(Wtk(Nktot))
  Wtk=1d0/Nktot

  allocate(Hloc(Nlso,Nlso))
  Hloc= sum(Hk(:,:,:),dim=3)/size(Hk,3)
  where(abs(Hloc)<1d-6)Hloc=zero
  if(master)call TB_write_Hloc(Hloc,"w90Hloc.dat")

  write(*,*)"Using Nk_total="//str(size(Hk,3))


  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));      Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));      Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));      Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));      Greal=zero
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));      Weiss=zero
  allocate(Weiss_prev(Nspin,Nspin,Norb,Norb,Lmats)); Weiss_prev=zero
  allocate(dens(Norb))

  Nb=get_bath_dimension(so2nn(Hloc))
  allocate(Bath(Nb));    Bath=0.0d0

  call ed_init_solver(comm,Bath,so2nn(Hloc))
  !
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !--------  solve impurity (diagonal basis)  --------
     call ed_solve(comm,Bath,so2nn(Hloc))
     !
     !--------    get sigmas   (diagonal basis)  --------
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)

     !------  get local Gf's  (Wannier90 basis)  ------
     call dmft_gloc_matsubara(comm,Hk,Wtk,Gmats,Smats,mpi_split='k')
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     !------    get Weiss     (Wannier90 basis)  --------
     call dmft_self_consistency(comm,Gmats,Smats,Weiss,so2nn(Hloc),cg_scheme)
     if(master)call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)


     !------    mix Weiss     ------
     if(iloop>1)then
        Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_prev
     endif
     Weiss_prev=Weiss
     !
     !------    fit Weiss     ------
     call ed_chi2_fitgf(comm,Weiss,Bath,ispin=1)
     !
     if(master)then
        if(nread/=0d0)then
           call ed_get_dens(dens)
           call search_chemical_potential(xmu,sum(dens),converged)
        endif
        !
        !Check convergence (if required change chemical potential)
        converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)
     endif
     !
     call Bcast_MPI(comm,xmu)
     call Bcast_MPI(comm,converged)
     !
     if(master)call end_loop
     !
  enddo


  !Compute the local gfs:
  call dmft_gloc_realaxis(comm,Hk,Wtk,Greal,Sreal,mpi_split='k')
  if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)

  !Compute the Kinetic Energy:
  call dmft_kinetic_energy(comm,Hk,Wtk,Smats)


  call finalize_MPI()

contains


  function so2nn(Hlso) result(Hnnn)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hlso
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                     :: iorb,jorb
    integer                                     :: ispin,jspin
    integer                                     :: is,js
    Hnnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Hnnn(ispin,jspin,iorb,jorb) = Hlso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function so2nn


  function nn2so(Hnnn) result(Hlso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnnn
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hlso
    integer                                     :: iorb,jorb
    integer                                     :: ispin,jspin
    integer                                     :: is,js
    Hlso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Hlso(is,js) = Hnnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function nn2so

end program ed_W90
