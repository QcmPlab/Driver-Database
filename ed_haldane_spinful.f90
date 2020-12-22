program ed_haldane
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso,Nlat,Nineq
  logical                                       :: converged
  integer                                       :: ilat

  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: kmHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk

  real(8),dimension(2)                          :: e1,e2   !real-space lattice basis
  real(8),dimension(2)                          :: bk1,bk2 !reciprocal space lattice basis
  real(8),dimension(2)                          :: d1,d2,d3
  real(8),dimension(2)                          :: bklen

  !variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: t1,t2,phi,Mh,wmixing
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym,timerevsym,afmkick
  !
  integer                                       :: comm,rank
  logical                                       :: master

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputHALDANE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in",comment='Hk will be written here')
  call parse_input_variable(nk,"NK",finput,default=100,comment='Number of kpoints per directoin')
  call parse_input_variable(nkpath,"NKPATH",finput,default=500,comment='Number of kpoints per intervall on kpath')
  call parse_input_variable(t1,"T1",finput,default=1d0,comment='Nearest neighbour hopping')
  call parse_input_variable(t2,"T2",finput,default=0.2d0,comment='Next-nearest neighbour hopping amplitude')
  call parse_input_variable(phi,"PHI",finput,default=pi/2d0,comment='Next-nearest neighbour hpping phase. With TIMEREVSYM=T and PHI=pi/2 this corrensponds to the Kane-Mele Model without Rashba terms.')
  call parse_input_variable(mh,"MH",finput,default=0d0,comment='On-site staggering i.e. Semenoff-Mass')
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.6d0,comment='Mixing parameter: 0 means 100% of the old bath (no update at all), 1 means 100% of the new bath (pure update). Use a large value with care, decrease a lot near transitions.')
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.,comment='SU(2) constraint. Incompatible with magnetic order. Incompatible with SO interaction. Be careful.')
  call parse_input_variable(afmkick,"AFMKICK",finput,default=.false.)
  call parse_input_variable(timerevsym,"TIMEREVSYM",finput,default=.false.,comment='Use time-reversal symmetric version of Hamiltonian')
  !
  call ed_read_input(trim(finput),comm)
  !
  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Norb/=1.OR.Nspin/=2)stop "Wrong setup from input file: Norb!=1 OR Nspin!=2"
  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso                 !=4 = 2(ineq sites)*2(spin)*1(orb)

  !Lattice basis (a=1; a0=sqrt3*a) is:
  !e_1 = a0 [ sqrt3/2 , 1/2 ] = 3/2a[1, 1/sqrt3]
  !e_2 = a0 [ sqrt3/2 ,-1/2 ] = 3/2a[1,-1/sqrt3]
  e1 = 3d0/2d0*[1d0, 1d0/sqrt(3d0)]
  e2 = 3d0/2d0*[1d0,-1d0/sqrt(3d0)]
  !call TB_set_ei(eix=e1,eiy=e2)

  !LATTICE BASIS: nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]

  !RECIPROCAL LATTICE VECTORS:
  bklen=2d0*pi/3d0
  bk1=bklen*[ 1d0, sqrt(3d0)]
  bk2=bklen*[ 1d0,-sqrt(3d0)]
  call TB_set_bk(bkx=bk1,bky=bk2)

  !Build the Hamiltonian on a grid or on path
  call build_hk(trim(hkfile))
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  Hloc = lso2nnn_reshape(kmHloc,Nlat,Nspin,Norb)

  !Allocate Weiss Field:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero

  !Setup solver
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver(comm,Bath,Hloc)
  if(afmkick.AND.spinsym)stop "Wrong setup from input file: AFMKICK=T not with SPINSYM=T"
  if((.not.spinsym).and.afmkick)then
     if(master)write(*,*) "**********************************************"
     if(master)write(*,*) "*  Applying AFM kick to initial Bath. This   *"
     if(master)write(*,*) "*  calculation will tend to an AFM solution! *"
     if(master)write(*,*) "**********************************************"
     do ilat=1,Nlat
        call ed_break_symmetry_bath(Bath(ilat,:),sb_field,(-1d0)**(ilat+1))
     enddo
  endif

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,*) "Rank ", rank, " reporting loop number: ", iloop
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     !Solve separately the two atoms:
     do ilat=1,Nlat
        write(*,*) "Rank ", rank, " reporting solving atom", ilat
        call ed_set_suffix(ilat) !this is need to print different files for different sites
        call ed_solve(comm,Bath(ilat,:),Hloc(ilat,:,:,:,:))
        write(*,*) "Rank ", rank, " reporting solved atom", ilat
        !get the important functions to be used later: next call they will be overwritten.
        write(*,*) "Rank ", rank, " reporting get sig_mat atom", ilat
        call ed_get_sigma_matsubara(Smats(ilat,:,:,:,:,:))
        write(*,*) "Rank ", rank, " reporting get sig_real atom", ilat
        call ed_get_sigma_realaxis(Sreal(ilat,:,:,:,:,:))
     enddo
     call ed_reset_suffix()
     !
     ! compute the local gf:
     write(*,*) "Rank ", rank, " reporting gloc_mats computation"
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     write(*,*) "Rank ", rank, " reporting gloc_mats computation finished"
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     !
     ! compute the Weiss field (only the Nineq ones)
     write(*,*) "Rank ", rank, " reporting weiss computation"
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
     write(*,*) "Rank ", rank, " reporting dmft_weiss"

     !Fit the new bath, starting from the old bath + the supplied delta
     !Behaves differently depending on the ed_mode input:
     !IF(NORMAL): normal/magZ phase is solved, so either fit spin1 or spin1&2
     !IF(NONSU2): Sz-conservation broken, SOC or magXY, fit all components
     !IF(SUPERC): error, superconductivity is not allowed here
     select case(ed_mode)
     case default
        stop "ed_mode!=Normal/Nonsu2"
     case("normal")
        call ed_chi2_fitgf(comm,Bath,Weiss,Hloc,ispin=1)
        if(.not.spinsym)then
           call ed_chi2_fitgf(comm,Bath,Weiss,Hloc,ispin=2)
        else
           call ed_spin_symmetrize_bath(bath,save=.true.)
        endif
     case("nonsu2")
        call ed_chi2_fitgf(comm,Bath,Weiss,Hloc)
     end select

     !MIXING:
     write(*,*) "Rank ", rank, " calculating new bath"
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     !Check convergence. This is now entirely MPI-aware:
     converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  write(*,*) "Rank ", rank, " exited DMFT loop"
  if(master)write(*,*) "Done with DMFT loops"
  if(master)write(*,*) "print G_matsubara"
  call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=4)

  !Extract and print retarded self-energy and Green's function
  if(master)write(*,*) "extracting real axis local Green's function"
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  if(master)write(*,*) "printing real axis local Green's function"
  call dmft_print_gf_realaxis(Greal,"Greal",iprint=4)

  if(master)write(*,*) "getting kinetic energy"
  call dmft_kinetic_energy(Hk,Smats)

  if(master) then
     write(*,*) "!***************************!"
     write(*,*) "!*                         *!"
     write(*,*) "!*   !!!  FINISHED  !!!    *!"
     write(*,*) "!*                         *!"
     write(*,*) "!***************************!"
  endif

  call finalize_MPI()

contains



  !---------------------------------------------------------------------
  !PURPOSE: Get Haldane Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional          :: file
    real(8),dimension(2)               :: pointK,pointKp
    real(8),dimension(:,:),allocatable :: KPath
    real(8),dimension(:,:),allocatable :: kgrid
    real(8),dimension(:),allocatable   :: gridx,gridy
    integer                            :: i,j,ik
    !
    Lk= Nk*Nk
    if(master)write(*,*)"Build H(k) Haldane:",Lk
    if(master)write(*,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    !
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0
    !

    if(.not.timerevsym)then
        call TB_build_model(Hk,hk_haldane_model,Nlso,[Nk,Nk])
    else
        call TB_build_model(Hk,hk_haldane_timerevsym_model,Nlso,[Nk,Nk])
    endif
    if(master)call TB_write_Hk(Hk,file,Nlat,Nspin,Norb,[Nk,Nk])
    !
    allocate(kmHloc(Nlso,Nlso))
    kmHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(kmHloc))<1.d-4)kmHloc=0d0
    if(master)call TB_write_Hloc(kmHloc,'Hloc.txt')
    !
    pointK  = 1d0/3d0*bk1 - 1d0/3d0*bk2
    pointKp = 2d0/3d0*bk1 + 1d0/3d0*bk2
    write(*,*) pointK
    write(*,*) pointKp
    allocate(Kpath(4,2))
    KPath(1,:)=[0,0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointKp
    KPath(4,:)=[0d0,0d0]
    if(.not.timerevsym)then
        call TB_Solve_model(hk_haldane_model,Nlso,KPath,Nkpath,&
             colors_name=[red1,blue1,red1,blue1],&
             points_name=[character(len=10) :: "G","K","K`","G"],&
             file="Eigenbands.nint",iproject=.false.)
    else
        call TB_Solve_model(hk_haldane_timerevsym_model,Nlso,KPath,Nkpath,&
             colors_name=[red1,blue1,red1,blue1],&
             points_name=[character(len=10) :: "G","K","K`","G"],&
             file="Eigenbands.nint",iproject=.false.)
    endif

  end subroutine build_hk


  !--------------------------------------------------------------------!
  !Haldane HAMILTONIAN:                                                !
  !--------------------------------------------------------------------!

  function hk_haldane_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(2,2)       :: hk11,hk22
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz, kdote1, kdote2
    !
    !----- consistent with Bernvig ---------
    ! This way of creating the Hamiltonian we get a Bloch structure and can obtain local
    ! (i.e. in the unitcel) quantities via averaging over all k-points
    kdote1 = dot_product(kpoint,e1)
    kdote2 = dot_product(kpoint,e2)
    !
    h0 = 2*t2*cos(phi)*( cos(kdote1) + cos(kdote2) + cos(kdote1-kdote2) )
    hx = t1*( cos(kdote1) + cos(kdote2) + 1)
    hy = t1*( sin(kdote1) + sin(kdote2) )
    hz = 2*t2*sin(phi)*( sin(kdote1) - sin(kdote2) - sin(kdote1-kdote2) )
    !
    hk11 = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Mh*pauli_z
    hk22 = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Mh*pauli_z
    !
    !create hk such that it works well with lso2nnn_reshape(kmHloc,Nlat,Nspin,Norb)
    hk          = zero
    hk(1:3:2,1:3:2) = hk11
    hk(2:4:2,2:4:2) = hk22
    !
  end function hk_haldane_model

  function hk_haldane_timerevsym_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(2,2)       :: hk11,hk22
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz, kdote1, kdote2
    !
    !----- consistent with Bernvig ---------
    ! This way of creating the Hamiltonian we get a Bloch structure and can obtain local
    ! (i.e. in the unitcel) quantities via averaging over all k-points. This version is
    ! time-reversal symmetric due to the extra sign infront of hz
    kdote1 = dot_product(kpoint,e1)
    kdote2 = dot_product(kpoint,e2)
    !
    h0 = 2*t2*cos(phi)*( cos(kdote1) + cos(kdote2) + cos(kdote1-kdote2) )
    hx = t1*( cos(kdote1) + cos(kdote2) + 1)
    hy = t1*( sin(kdote1) + sin(kdote2) )
    hz = 2*t2*sin(phi)*( sin(kdote1) - sin(kdote2) - sin(kdote1-kdote2) )
    !
    hk11 = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Mh*pauli_z
    hk22 = h0*pauli_0 + hx*pauli_x + hy*pauli_y - hz*pauli_z + Mh*pauli_z
    !
    !create hk such that it works well with lso2nnn_reshape(kmHloc,Nlat,Nspin,Norb)
    hk          = zero
    hk(1:3:2,1:3:2) = hk11
    hk(2:4:2,2:4:2) = hk22
    !
  end function hk_haldane_timerevsym_model


  !--------------------------------------------------------------------!
  !reshaping functions:                                                !
  !--------------------------------------------------------------------!

  function nnn2nlso_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fin
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Fout(is,js) = Fin(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function nnn2nlso_reshape

  function lso2nnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fin
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Fout(ilat,ispin,jspin,iorb,jorb) = Fin(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn_reshape

end program ed_haldane



