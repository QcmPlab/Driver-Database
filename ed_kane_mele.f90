program ed_kanemele
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso,Nlat,Nineq
  logical                                       :: converged
  integer                                       :: ispin,ilat!,i,j

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

  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: e1,e2   !real-space lattice basis
  real(8),dimension(2)                          :: bk1,bk2 !reciprocal space lattice basis
  real(8),dimension(2)                          :: d1,d2,d3
  real(8),dimension(2)                          :: a1,a2,a3
  real(8),dimension(2)                          :: bklen

  !variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: t1,t2,phi,Mh,wmixing
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym,neelsym
  !
  real(8),allocatable,dimension(:)              :: dens
  integer                                     :: comm,rank
  logical                                     :: master

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputKANEMELE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(t1,"T1",finput,default=2d0,comment='NN hopping, fixes noninteracting bandwidth')
  call parse_input_variable(t2,"T2",finput,default=0d0,comment='Haldane-like NNN hopping-strenght, corresponds to lambda_SO in KM notation')
  call parse_input_variable(phi,"PHI",finput,default=pi/2d0,comment='Haldane-like flux for the SOI term, KM model corresponds to a pi/2 flux')
  call parse_input_variable(mh,"MH",finput,default=0d0, comment='On-site staggering, aka Semenoff-Mass term')
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0, comment='Mixing parameter: 0 means 100% of the old bath (no update at all), 1 means 100% of the new bath (pure update)')
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.,comment='T fits just one spin bath (and copies on the other); F fits both. For full spin freedom set to false and use ed_mode=nonsu2.')
  call parse_input_variable(neelsym,"NEELSYM",finput,default=.false.,comment='Refers to AFM ordering; currently has no real consequence on this driver.')
  !
  call ed_read_input(trim(finput),comm)
  !
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  if(neelsym.AND.spinsym)stop "Wrong setup from input file: NEELSYM=T not with SPINSYM=T"

  if(Norb/=1.OR.Nspin/=2)stop "Wrong setup from input file: Norb=1 AND Nspin=2 is the correct configuration for the model."
  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso                 !=4 = 2(ineq sites)*2(spin)*1(orb)


  !Lattice basis (a=1; a0=sqrt3*a) is:
  !e_1 = a0 [ sqrt3/2 , 1/2 ] = 3/2a[1, 1/sqrt3]
  !e_2 = a0 [ sqrt3/2 ,-1/2 ] = 3/2a[1,-1/sqrt3]
  e1 = 3d0/2d0*[1d0, 1d0/sqrt(3d0)]
  e2 = 3d0/2d0*[1d0,-1d0/sqrt(3d0)]

  !LATTICE BASIS: nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]

  !next nearest-neighbor displacements: A-->A, B-->B, cell basis
  a1 = d1-d3                    !3/2*a[1,1/sqrt3]
  a2 = d2-d3                    !3/2*a[1,-1/sqrt3]
  a3 = d1-d2



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


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     !Solve separately the two atoms:
     do ilat=1,Nlat
        call ed_solve(comm,Bath(ilat,:),Hloc(ilat,:,:,:,:))
        call ed_get_sigma_matsubara(Smats(ilat,:,:,:,:,:))
        call ed_get_sigma_realaxis(Sreal(ilat,:,:,:,:,:))
     enddo
     !
     ! Smats(2,2,2,:,:,:) = -Smats(1,1,1,:,:,:) !sub_B(dw,dw) = -sub_A(up,up)
     ! Smats(2,1,1,:,:,:) = -Smats(1,2,2,:,:,:) !sub_B(up,up) = -sub_A(dw,dw)
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     !
     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc)
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc)
     endif
     !
     !Fit the new bath, starting from the old bath + the supplied delta
     !Behaves differently depending on the ed_mode input:
     !IF(NORMAL): normal/magZ phase is solved, so either fit spin1 or spin1&2 -> SPINSYM of choice
     !IF(NONSU2): Sz-conservation broken, SOI or magXY, fit both spins -> SPINSYM has to be false
     !IF(SUPERC): gives error, superconductivity is not allowed here!
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
     !
     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !
     !Check convergence. This is now entirely MPI-aware:
     converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     !
     call end_loop
  enddo

  call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=4)


  !Extract and print retarded self-energy and Green's function 
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
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
  !PURPOSE: Get Kane Mele Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional          :: file
    integer                            :: i,j,ik
    integer                            :: ix,iy
    real(8)                            :: kx,ky
    real(8),dimension(2)               :: pointK,pointKp
    integer                            :: iorb,jorb
    integer                            :: isporb,jsporb
    integer                            :: ispin,jspin
    integer                            :: unit
    real(8),dimension(:,:),allocatable :: KPath
    !
    Lk= Nk*Nk
    write(LOGfile,*)"Build H(k) Kane-Mele:",Lk
    write(LOGfile,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    !
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    !
    !
    call TB_build_model(Hk,hk_kanemele_model,Nlso,[Nk,Nk],wdos=.false.)
    !
    !
    allocate(kmHloc(Nlso,Nlso))
    kmHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(kmHloc))<1.d-4)kmHloc=0d0
    if(master)call TB_write_Hloc(kmHloc)
    if(master)call TB_write_Hloc(kmHloc,'Hloc.txt')
    !
    !
    ! pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
    ! pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]
    ! if(master)then
    !    allocate(Kpath(4,2))
    !    KPath(1,:)=[0,0]
    !    KPath(2,:)=pointK
    !    Kpath(3,:)=pointKp
    !    KPath(4,:)=[0d0,0d0]
    !    call TB_Solve_model(hk_kanemele_model,Nlso,KPath,Nkpath,&
    !         colors_name=[red1,blue1,red1,blue1],&
    !         points_name=[character(len=10) :: "G","K","K`","G"],&
    !         file="Eigenbands.nint",iproject=.false.)
    ! endif
  end subroutine build_hk



  !--------------------------------------------------------------------!
  !Kane-Mele HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_kanemele_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(2,2)       :: hk11,hk22
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz
    !----- Old Version ---------
    !real(8)                         :: kdotd(3),kdota(3)
    !!(k.d_j)
    !kdotd(1) = dot_product(kpoint,d1)
    !kdotd(2) = dot_product(kpoint,d2)
    !kdotd(3) = dot_product(kpoint,d3)
    !!(k.a_j)
    !kdota(1) = dot_product(kpoint,a1)
    !kdota(2) = dot_product(kpoint,a2)
    !kdota(3) = dot_product(kpoint,a3)
    !
    !h0 = 2*t2*cos(phi)*sum( cos(kdota(:)) )
    !hx =-t1*sum( cos(kdotd(:)) )
    !hy =-t1*sum( sin(kdotd(:)) )
    !hz = 2*t2*sin(phi)*sum( sin(kdota(:)) )
    !
    !
    !----- New Version [Consistent with Bernevig&Hughes book] ---------
    ! This way of creating the Hamiltonian we get a Bloch structure and can obtain local
    ! (i.e. in the unit-cell) quantities by averaging over all k-points. 
    real(8)			     :: kdote1, kdote2
    !
    kdote1 = dot_product(kpoint,e1)
    kdote2 = dot_product(kpoint,e2)
    !
    h0 = 2*t2*cos(phi)*( cos(kdote1) + cos(kdote2) + cos(kdote1-kdote2) )
    hx = t1*( cos(kdote1) + cos(kdote2) + 1)
    hy = t1*( sin(kdote1) + sin(kdote2) )
    hz = 2*t2*sin(phi)*( sin(kdote1) - sin(kdote2) - sin(kdote1-kdote2) )
    !
    hk11 = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Mh*pauli_z
    !
    hk22 = h0*pauli_0 + hx*pauli_x + hy*pauli_y - hz*pauli_z + Mh*pauli_z
    !
    hk          = zero
    !hk(1:2,1:2) = hk11
    !hk(3:4,3:4) = hk22
    !create hk such that it works well with lso2nnn_reshape(kmHloc,Nlat,Nspin,Norb)
    hk(1:3:2,1:3:2) = hk11
    hk(2:4:2,2:4:2) = hk22
    !
  end function hk_kanemele_model



  !--------------------------------------------------------------------!
  !Reshaping functions:                                                !
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







end program ed_kanemele




