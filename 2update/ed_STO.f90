
program ed_TEST_REPLICA
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
  integer                :: Nb,unit
  real(8),allocatable    :: Bath(:),Bath_(:)
  !dmft functions:
  complex(8),allocatable :: Delta(:,:,:,:,:),Delta_old(:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:),SigmaHk(:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  complex(8),allocatable :: impGw(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  complex(8),allocatable :: Ti3dt2g_Hloc(:,:),Ti3dt2g_Hloc_nn(:,:,:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:),kzgrid(:)
  !variables for the model:
  integer                :: Nk,Nkpath,i,j,iorb,jorb,io,jo,ispin,jspin
  real(8)                :: soc,ivb,wmixing,sumdens,xmu_old
  logical                :: surface,Hk_test,rotateG0loc
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  !custom variables:
  integer                :: conv_n_loop=0,shift_n_loop=0
  real(8)                :: Alvl
  logical                :: converged_n,upprshft
  !convergence functions:
  complex(8),allocatable :: delta_conv(:,:,:),delta_conv_avrg(:)
  !density matrix:
  real(8),allocatable    :: dm_eig(:)
  complex(8),allocatable :: density_matrix(:,:),dm_rot(:,:)
  !SOC expectations:
  complex(8),allocatable :: Stot(:,:,:),Ltot(:,:,:),jz(:)
  !rotated chempot shift:
  real(8)                :: bottom,top,shift,dw,nread_old
  real(8),allocatable    :: w(:)
  !
#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(LOGfile,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif
  !
  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,   "FINPUT",           default='inputED_REPL.in')
  call parse_input_variable(hkfile, "HKFILE",finput,    default="hkfile.in")
  call parse_input_variable(nk,     "NK",finput,        default=10)
  call parse_input_variable(nkpath, "NKPATH",finput,    default=500)
  call parse_input_variable(wmixing,"WMIXING",finput,   default=0.5d0)
  call parse_input_variable(soc,    "SOC",finput,       default=0.0d0)
  call parse_input_variable(ivb,    "IVB",finput,       default=0.0d0)
  call parse_input_variable(surface,"SURFACE",finput,   default=.false.)
  call parse_input_variable(Hk_test,"HK_TEST",finput,   default=.true.)
  call parse_input_variable(upprshft,    "upprshft",finput,   default=.false.)
  call parse_input_variable(rotateG0loc,"ROTATEG0loc",finput, default=.false.)
  !
  call ed_read_input(trim(finput))
  !
  Nso=Nspin*Norb
  !
  !Allocate dmft functions:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats));delta=zero
  allocate(delta_old(Nspin,Nspin,Norb,Norb,Lmats));delta_old=zero
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(SigmaHk(Nspin,Nspin,Norb,Norb,Lmats));SigmaHk=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(impGw(Nspin,Nspin,Norb,Norb,Lmats));impGw=zero
  !Allocate convergence functions:
  allocate(delta_conv(Nso,Nso,Lmats));delta_conv=zero
  allocate(delta_conv_avrg(Lmats));delta_conv_avrg=zero
  !Allocate density matrix:
  allocate(density_matrix(Nspin*Norb,Nspin*Norb));density_matrix=zero
  allocate(dm_eig(Nspin*Norb));dm_eig=zero
  allocate(dm_rot(Nspin*Norb,Nspin*Norb));dm_rot=zero
  !
  !Allocate SOC expectations:
  allocate(Stot(3,Norb,Norb));Stot=zero
  allocate(Ltot(3,Nspin,Nspin));Ltot=zero
  allocate(jz(3));jz=zero
  !
  !Buil the non interacting Hamiltonian:
  call build_hk(trim(hkfile))
  !
  !Allocate bath:
  if (bath_type/="replica") then
     Nb=get_bath_size()
  else
     Nb=get_bath_size(Ti3dt2g_Hloc_nn)
  endif
  if(ED_MPI_ID==0)write(LOGfile,*)"Bath_size:",Nb
  allocate(Bath(Nb));Bath=0.0d0
  allocate(Bath_(Nb));Bath_=0.0d0
  !
  !Setup solver:
  call ed_init_solver(bath,Ti3dt2g_Hloc_nn)
  !
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !MAIN LOOP:
     !
     call ed_solve(bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     call ed_get_gloc(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=3)
     call ed_get_weiss(Gmats,Smats,Delta,Ti3dt2g_Hloc_nn,iprint=3)
     if(ED_MPI_ID==0)write(LOGfile,'(a10,F10.5,a10,i3)') " wmixing",wmixing,"cg_weight",cg_weight
     if(iloop>1)then
        delta = wmixing*delta + (1.d0-wmixing)*delta_old
        !bath = wmixing*bath_ + (1.d0-wmixing)*bath_
     endif
     if (ed_mode=="normal") then
        call ed_chi2_fitgf(delta,bath,ispin=1)
        call spin_symmetrize_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(delta,bath)
     endif
     delta_old = delta
     bath_=bath
     !
     !operations:
     !
     if(bath_type=="replica")then
        Alvl=0.8d0
        if(ED_MPI_ID==0)then
           call build_hk_path
           call ed_get_density_matrix(density_matrix,dm_eig,dm_rot)
           call ed_get_quantum_SOC_operators(Stot,Ltot,jz)
           Stot=zero;Ltot=zero;jz=zero
           call Quantum_operator()
        endif
        call Jz_rotate(Greal,"Gr","A",bottom,top,pi*Alvl)
        call Jz_rotate(Smats,"Sw","A")
        call ed_get_gimp_matsubara(impGw)
        call Jz_rotate(impGw,"Gw","A")
     endif
     !
     !chemical potential find:
     !
     converged_n=.true.
     xmu_old=xmu
     sumdens=sum(ed_get_dens())
     if(ED_MPI_ID==0)write(*,'(3(a10,F10.5))') "sumdens",sumdens,"diffdens",abs(nread-sumdens),"nread",nread
     if(nread/=0.d0)then
        converged_n=.false.
        if(iloop>3)call search_chempot(xmu,sumdens,converged_n,Bath)
        if(ED_MPI_ID==0)write(*,'(2(a10,F10.5))') "xmu_old",xmu_old,"xmu_new",xmu
     endif
     if(converged_n)then
        conv_n_loop=conv_n_loop+1
     else
        conv_n_loop=0
     endif
     !
     !convergence:
     !
     do i=1,Lmats
        delta_conv(:,:,i)=nn2so_reshape(delta(:,:,:,:,i))
        delta_conv_avrg(i)=sum(delta_conv(:,:,i))
     enddo
     if(ED_MPI_ID==0) then
        converged = check_convergence(delta_conv_avrg,dmft_error,nsuccess,nloop)
        write(LOGfile,'(2(a15,L3))') "converged",converged,"converged(n)",converged_n
        converged = converged .and. converged_n
        write(LOGfile,'(a25,L3)') "total converged",converged
        write(LOGfile,'(a25,I3)') "iloop",iloop
        write(LOGfile,'(a25,I3)') "conv_n_loop",conv_n_loop
        write(LOGfile,'(a25,I3)') "shift_n_loop",shift_n_loop
     endif
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
     call mpi_barrier(MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     !
     !final mu shift:
     !
     shift_n_loop=0
     if(converged_n.and.upprshft.and.((nread==5.d0).or.(nread==2.d0)))then
        shift_n_loop=shift_n_loop+1
        if(bath_type/="replica")then
           if(allocated(w))deallocate(w);allocate(w(Lreal));w=0.0d0
           w = linspace(wini,wfin,Lreal,mesh=dw)
           loop1: do i=1,Lreal
              if(abs(aimag(Greal(1,1,1,1,i))).gt.0.8d0)then
                 bottom=w(i)
                 exit loop1
              endif
           enddo loop1
           loop2: do i=1,Lreal
              if(abs(aimag(Greal(1,1,1,1,Lreal-i+1))).gt.0.8d0)then
                 top=w(Lreal-i+1)
                 exit loop2
              endif
           enddo loop2
        endif
        if(ED_MPI_ID==0)write(LOGfile,*)"top",top,"bottom",bottom
        shift      = bottom + ( top - bottom ) / 2.d0
        xmu_old    = xmu
        if(abs(shift)>=0.005)then
           xmu        = xmu_old + shift
           converged  = .false.
           nread  = 0.0d0!con questo una volta che comincio a shiftare rigidamente la densità non la controllo piu
        endif
        if(ED_MPI_ID==0)then
           write(LOGfile,'(5(a10,F10.5))') "shift",shift,"xmu_old",xmu_old,"xmu_new",xmu
           unit=free_unit()
           open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
           write(unit,*)xmu,sumdens,sumdens-nerr,"shift"
           close(unit)
        endif
     endif
     !
     !
     if(ED_MPI_ID==0)call end_loop
  enddo
#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif
contains


!_______________________________________________________________________
!                      NON-INTERACTING HAMILTONIAN
!_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: H(k) file for main program and write G0_loc
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky,kz    
    integer                             :: io,jo
    integer                             :: iorb,jorb,ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Nso,Nso,Lmats) :: Gmats
    complex(8),dimension(Nso,Nso,Lreal) :: Greal
    complex(8),allocatable              :: G_in(:,:,:,:,:)
    real(8)                             :: wm(Lmats),wr(Lreal),dw
    complex(8),dimension(Nso,Nso)       :: inv_impHloc
    !
    if(allocated(Bath))stop " H(K) must be build before bath allocation, errors shall come otherwise"
    !
    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) for STO:"
    !
    !Alloco quantità
    !
    if(surface) then
       Lk=Nk**2
    else
       Lk=Nk**3
    endif
    if(ED_MPI_ID==0)write(LOGfile,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(LOGfile,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Nso,Nso,Lk));allocate(wtk(Lk));allocate(kxgrid(Nk),kygrid(Nk),kzgrid(Nk))
    wtk = 1.0d0/Lk
    kxgrid=0.0d0;kygrid=0.0d0;kzgrid=0.0d0
    SigmaHk=zero
    kxgrid = kgrid(Nk)
    kygrid = kgrid(Nk)
    if(.not.surface) kzgrid = kgrid(Nk)
    !
    !scrivo H(k)
    !
    if(surface) then
       Hk = build_hk_model(hk_Ti3dt2g,Nso,kxgrid,kygrid,[0d0])
       if(ED_MPI_ID==0.AND.present(file)) call write_hk_w90(file,Nso,Nd=Norb,Np=1,Nineq=1,hk=Hk,kxgrid=kxgrid,kygrid=kxgrid,kzgrid=[0d0] )
    else
       Hk = build_hk_model(hk_Ti3dt2g,Nso,kxgrid,kygrid,kzgrid)
       if(ED_MPI_ID==0.AND.present(file)) call write_hk_w90(file,Nso,Nd=Norb,Np=1,Nineq=1,hk=Hk,kxgrid=kxgrid,kygrid=kxgrid,kzgrid=kzgrid)
    endif
    !
    !calcolo impHloc = Sum_k [ H(k) ]
    !
    allocate(Ti3dt2g_Hloc(Nso,Nso));                 Ti3dt2g_Hloc=zero
    allocate(Ti3dt2g_Hloc_nn(Nspin,Nspin,Norb,Norb));Ti3dt2g_Hloc_nn=zero
    Ti3dt2g_Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs((Ti3dt2g_Hloc))<1.d-9)Ti3dt2g_Hloc=0d0
    Ti3dt2g_Hloc_nn=so2nn_reshape(Ti3dt2g_Hloc)
    !
    !scrivo impHloc
    !
    if(ED_MPI_ID==0) then
       call write_Hloc(Ti3dt2g_Hloc)
       write(LOGfile,*)
       write(LOGfile,*) "Sum over k of H(k) nella versione A1"
       write(LOGfile,*) "real"
       do i=1,Nso
          write(LOGfile,'(6F10.4)') (real(Ti3dt2g_Hloc(i,j)),j=1,Nso)
       enddo
       write(LOGfile,*) "complex"
       do i=1,Nso
          write(LOGfile,'(6F10.4)') (dimag(Ti3dt2g_Hloc(i,j)),j=1,Nso)
       enddo
       write(LOGfile,*)
    endif
    !
    !Build the local GF in the spin-orbital Basis:
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i) , Hk(:,:,ik) )/Lk
       enddo
       do i=1,Lreal
          Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps),Hk(:,:,ik))/Lk
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed",wm,Gmats(io,jo,:))
                call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,dreal(Greal(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    if(rotateG0loc) then
       allocate(G_in(Nspin,Nspin,Norb,Norb,Lreal))
       G_in=zero
       do i=1,Lreal
          G_in(:,:,:,:,i)=so2nn_reshape(Greal(:,:,i))
       enddo
       call Jz_rotate(G_in,"G0","A")
       deallocate(G_in)
    endif
  end subroutine build_hk


  !---------------------------------------------------------------------
  !PURPOSE: GET STO HAMILTONIAN eventually corrected with Sigma(iw=0)
  !---------------------------------------------------------------------
  function hk_Ti3dt2g(kvec,N) result(hk)
    real(8),dimension(:)        :: kvec
    complex(8),dimension(N,N)   :: hk
    complex(8),dimension(2,2)   :: s_x,s_y,s_z,Hk_temp
    complex(8),dimension(2,2)   :: t_inter
    real(8)                     :: kx,ky,kz
    integer                     :: N,ndx
    real(8),allocatable         :: HoppingMatrix(:,:)
    real(8),dimension(0:6)      :: orbitalHopping
    !
    kx=kvec(1);ky=kvec(2);kz=kvec(3)
    !
    allocate(HoppingMatrix(Norb,0:6));HoppingMatrix=0.0d0
    call get_hopping(HoppingMatrix)
    !
    Hk=zero
    do i=1,Norb
       ndx=2*i-1
       if(Hk_test)then
          !spin degenerate H(k)
          Hk(ndx:ndx+1,ndx:ndx+1) = diagonal_orbital_dispersion(kx,ky,kz,HoppingMatrix(i,:),0.1d0)
       else
          Hk(ndx:ndx+1,ndx:ndx+1) = diagonal_orbital_dispersion(kx,ky,kz,HoppingMatrix(i,:),1.0d0)
       endif
    enddo
    !
    if(SOC/=zero.or.IVB/=zero)then
       Hk = Hk + SOC*H_LS(kx,ky,kz) + IVB*H_IVB(kx,ky,kz)
    endif
    !
    !correction with Sigma(iw=0)
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = ispin + (iorb-1)*Nspin
                jo = jspin + (jorb-1)*Nspin
                Hk(io,jo) = Hk(io,jo) + SigmaHk(ispin,jspin,iorb,jorb,1)
             enddo
          enddo
       enddo
    enddo
    !
    !A1 shape: [Norb*Norb]*Nspin
    Hk = Z2so_reshape(Hk)
    !
  end function hk_Ti3dt2g

  !---------------------------------------------------------------------
  !PURPOSE: 2x2 band structures
  !---------------------------------------------------------------------
  function diagonal_orbital_dispersion(kx,ky,kz,t,t0) result(hk)
    real(8),intent(in)                :: t0!scaling factor for bandwidth
    real(8),intent(in)                :: kx,ky,kz
    real(8),intent(in),dimension(0:6) :: t
    complex(8),dimension(2,2)         :: hk
    !perovskite dispersion
    hk = zero
    if (surface) then
       if(Hk_test)then
          hk(1,1) = t(0)-2.*t0*(cos(kx)+cos(ky)+1.d0)
          hk(2,2) = hk(1,1)
       else
          hk(1,1) = t(0)+(                       & !onsite
                    -2.*t(1)*cos(kx)             & !t_100
                    -2.*t(2)*cos(ky)             & !t_010
                    -1.*t(3)                     & !t_001
                    -2.*t(4)*cos(ky)             & !t_011
                    -2.*t(5)*cos(kx)             & !t_101
                    -4.*t(6)*cos(kx)*cos(ky))*t0   !t_110
          hk(2,2) = hk(1,1)
       endif
    else
       if(Hk_test)then
          hk(1,1) = t(0)+(                       & !onsite
                    -2.*t(1)*cos(kx)             & !t_100
                    -2.*t(2)*cos(ky)             & !t_010
                    -2.*t(3)*cos(kz))*t0           !t_001
          hk(2,2) = hk(1,1)
       else
          hk(1,1) = t(0)+(                       & !onsite
                    -2.*t(1)*cos(kx)             & !t_100
                    -2.*t(2)*cos(ky)             & !t_010
                    -2.*t(3)*cos(kz)             & !t_001
                    -4.*t(4)*cos(ky)*cos(kz)     & !t_011
                    -4.*t(5)*cos(kx)*cos(kz)     & !t_101
                    -4.*t(6)*cos(kx)*cos(ky))*t0   !t_110
          hk(2,2) = hk(1,1)
       endif
    endif
  end function diagonal_orbital_dispersion

  function H_LS(kx,ky,kz) result(hls)
    real(8),intent(in)                           :: kx,ky,kz
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: hls
    hls=zero
    hls(1:2,3:4) = -Xi * pauli_z / 2.
    hls(1:2,5:6) = +Xi * pauli_y / 2.
    hls(3:4,5:6) = -Xi * pauli_x / 2.
    !hermiticity
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          hls(j,i)=conjg(hls(i,j))
       enddo
    enddo
  end function H_LS

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




  !---------------------------------------------------------------------
  !PURPOSE: Build the hopping integrals in k-space for realistic bandstructure
  !---------------------------------------------------------------------
  subroutine get_hopping(T)
  real(8),dimension(Norb,0:6),intent(out)      ::  T
  real(8),dimension(3,0:6)                     ::  T_bulk,T_VACSTO,T_LAOSTO
  real(8)                                      ::  Eo,t1,t2,t3
  real(8)                                      ::  t_010_yz,t_001_yz
  real(8)                                      ::  t_100_zx,t_001_zx
  real(8)                                      ::  t_100_xy,t_010_xy,t_001_xy

  !pristine lattice
  Eo = 3.31
  t1 = 0.276536
  t2 = 0.031329
  t3 = 0.076842

  !lattice distortion
  t_010_yz = 0.232 !se c'è solo l'abbassamento del Ti questo dovrebbe essere uguale a t1, magari c'è anche altro dovuto ad LAO
  t_001_yz = 0.475

  t_100_zx = 0.232
  t_001_zx = 0.475

  t_100_xy = 0.286
  t_010_xy = 0.286
  t_001_xy = 0.03

  !  BULK STO
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
  ! VAC/STO (not used)
  T_VACSTO=T_bulk
  !
  !  LAO/STO
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
  if(surface) then
     T=T_LAOSTO(1:Norb,:)
  else
     if(.not.Hk_test)then
        T=0.d0
        T=T_bulk(1:Norb,:)
     elseif(bath_type=="replica".and.Hk_test)then
        T=1.0d0
        T(1,0) = 0.0d0
        T(2,0) = 0.0d0
        T(3,0) = 0.0d0
     elseif(bath_type/="replica".and.Hk_test)then
        T=1.0d0
        T(1,0) = +0.30d0
        T(2,0) = -0.15d0
        T(3,0) = -0.15d0
     endif
  endif
  !
  end subroutine get_hopping

  !---------------------------------------------------------------------
  !PURPOSE: Get the STO(bulk) Hamiltonian along path
  !---------------------------------------------------------------------
  subroutine build_hk_path()
    integer                            :: i,j
    integer                            :: Npts
    real(8),dimension(:,:),allocatable :: kpath
    if(.not.surface)then
       if(ED_MPI_ID==0)then
          write(LOGfile,*)"Build bulk H(k) along the path M-R-G-M-X-G-X"
          write(LOGfile,*)
       endif
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
       if(ED_MPI_ID==0) then
          call solve_Hk_along_BZpath(hk_Ti3dt2g,Nso,kpath,Lk,colors_name=[character(len=20) :: 'red','green','blue','red','green','blue'],&
               points_name=[character(len=20) :: 'M', 'R', 'G', 'M', 'X', 'G', 'X'],file="Eigenband_bulk.nint")
          SigmaHk=Smats
          call solve_Hk_along_BZpath(hk_Ti3dt2g,Nso,kpath,Lk,colors_name=[character(len=20) :: 'red','green','blue','red','green','blue'],&
               points_name=[character(len=20) :: 'M', 'R', 'G', 'M', 'X', 'G', 'X'],file="Eigenband_bulk.Sigma")
          SigmaHk=zero
       endif
    else
       if(ED_MPI_ID==0)then
          if(ED_MPI_ID==0)write(LOGfile,*)"Build surface H(k) along the path M-X-G-X"
          write(LOGfile,*)
       endif
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_M1
       kpath(2,:)=kpoint_X1
       kpath(3,:)=kpoint_Gamma
       kpath(4,:)=kpoint_X1
       if(ED_MPI_ID==0) then
          call solve_Hk_along_BZpath(hk_Ti3dt2g,Nso,kpath,Lk,colors_name=[character(len=20) :: 'red','green','blue','red','green','blue'],&
               points_name=[character(len=20) :: 'M', 'X', 'G', 'X'],file="Eigenband_surf.nint")
          SigmaHk=Smats
          call solve_Hk_along_BZpath(hk_Ti3dt2g,Nso,kpath,Lk,colors_name=[character(len=20) :: 'red','green','blue','red','green','blue'],&
               points_name=[character(len=20) :: 'M', 'X', 'G', 'X'],file="Eigenband_surf.Sigma")
          SigmaHk=zero
       endif
    endif
  end subroutine build_hk_path



!_______________________________________________________________________
!                                    Gfs
!_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: G0_loc functions
  !---------------------------------------------------------------------
  function inverse_g0k(iw,hk) result(g0k)
    implicit none
    complex(8)                                    :: iw
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: hk
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: g0k,g0k_tmp
    integer                                       :: i,ndx
    integer (kind=4), dimension(6)                :: ipiv
    integer (kind=1)                              :: ok
    integer (kind=4), parameter                   :: lwork=2000
    complex (kind=8), dimension(lwork)            :: work
    real    (kind=8), dimension(lwork)            :: rwork
    !
    g0k=zero;g0k_tmp=zero
    !
    g0k=iw*eye(Nspin*Norb)-hk
    g0k_tmp=g0k
    !
    call inv(g0k)
    call inversion_test(g0k,g0k_tmp,1.e-9)
  end function inverse_g0k

  !---------------------------------------------------------------------
  !PURPOSE: rotations on G0loc/Gloc
  !---------------------------------------------------------------------
  subroutine Jz_rotate(Gsowr,type_funct,type_rot,bottom_,top_,lvl_)
    implicit none
    complex(8),allocatable,intent(in)             ::   Gsowr(:,:,:,:,:)
    character(len=2),intent(in)                   ::   type_funct
    character(len=1),intent(in),optional          ::   type_rot
    real(8),intent(out),optional                  ::   bottom_,top_
    real(8),intent(in),optional                   ::   lvl_
    complex(8),allocatable                        ::   G_in(:,:,:),G_out(:,:,:),Gimp(:,:,:,:,:)
    complex(8),allocatable                        ::   Luttinger(:)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   ::   theta_C,theta_R,impHloc_rot
    integer                                       ::   io,jo,ndx
    integer                                       ::   ispin,jspin,iorb,jorb
    integer                                       ::   Lfreq
    real(8)                                       ::   wr(Lreal),wm(Lmats),z_rot(Nspin*Norb)
    real(8)                                       ::   dw,bttm,tp,lvl,fact
    character(len=13)                             ::   file_rotation
    integer                                       ::   isetup=0
    !
    if(type_funct=="G0") isetup=1
    if(type_funct=="Gr") isetup=2
    if(type_funct=="Sw") isetup=3
    if(type_funct=="Gw") isetup=4
    lvl=1.0d0;if(present(lvl_))lvl=lvl_
    !
    call build_rotation(theta_C,impHloc_rot)
    !
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Lfreq=size(Gsowr,dim=5)
    if(allocated( G_in))deallocate( G_in);allocate( G_in(Nspin*Norb,Nspin*Norb,Lfreq));G_in=zero
    if(allocated(G_out))deallocate(G_out);allocate(G_out(Nspin*Norb,Nspin*Norb,Lfreq));G_out=zero
    !
    if(isetup==1) then
       if(ED_MPI_ID==0)write(LOGfile,*) "  A(w) rotation - non interacting system"
       fact=1.d0
    elseif(isetup==2) then
       if(ED_MPI_ID==0)write(LOGfile,*) "  A(w) rotation - interacting system"
       fact=1.d0
    elseif(isetup==3) then
       if(ED_MPI_ID==0)write(LOGfile,*) "  Sigma(iw) rotation"
       fact=1.d0
    elseif(isetup==4) then
       if(ED_MPI_ID==0)write(LOGfile,*) "  Gimp(iw) rotation"
       if(allocated(Luttinger))deallocate(Luttinger);allocate(Luttinger(2*Nspin*Norb));Luttinger=zero
       if(allocated(Gimp))deallocate(Gimp);allocate(Gimp(Nspin,Nspin,Norb,Norb,Lfreq));Gimp=zero
       call ed_get_gimp_matsubara(Gimp)
    endif
    !
    !
    !###############################################################
    !#                                                             #
    !#                      REFERENCE SYSTEM                       #
    !#                                                             #
    !###############################################################
    !
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                if(isetup/=4)then
                   G_in(io,jo,:)=Gsowr(ispin,jspin,iorb,jorb,:)
                else
                   G_in(io,jo,:)=Gimp(ispin,jspin,iorb,jorb,:)
                   if(io==jo)Luttinger(io)=Gimp(ispin,jspin,iorb,jorb,1)
                endif
             enddo
          enddo
       enddo
    enddo
    !save the integral
    if(isetup==1) then
       open(unit=106,file='sum_w_G0loc.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==2) then
       open(unit=106,file='sum_w_Gloc.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==3) then
       open(unit=106,file='sum_w_Smats.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==4) then
       open(unit=106,file='sum_w_Gmats.dat',status='unknown',action='write',position='rewind')
    endif
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(G_in(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    !
    !###############################################################
    !#                                                             #
    !#                   ROTATION WITH impHloc                     #
    !#                                                             #
    !###############################################################
    !
    if(present(type_rot).and.(type_rot=="H"))then
       if(ED_MPI_ID==0)write(LOGfile,*) "  Rotation with impHloc"
    !
    !1)rotation
    G_out=zero
    do i=1,Lfreq
       G_out(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(G_in(:,:,i),impHloc_rot))
    enddo
    !
    !Zqp save
    if(isetup==3)then
       open(unit=106,file='Zqp_rot.dat',status='unknown',action='write',position='append')
       do io=1,Nspin*Norb
          z_rot(io)   = 1.d0/( 1.d0 + abs( dimag(G_out(io,io,1))/(pi/beta) ))
       enddo
       rewind(106)
       write(106,*) "# impHloc rotation"
       write(106,'(90F15.9,1X)')(z_rot(io),io=1,Nspin*Norb)
       write(106,*)
       write(106,*)
       write(106,*)
       close(106)
    endif
    !
    !2)output save
    if(isetup==1) then
       file_rotation="G0loc_rot_H_l"
    elseif(isetup==2) then
       file_rotation="Giloc_rot_H_l"
    elseif(isetup==3) then
       file_rotation="Smats_rot_H_l"
    elseif(isetup==4) then
       file_rotation="Gmats_rot_H_l"
    endif
    if(isetup<=2)then
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   call splot(file_rotation//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(G_out(io,jo,:))/pi,dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    else
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   call splot(file_rotation//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed",wm,dimag(G_out(io,jo,:)),dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    endif
    !3)save the integral
    if(isetup==1) then
       open(unit=106,file='sum_w_G0loc_rot_H.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==2) then
       open(unit=106,file='sum_w_Gloc_rot_H.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==3) then
       open(unit=106,file='sum_w_Smats_rot_H.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==4) then
       open(unit=106,file='sum_w_Gmats_rot_H.dat',status='unknown',action='write',position='rewind')
    endif
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(G_out(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    close(106)
    endif
    !
    !
    !###############################################################
    !#                                                             #
    !#                    ROTATION WITH theta_C                    #
    !#                                                             #
    !###############################################################
    !
    if(present(type_rot).and.(type_rot=="A"))then
       if(ED_MPI_ID==0)write(LOGfile,*) "  Rotation with LS(C) Martins"
    !
    !1)rotation
    G_out=zero
    do i=1,Lfreq
       G_out(:,:,i)=matmul(transpose(conjg(theta_C)),matmul(G_in(:,:,i),theta_C))
    enddo
    !
    !Zqp save
    if(isetup==3)then
       open(unit=106,file='Zqp_rot.dat',status='unknown',action='write',position='append')
       do io=1,Nspin*Norb
          z_rot(io)   = 1.d0/( 1.d0 + abs( dimag(G_out(io,io,1))/(pi/beta) ))
       enddo
       rewind(106)
       write(106,'(90A15,1X)') "#J=1/2,jz=-1/2","#J=3/2,jz=+1/2","#J=3/2,jz=-3/2","#J=1/2,jz=+1/2","#J=3/2,jz=+3/2","#J=3/2,jz=-1/2"
       write(106,'(90F15.9,1X)')(z_rot(io),io=1,Nspin*Norb)
       write(106,*)
       write(106,*)
       write(106,*)
       close(106)
    endif
    !
    !Luttinger save
    if(isetup==4)then
       open(unit=106,file='Luttinger.dat',status='unknown',action='write',position='append')
       do io=1,Nspin*Norb
          Luttinger(Nspin*Norb+io)=G_out(io,io,1)  
       enddo
       rewind(106)
       write(106,'(90F15.9,1X)') (real(Luttinger(io)),io=1,2*Nspin*Norb),(aimag(luttinger(io)),io=1,2*Nspin*Norb)
       close(106)
    endif
    !
    !top-bottom find of the upper band
    if(nread==5.d0)ndx=1
    if(nread==2.d0)ndx=2
    if(present(top_).and.present(bottom_))then
       outerloop1:do i=1,Lfreq
          if(abs(aimag(G_out(ndx,ndx,i))).gt.lvl)then
             bottom_=wr(i)
             exit outerloop1
          endif
       enddo outerloop1
       outerloop2:do i=1,Lfreq
          if(abs(aimag(G_out(ndx,ndx,Lfreq-i+1))).gt.lvl)then
             top_=wr(Lfreq-i+1)
             exit outerloop2
          endif
       enddo outerloop2
    endif
    !
    !2)output save
    if(isetup==1) then
       file_rotation="G0loc_rot_A_l"
    elseif(isetup==2) then
       file_rotation="Giloc_rot_A_l"
    elseif(isetup==3) then
       file_rotation="Smats_rot_A_l"
    elseif(isetup==4) then
       file_rotation="Gmats_rot_A_l"
    endif
    if(isetup<=2)then
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   call splot(file_rotation//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(G_out(io,jo,:))/pi,dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    else
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   call splot(file_rotation//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed",wm,dimag(G_out(io,jo,:)),dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    endif
    !3)save the integral
    if(isetup==1) then
       open(unit=106,file='sum_w_G0loc_rot_A.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==2) then
       open(unit=106,file='sum_w_Gloc_rot_A.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==3) then
       open(unit=106,file='sum_w_Smats_rot_A.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==4) then
       open(unit=106,file='sum_w_Gmats_rot_A.dat',status='unknown',action='write',position='rewind')
    endif
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(G_out(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    close(106)
    endif
    !
    !
    !###############################################################
    !#                                                             #
    !#                  ROTATION WITH rot_rho                      #
    !#                                                             #
    !###############################################################
    !
    if(sum(abs(dm_rot))>=1e-6)then
       if(ED_MPI_ID==0)write(LOGfile,*) "  Rotation with impurity density matrix"
    !
    !1)rotation
    G_out=zero
    do i=1,Lfreq
       G_out(:,:,i)=matmul(transpose(conjg(dm_rot)),matmul(G_in(:,:,i),dm_rot))
    enddo
    !
    !Zqp save
    if(isetup==3)then
       open(unit=106,file='Zqp_rot.dat',status='unknown',action='write',position='append')
       do io=1,Nspin*Norb
          z_rot(io)   = 1.d0/( 1.d0 + abs( dimag(G_out(io,io,1))/(pi/beta) ))
       enddo
       write(106,*) "# density matrix rotation"
       write(106,'(90F15.9,1X)')(z_rot(io),io=1,Nspin*Norb)
       write(106,*)
       write(106,*)
       write(106,*)
       close(106)
    endif
    !
    !2)output save
    if(isetup==1) then
       file_rotation="G0loc_rot_R_l"
    elseif(isetup==2) then
       file_rotation="Giloc_rot_R_l"
    elseif(isetup==3) then
       file_rotation="Smats_rot_R_l"
    elseif(isetup==4) then
       file_rotation="Gmats_rot_R_l"
    endif
    if(isetup<=2)then
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   call splot(file_rotation//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(G_out(io,jo,:))/pi,dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    else
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   call splot(file_rotation//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed",wm,dimag(G_out(io,jo,:)),dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    endif
    !3)save the integral
    if(isetup==1) then
       open(unit=106,file='sum_w_G0loc_rot_R.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==2) then
       open(unit=106,file='sum_w_Gloc_rot_R.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==3) then
       open(unit=106,file='sum_w_Smats_rot_R.dat',status='unknown',action='write',position='rewind')
    elseif(isetup==4) then
       open(unit=106,file='sum_w_Gmats_rot_R.dat',status='unknown',action='write',position='rewind')
    endif
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(G_out(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    close(106)
    endif
    !
  end subroutine Jz_rotate


!_______________________________________________________________________
!                    Operators & Operations related to SOC
!_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: Build the rotations
  !---------------------------------------------------------------------
  subroutine build_rotation(theta_C_,impHloc_rot_)
    complex(8),dimension(6,6),intent(out)            ::   theta_C_
    complex(8),dimension(6,6),intent(out)            ::   impHloc_rot_
    real(8),dimension(6)                             ::   impHloc_eig
    theta_C_=zero
    !J=1/2 jz=-1/2
    theta_C_(1,1)=-Xi
    theta_C_(3,1)=-1.0d0
    theta_C_(6,1)=+Xi
    theta_C_(:,1)=theta_C_(:,1)/sqrt(3.)
    !J=1/2 jz=+1/2
    theta_C_(2,2)=-Xi
    theta_C_(4,2)=+1.0d0
    theta_C_(5,2)=-Xi
    theta_C_(:,2)=theta_C_(:,2)/sqrt(3.)
    !J=3/2 jz=-3/2
    theta_C_(2,3)=-Xi
    theta_C_(4,3)=+1.0d0
    theta_C_(5,3)=+2.0d0*Xi
    theta_C_(:,3)=theta_C_(:,3)/sqrt(6.)
    !J=3/2 jz=-1/2
    theta_C_(1,4)=+Xi
    theta_C_(3,4)=-1.0d0
    theta_C_(:,4)=theta_C_(:,4)/sqrt(2.)
    !J=3/2 jz=+1/2
    theta_C_(2,5)=-Xi 
    theta_C_(4,5)=-1.0d0
    theta_C_(:,5)=theta_C_(:,5)/sqrt(2.)
    !J=3/2 jz=+3/2
    theta_C_(1,6)=+Xi
    theta_C_(3,6)=+1.0d0
    theta_C_(6,6)=+2.0d0*Xi
    theta_C_(:,6)=theta_C_(:,6)/sqrt(6.)
    theta_C_=Z2so_reshape(theta_C_)
    !
    impHloc_rot_=zero
    impHloc_rot_=Ti3dt2g_Hloc
    call matrix_diagonalize(impHloc_rot_,impHloc_eig,'V','U')
    !
  end subroutine build_rotation




!_______________________________________________________________________
!                            reshape functions
!_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: reshape functions
  !  Z  = [Nspin,Nspin]*Norb
  !  A1 = [Norb*Norb]*Nspin
  !  A2 = [Nspin,Nspin,Norb,Norb]
  !---------------------------------------------------------------------
  function Z2so_reshape(fg) result(g)
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: fg
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin
    integer                                         :: io1,jo1,io2,jo2
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   !O-index
                   io1 = iorb + (ispin-1)*Norb
                   jo1 = jorb + (jspin-1)*Norb
                   !I-index
                   io2 = ispin + (iorb-1)*Nspin
                   jo2 = jspin + (jorb-1)*Nspin
                   !switch
                   g(io1,jo1)  = fg(io2,jo2)
                   !
                enddo
             enddo
          enddo
       enddo
  end function Z2so_reshape

  function so2Z_reshape(fg) result(g)
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: fg
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin
    integer                                         :: io1,jo1,io2,jo2
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   !O-index
                   io1 = ispin + (iorb-1)*Nspin
                   jo1 = jspin + (jorb-1)*Nspin
                   !I-index
                   io2 = iorb + (ispin-1)*Norb
                   jo2 = jorb + (jspin-1)*Norb
                   !switch
                   g(io1,jo1)  = fg(io2,jo2)
                   !
                enddo
             enddo
          enddo
       enddo
  end function so2Z_reshape

  function so2nn_reshape(fg) result(g)
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb)     :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin,io,jo
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   g(ispin,jspin,iorb,jorb)  = fg(io,jo)
                enddo
             enddo
          enddo
       enddo
  end function so2nn_reshape

  function nn2so_reshape(fg) result(g)
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: g
    complex(8),dimension(Nspin,Nspin,Norb,Norb)     :: fg
    integer                                         :: i,j,iorb,jorb,ispin,jspin,io,jo
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   g(io,jo) = fg(ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
  end function nn2so_reshape











































  !---------------------------------------------------------------------
  !PURPOSE: 
  !---------------------------------------------------------------------
  subroutine Quantum_operator()
    implicit none
    complex(8),allocatable             :: Gso(:,:,:,:,:)
    complex(8),allocatable             :: Stot(:,:,:),Ltot(:,:,:)
    complex(8)                         :: LdotS,intG,intP
    complex(8)                         :: Sx,Lx,Sy,Ly,Sz,Lz
    complex(8)                         :: j_aplha(3),jz
    integer                            :: ilat,io,jo
    integer                            :: ispin,jspin
    integer                            :: iorb,jorb
    real(8)                            :: wm(Lmats),wr(Lreal),dw
    real(8)                            :: site_mag(Norb)

    real(8)                                      ::  LS_R,LS_I
    real(8)                                      ::  Jzsq_R,Jzsq_I

    real(8)                                      ::  Lxsq,Lysq,Lzsq,Lsq
    real(8)                                      ::  Sxsq,Sysq,Szsq,Ssq
    real(8)                                      ::  jxsq,jysq,jzsq,Jsq
    integer :: unit_
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    allocate( Gso(Nspin,Nspin,Norb,Norb,Lmats)); Gso=zero
    call ed_get_gimp_matsubara(Gso)
    !
    !##############################################################
    !
    !                              S
    !
    !##############################################################
    !
    write(LOGfile,*) "Computing total Spin operator per orbital"
    write(LOGfile,*) "Lmats used:",Lmats
    allocate(Stot(3,Norb,Norb));Stot=zero
    !
    do iorb=1,Norb
       do jorb=1,Norb
          !Sx
          Stot(1,iorb,jorb)=0.5d0*sum(    (Gso(1,2,iorb,jorb,:)+Gso(2,1,iorb,jorb,:) ))/beta
          !Sy
          Stot(2,iorb,jorb)=0.5d0*sum( xi*(Gso(2,1,iorb,jorb,:)-Gso(1,2,iorb,jorb,:) ))/beta
          !Sz
          Stot(3,iorb,jorb)=0.5d0*sum(    (Gso(1,1,iorb,jorb,:)-Gso(2,2,iorb,jorb,:) ))/beta
       enddo
    enddo
    !
    site_mag=0.d0
    site_mag=ed_get_mag()
    !
    !##############################################################
    !
    !                              L
    !
    !##############################################################
    !
    write(LOGfile,*) "Computing total Orbital operator per spin"
    write(LOGfile,*) "Lmats used:",Lmats
    allocate(Ltot(3,Nspin,Nspin));Ltot=zero
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          !Lx
          Ltot(1,ispin,jspin)=sum(  xi*(Gso(ispin,jspin,3,2,:)-Gso(ispin,jspin,2,3,:))  )/beta
          !Ly
          Ltot(2,ispin,jspin)=sum(  xi*(Gso(ispin,jspin,1,3,:)-Gso(ispin,jspin,3,1,:))  )/beta
          !Lz
          Ltot(3,ispin,jspin)=sum( -xi*(Gso(ispin,jspin,1,2,:)-Gso(ispin,jspin,2,1,:))  )/beta
       enddo
    enddo
    !
    !##############################################################
    !
    !                              L.dot.S
    !
    !##############################################################
    !
    write(LOGfile,*) "Computing total L dot S operator per site"
    write(LOGfile,*) "Lmats used:",Lmats
    LdotS=zero
    LdotS =sum ( (Gso(1,1,3,1,:)-Gso(1,1,1,3,:))*xi &
                +(Gso(2,2,2,1,:)-Gso(2,2,1,2,:))*xi &
                -(Gso(1,2,2,2,:)+Gso(2,1,2,2,:))    &
                +(Gso(2,1,3,3,:)-Gso(1,2,3,3,:))*xi &
                -(Gso(2,1,3,1,:)+Gso(1,2,1,3,:))    &
                +(Gso(1,2,2,1,:)-Gso(2,1,1,2,:))*xi ) 

    !
    j_aplha(1)=trace(Stot(1,:,:))+trace(Ltot(1,:,:))
    j_aplha(2)=trace(Stot(2,:,:))+trace(Ltot(2,:,:))
    j_aplha(3)=trace(Stot(3,:,:))+trace(Ltot(3,:,:))
    !
    Sxsq = (trace(Stot(1,:,:)))*conjg(trace(Stot(1,:,:)))
    Sysq = (trace(Stot(2,:,:)))*conjg(trace(Stot(2,:,:)))
    Szsq = (trace(Stot(3,:,:)))*conjg(trace(Stot(3,:,:)))
    Ssq  = Sxsq + Sysq + Szsq
    !
    Lxsq = (trace(Ltot(1,:,:)))*conjg(trace(Stot(1,:,:)))
    Lysq = (trace(Ltot(2,:,:)))*conjg(trace(Ltot(2,:,:)))
    Lzsq = (trace(Ltot(3,:,:)))*conjg(trace(Ltot(3,:,:)))
    Lsq  = Lxsq + Lysq + Lzsq
    !
    jxsq = (j_aplha(1))*conjg(j_aplha(1))
    jysq = (j_aplha(2))*conjg(j_aplha(2))
    jzsq = (j_aplha(3))*conjg(j_aplha(3))
    Jsq  = jxsq + jysq + jzsq
    !

    !
    !   IMPURITY SPIN OPERATOR - (S)
    !
    unit_ = free_unit()
    open(unit=unit_,file='S_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(A)')"# Re{Sx_(iorb,jorb)}, Im{Sx_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))')(real(Stot(1,iorb,jorb)),jorb=1,Norb),(aimag(Stot(1,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Sy_(iorb,jorb)}, Im{Sy_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))')(real(Stot(2,iorb,jorb)),jorb=1,Norb),(aimag(Stot(2,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Sz_(iorb,jorb)}, Im{Sz_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))')(real(Stot(3,iorb,jorb)),jorb=1,Norb),(aimag(Stot(3,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Sx]}","2-Im{Tr[Sx]}","3-Re{Tr[Sy]}","4-Im{Tr[Sy]}","5-Re{Tr[Sz]}","6-Im{Tr[Sz]}","7-|Sx|^2","8-|Sy|^2","9-|Sz|^2","10-|S|^2"
    write(unit_,'(30(F20.12,1X))')real(trace(Stot(1,:,:))),aimag(trace(Stot(1,:,:))),&
                                   real(trace(Stot(2,:,:))),aimag(trace(Stot(2,:,:))),&
                                   real(trace(Stot(3,:,:))),aimag(trace(Stot(3,:,:))),&
                                   Sxsq,Sysq,Szsq,Ssq
    close(unit_)
    !
    !   IMPURITY ORBITAL ANGULAR MOMENTUM OPERATOR - (L)
    !
    unit_ = free_unit()
    open(unit=unit_,file='L_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(A)')"# Re{Lx_(ipin,jspin)}, Im{Lx_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))')(real(Ltot(1,ispin,jspin)),jspin=1,Nspin),(aimag(Ltot(1,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Ly_(ipin,jspin)}, Im{Ly_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))')(real(Ltot(2,ispin,jspin)),jspin=1,Nspin),(aimag(Ltot(2,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Lz_(ipin,jspin)}, Im{Lz_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))')(real(Ltot(3,ispin,jspin)),jspin=1,Nspin),(aimag(Ltot(3,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Lx]}","2-Im{Tr[Lx]}","3-Re{Tr[Ly]}","4-Im{Tr[ly]}","5-Re{Tr[Lz]}","6-Im{Tr[Lz]}","7-|Lx|^2","8-|Ly|^2","9-|Lz|^2","10-|L|^2"
    write(unit_,'(30(F20.12,1X))')real(trace(Ltot(1,:,:))),aimag(trace(Ltot(1,:,:))),&
                                   real(trace(Ltot(2,:,:))),aimag(trace(Ltot(2,:,:))),&
                                   real(trace(Ltot(3,:,:))),aimag(trace(Ltot(3,:,:))),&
                                   Lxsq,Lysq,Lzsq,Lsq
    close(unit_)
    !
    !   IMPURITY TOTAL ANGULAR MOMENTUM OPERATOR - (J = S + L)
    !
    unit_ = free_unit()
    open(unit=unit_,file='J_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))')"#1-Re{jx}","2-Im{jx}","3-Re{jy}","4-Im{jy}","5-Re{jz}","6-Im{jz}","7-|jx|^2","8-|jy|^2","9-|jz|^2","10-|j|^2","11-Re{L.S}","12-Im{L.S}"
    write(unit_,'(30(F20.12,1X))')real(j_aplha(1)),aimag(j_aplha(1)),real(j_aplha(2)),aimag(j_aplha(2)),real(j_aplha(3)),aimag(j_aplha(3))&
                                  ,jxsq,jysq,jzsq,Jsq,real(impLdotS),aimag(impLdotS)
    close(unit_)
    !




    !       TRIAL CALCULATION OF SOC RENORMALIZTION
    Gso=zero
    call ed_get_gimp_matsubara(Gso)
    intG=sum(Gso(1,1,1,1,:))
    intP=sum(Gso(1,2,1,3,:))
    Jzsq_R=real( 5.d0*intG-4.d0*intP)/beta
    Jzsq_I=aimag(5.d0*intG-4.d0*intP)/beta
    LS_R=real( 6.d0*intP)/beta
    LS_I=aimag(6.d0*intP)/beta
    open(unit=unit_,file='J_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))')"# 1-Re[jz_sq]", "2-Im[jz_sq]","3-Re[LS]","4-Im[LS]"
    write(unit_,'(30(F20.12,1X))')Jzsq_R,Jzsq_I,LS_R,LS_I
    close(unit_)

    !
    deallocate(Ltot,Stot)
    !
  end subroutine Quantum_operator







  !---------------------------------------------------------------------
  !PURPOSE: Inversion test
  !---------------------------------------------------------------------
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


end program ed_TEST_REPLICA

