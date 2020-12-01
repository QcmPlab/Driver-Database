
program ed_SIO
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI_INEQ
  USE MPI
#endif
  implicit none
  integer                :: iloop,Lk,Nso
  logical                :: converged
  !Bath:
  integer                :: Nb,unit,Lstart
  real(8),allocatable    :: Bath(:,:),Bath_(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:,:),Sreal(:,:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:,:),Greal(:,:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  complex(8),allocatable :: Ti3dt2g_Hloc(:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:),kzgrid(:)
  !variables for the model:
  integer                :: Nk,Nkpath,i,j,iorb,jorb,io,jo,ispin,jspin,ilat
  real(8)                :: wmixing,dens_per_site,soc
  real(8),allocatable    :: orb_dens(:,:)
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  character(len=32)      :: HlocFILE
  logical                :: spinsym,IOfile,paraJz,intersite,rotateG0loc
  !convergence function
  complex(8),allocatable :: delta_conv(:,:,:,:),delta_conv_avrg(:)
  !rotation on impHloc
  complex(8),allocatable     :: impHloc_rot(:,:),analytic_rot(:,:)
  real(8),allocatable        :: impHloc_eig(:)

#ifdef _MPI_INEQ
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,       "FINPUT",             default='inputED_SIO.in')
  call parse_input_variable(HlocFILE,   "HlocFILE",finput,    default="Hloc.dat")
  call parse_input_variable(nk,         "NK",finput,          default=10)
  call parse_input_variable(wmixing,    "WMIXING",finput,     default=0.5d0)
  call parse_input_variable(soc,        "SOC",finput,         default=1.d0)
  call parse_input_variable(Nlat,       "NLAT",finput,        default=8)
  call parse_input_variable(paraJz,     "PARAJz",finput,      default=.false.)
  call parse_input_variable(intersite,  "INTERSITE",finput,   default=.true.)
  call parse_input_variable(rotateG0loc,"ROTATEG0loc",finput, default=.false.)

  call ed_read_input(trim(finput))

  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  !Allocate convergence funct
  allocate(delta_conv(Nlat,Nso,Nso,Lmats),delta_conv_avrg(Lmats))
  allocate(orb_dens(Nlat,Norb))
  !
  !Read the Hamiltonian
  call non_interacting_setup()
  !
  !stop
  !
  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nlat,Nb),Bath_(Nlat,Nb))
  call ed_init_solver_lattice(bath)
  !
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(mpiID==0)call start_loop(iloop,nloop,"DMFT-loop")
     !
     call ed_solve_lattice(bath,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc),iprint=3)
     !
     if(mpiID==1)call Quantum_operator()
     !
     call ed_get_sigma_matsubara_lattice(Smats,Nlat)
     call ed_get_sigma_real_lattice(Sreal,Nlat)
     !
     call ed_get_gloc_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=3)
     !
     call rotate_Gloc()
     !
     call ed_get_weiss_lattice(Gmats,Smats,Delta,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc),iprint=3)
     !
     if(paraJz) then
        delta(1,:,:,:,:,:)=(delta(1,:,:,:,:,:)+delta(2,:,:,:,:,:))/2
        delta(2,:,:,:,:,:)=(delta(1,:,:,:,:,:)+delta(2,:,:,:,:,:))/2
        delta(3,:,:,:,:,:)=(delta(3,:,:,:,:,:)+delta(4,:,:,:,:,:))/2
        delta(4,:,:,:,:,:)=(delta(3,:,:,:,:,:)+delta(4,:,:,:,:,:))/2
        delta(5,:,:,:,:,:)=(delta(5,:,:,:,:,:)+delta(6,:,:,:,:,:))/2
        delta(6,:,:,:,:,:)=(delta(5,:,:,:,:,:)+delta(6,:,:,:,:,:))/2
        delta(7,:,:,:,:,:)=(delta(7,:,:,:,:,:)+delta(8,:,:,:,:,:))/2
        delta(8,:,:,:,:,:)=(delta(7,:,:,:,:,:)+delta(8,:,:,:,:,:))/2
     endif
     !
     Bath_=bath
     call ed_chi2_fitgf_lattice(bath,delta,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc))
     if(iloop>1) Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath
     delta_conv=zero
     delta_conv_avrg=zero
     do i=1,Lmats
        do ilat=1,Nlat
           do ispin=1,Nspin
              do jspin=1,Nspin
                 do iorb=1,Norb
                    do jorb=1,Norb
                       if((ispin.eq.jspin).and.(iorb.eq.jorb)) then
                          io = iorb + (ispin-1)*Norb
                          jo = jorb + (jspin-1)*Norb
                          delta_conv(ilat,io,jo,i)=delta(ilat,ispin,jspin,iorb,jorb,i)
                          delta_conv_avrg(i)=delta_conv_avrg(i)+delta_conv(ilat,io,jo,i)
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     delta_conv_avrg=delta_conv_avrg/(Nso*Nlat)
     if(mpiID==0) converged = check_convergence(delta_conv_avrg,dmft_error,nsuccess,nloop)
     !if(mpiID==0) converged = check_convergence_global(delta_conv_avrg,dmft_error,nsuccess,nloop)
     !if(mpiID==0) converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
     !if(mpiID==0)converged = check_convergence_global(delta_conv(:,:,:),dmft_error,nsuccess,nloop)
     !
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif
     !
     orb_dens=ed_get_dens_lattice(Nlat)
     dens_per_site=sum(orb_dens)/Nlat
     if (mpiID==0) write(*,*) "dens_per_site",dens_per_site,"xmu",xmu,"converged",converged
     if(nread/=0.d0)call search_chemical_potential(xmu,dens_per_site,converged)
     if (mpiID==0) write(*,*) "dens_per_site",dens_per_site,"xmu",xmu,"converged",converged
     !
     if(mpiID==0)call end_loop
  enddo
#ifdef _MPI_INEQ
  call MPI_FINALIZE(mpiERR)
#endif
contains



  !_______________________________________________________________________
  !                            HAMILTONIAN
  !_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: H(k) file for main program and write G0_loc
  !---------------------------------------------------------------------
  subroutine non_interacting_setup(file)
    character(len=*),optional                         :: file
    integer                                           :: i,j,ik=0
    integer                                           :: ix,iy
    real(8)                                           :: kx,ky,kz    
    integer                                           :: io,jo
    integer                                           :: iorb,jorb,ispin,jspin,ilat,jlat
    integer                                           :: unit
    complex(8),dimension(Nlat*Nso,Nlat*Nso,Lmats)     :: Gmats,Gmats_rot
    complex(8),dimension(Nlat*Nso,Nlat*Nso,Lreal)     :: Greal,Greal_rot
    complex(8)                                        :: Gso(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                                           :: wm(Lmats),wr(Lreal),dw,xmu_0
    real(8)                                           :: dumR(Nlat*Nso,Nlat*Nso),dumI(Nlat*Nso,Nlat*Nso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)  :: Hloc_dum
    complex(8),allocatable                            :: site_matrix_in(:,:),site_matrix_out(:,:)
    complex(8),allocatable                            :: impHloc_app(:,:)
    !
    complex(8)                                        :: Stot(Nlat,Nlat,3,Norb,Norb)
    complex(8)                                        :: S_print(Nlat*Norb,Nlat*Norb)
    complex(8)                                        :: Sz(Norb,Norb)
    !
    complex(8)                                        :: Ltot(Nlat,Nlat,3,Nspin,Nspin)
    complex(8)                                        :: L_print(Nlat*Nspin,Nlat*Nspin)
    complex(8)                                        :: Lz(Nspin,Nspin)
    !
    integer                                           :: ndx
    real(8)                                           :: alf
    complex(8)                                        :: Jz_avrg
    complex(8)                                        :: Jz(Nlat*Nspin*Norb,Nlat*Nspin*Norb)    !in funzione di alfa
    complex(8)                                        :: Jz_rot(Nlat*Nspin*Norb,Nlat*Nspin*Norb)!rotazione con impHloc_rot

    !#######################################################
    !
    !   H(k) reading
    !
    !#######################################################
    !
    if(mpiID==0)write(LOGfile,*)"Read H(k) for SIO:"
    Lk=Nk
    if(mpiID==0)write(*,*)"# of k-points     :",Lk
    if(mpiID==0)write(*,*)"# of sites        :",Nlat
    if(mpiID==0)write(*,*)"# of SO-bands     :",Nso
    if(mpiID==0)write(*,*)"# of SOC factor   :",soc
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Nlat*Nso,Nlat*Nso,Lk));allocate(wtk(Lk))
    wtk = 1.0d0/Lk
    !
    open(unit=123,file='hk_2_ED.in',status='old',action='read')
    do ik=1,Lk
       do io=1,Nlat*Nspin*Norb
          read(123,'(50F10.5)') (dumR(io,jo),jo=1,Nlat*Nspin*Norb)
       enddo
       do io=1,Nlat*Nspin*Norb
          read(123,'(50F10.5)') (dumI(io,jo),jo=1,Nlat*Nspin*Norb)
       enddo
       Hk(:,:,ik)=cmplx(dumR,dumI)
    enddo
    if (soc .ne. 1.0d0) then
       write(*,*)"rescaling SOC"
       do io=1,Nlat*Nspin*Norb
          do jo=1,Nlat*Nspin*Norb
             if(io.ne.jo) Hk(io,jo,:)=Hk(io,jo,:)/soc
          enddo
       enddo
    endif
    close(123)
    !
    allocate(Ti3dt2g_Hloc(Nlat*Nso,Nlat*Nso))
    Ti3dt2g_Hloc = sum(Hk,dim=3)/Lk
    where(abs((Ti3dt2g_Hloc))<1.d-9)Ti3dt2g_Hloc=zero
    if(mpiID==0) then
       call write_Hloc(Ti3dt2g_Hloc,HlocFILE)
    endif
    !
    Hloc_dum=reshape_A1_to_A2_L(Ti3dt2g_Hloc)
    open(unit=100,file='impHloc.dat',status='unknown',action='write',position='rewind')
    if(mpiID==0) then
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      write(100,'(5I3,10F15.10)')ilat,ispin,jspin,iorb,jorb,real(Hloc_dum(ilat,ispin,jspin,iorb,jorb))
                   enddo
                enddo
             enddo
          enddo
       enddo
       write(100,*)
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      write(100,'(5I3,10F15.10)')ilat,ispin,jspin,iorb,jorb,aimag(Hloc_dum(ilat,ispin,jspin,iorb,jorb))
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    close(100)
    !
    !#######################################################
    !
    !   rotation
    !
    !#######################################################
    !
    if(allocated(impHloc_rot)) deallocate(impHloc_rot)
    allocate(impHloc_rot(Nlat*Nspin*Norb,Nlat*Nspin*Norb));impHloc_rot=zero
    if(allocated(impHloc_eig)) deallocate(impHloc_eig)
    allocate(impHloc_eig(Nlat*Nspin*Norb));impHloc_eig=0.d0
    !
    if (intersite) then
       !
       impHloc_rot=zero
       impHloc_rot=Ti3dt2g_Hloc
       call matrix_diagonalize(impHloc_rot,impHloc_eig,'V','U')
       !
       if(allocated(site_matrix_in)) deallocate(site_matrix_in)
       allocate(site_matrix_in(Nlat*Nspin*Norb,Nlat*Nspin*Norb));site_matrix_in=zero
       if(allocated(site_matrix_out)) deallocate(site_matrix_out)
       allocate(site_matrix_out(Nlat*Nspin*Norb,Nlat*Nspin*Norb));site_matrix_out=zero
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                site_matrix_in(io,io)=cmplx(float(ilat),0.0d0)
             enddo
          enddo
       enddo
       !mappa per tenere traccia almeno del contributo dei siti
       site_matrix_out=matmul(transpose(conjg(impHloc_rot)),matmul(site_matrix_in,impHloc_rot))
       !
       open(unit=104,file='site_mixing.dat',status='unknown',action='write',position='rewind')
       write(104,*)"Sites before rotation"
       do io=1,Nlat*Nspin*Norb
          write(104,'(100I12)')(nint(real(site_matrix_in(io,jo))),jo=1,Nlat*Nspin*Norb)
       enddo
       write(104,*)"Sites after rotation"
       write(104,*)"R:"
       do io=1,Nlat*Nspin*Norb
          write(104,'(100F12.4)')(real(site_matrix_out(io,jo)),jo=1,Nlat*Nspin*Norb)
       enddo
       write(104,*)"I:"
       do io=1,Nlat*Nspin*Norb
          write(104,'(100F12.4)')(aimag(site_matrix_out(io,jo)),jo=1,Nlat*Nspin*Norb)
       enddo
       close(104)
       !
    else
       !
       if(allocated(impHloc_app)) deallocate(impHloc_app)
       allocate(impHloc_app(Nspin*Norb,Nspin*Norb));impHloc_app=zero
       !
       impHloc_rot=zero
       do ilat=1,Nlat
          io=1+(ilat-1)*Nso
          jo=Nso+(ilat-1)*Nso
          impHloc_app=Ti3dt2g_Hloc(io:jo,io:jo)
          call matrix_diagonalize(impHloc_app,impHloc_eig(io:jo),'V','U')
          impHloc_rot(io:jo,io:jo)=impHloc_app
       enddo
       !
    endif
    !
    open(unit=102,file='impHloc_eig.dat',status='unknown',action='write',position='rewind')
    do ilat=1,Nlat
       write(102,'(1I3,20F25.20)')ilat,(impHloc_eig(io),io=1+(ilat-1),Nspin*Norb+(ilat-1))
    enddo
    close(102)
    !
    open(unit=103,file='impHloc_rot.dat',status='unknown',action='write',position='rewind')
    write(103,*)"impHloc rotation, Real Part"
    do io=1,Nlat*Nspin*Norb
       write(103,'(100F12.4)')(real(impHloc_rot(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    write(103,*)"impHloc rotation, Iaginary Part"
    do io=1,Nlat*Nspin*Norb
       write(103,'(100F12.4)')(aimag(impHloc_rot(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    close(103)
    !
    !#######################################################
    !
    !   operation on non-interacting Gfs
    !
    !#######################################################
    !
    if (rotateG0loc) then
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    xmu_0=Ti3dt2g_Hloc(1,1)
    !wr=wr+xmu_0
    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i)+xmu_0,Hk(:,:,ik) )/Lk
       enddo
       do i=1,Lreal
          Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu_0,Hk(:,:,ik))/Lk
       enddo
    enddo
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_iw.ed",wm,dimag(Gmats(io,jo,:)),real(Gmats(io,jo,:)) )
                   call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,real(Greal(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    enddo
    open(unit=106,file='sum_w_G0loc.dat',status='unknown',action='write',position='rewind')
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      write(106,*) io,jo,"---",ilat,jlat,ispin,jspin,iorb,jorb,sum(abs(Greal(io,jo,:)))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    !
    do i=1,Lmats
       Gmats_rot(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(Gmats(:,:,i),impHloc_rot))
    enddo
    do i=1,Lreal
       Greal_rot(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(Greal(:,:,i),impHloc_rot))
    enddo
    !
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   call splot("G0loc_rot_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_iw.ed",wm,dimag(Gmats_rot(io,jo,:)),real(Gmats_rot(io,jo,:)) )
                   call splot("G0loc_rot_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_realw.ed",wr,-dimag(Greal_rot(io,jo,:))/pi,real(Greal_rot(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    enddo
    open(unit=105,file='sum_w_G0loc_rot.dat',status='unknown',action='write',position='rewind')
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      write(105,*) io,jo,"---",ilat,jlat,ispin,jspin,iorb,jorb,sum(abs(Greal_rot(io,jo,:)))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    close(105)
    endif
    !
  end subroutine non_interacting_setup


  !_______________________________________________________________________
  !                                    Gfs
  !_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: G0_loc functions DA RIFARE ATTENZIONE CHE H(k) è nella forma A1
  !---------------------------------------------------------------------
  function inverse_g0k(iw,hk) result(g0k)
    implicit none
    complex(8)                                              :: iw
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)   :: hk
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)   :: g0k
    !
    g0k=zero
    g0k=iw*eye(Nlat*Nspin*Norb)-hk
    !
    call inv(g0k)
    !
  end function inverse_g0k


  !---------------------------------------------------------------------
  !PURPOSE: rotaizione delle Gimp per portarle in una base simile a J,jz
  !---------------------------------------------------------------------
  subroutine rotate_Gloc()
    implicit none
    complex(8),allocatable             :: Gsowr(:,:,:,:,:,:,:),Gsoiw(:,:,:,:,:,:,:)
    complex(8),allocatable             :: G_in(:,:,:),G_out(:,:,:),analytic_rot_(:,:)
    integer                            :: ilat,jlat,io,jo
    integer                            :: ispin,jspin
    integer                            :: iorb,jorb
    real(8)                            :: wm(Lmats),wr(Lreal),dw
    !
    write(*,*) "A(w) rotation"
    !
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    if(allocated(Gsowr))deallocate(Gsowr);allocate(Gsowr(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gsowr=zero
    if(allocated(Gsoiw))deallocate(Gsoiw);allocate(Gsoiw(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gsoiw=zero
    if(allocated( G_in))deallocate( G_in);allocate( G_in(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal));G_in=zero
    if(allocated(G_out))deallocate(G_out);allocate(G_out(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal));G_out=zero
    !
    call  ed_get_gij_lattice(Hk,Wtk,Gsoiw,Gsowr,Smats,Sreal,iprint=0)
    !
    if(mpiID==0)then

    do ilat=1,Nlat
    do jlat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                   G_in(io,jo,:)=Gsowr(ilat,jlat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    enddo
    open(unit=106,file='sum_w_Gloc.dat',status='unknown',action='write',position='rewind')
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      write(106,*) io,jo,"---",ilat,jlat,ispin,jspin,iorb,jorb,sum(abs(G_in(io,jo,:)))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    do i=1,Lreal
       G_out(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(G_in(:,:,i),impHloc_rot))
    enddo
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   call splot("Glocrot_H_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_realw.ed",wr,-dimag(G_out(io,jo,:))/pi,dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    enddo
    open(unit=106,file='sum_w_Gloc_rot_Hloc.dat',status='unknown',action='write',position='rewind')
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      write(106,*) io,jo,"---",ilat,jlat,ispin,jspin,iorb,jorb,sum(abs(G_out(io,jo,:)))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    if(allocated(analytic_rot_))deallocate(analytic_rot_)
    if(allocated(analytic_rot))deallocate(analytic_rot)
    allocate(analytic_rot_(Nspin*Norb,Nspin*Norb))
    allocate(analytic_rot(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    !
    analytic_rot=zero
    analytic_rot_=zero
    !J=1/2 jz=-1/2
    analytic_rot_(1,1)=-Xi
    analytic_rot_(3,1)=-1.0d0
    analytic_rot_(6,1)=+Xi
    analytic_rot_(:,1)=analytic_rot_(:,1)/sqrt(3.)
    !J=1/2 jz=+1/2
    analytic_rot_(2,2)=-Xi
    analytic_rot_(4,2)=+1.0d0
    analytic_rot_(5,2)=-Xi
    analytic_rot_(:,2)=analytic_rot_(:,2)/sqrt(3.)
    !J=3/2 jz=-3/2
    analytic_rot_(2,3)=-Xi
    analytic_rot_(4,3)=+1.0d0
    analytic_rot_(5,3)=+2.0d0*Xi
    analytic_rot_(:,3)=analytic_rot_(:,3)/sqrt(6.)
    !J=3/2 jz=-1/2
    analytic_rot_(1,4)=+Xi
    analytic_rot_(3,4)=-1.0d0
    analytic_rot_(:,4)=analytic_rot_(:,4)/sqrt(2.)
    !J=3/2 jz=+1/2
    analytic_rot_(2,5)=-Xi 
    analytic_rot_(4,5)=-1.0d0
    analytic_rot_(:,5)=analytic_rot_(:,5)/sqrt(2.)
    !J=3/2 jz=+3/2
    analytic_rot_(1,6)=+Xi
    analytic_rot_(3,6)=+1.0d0
    analytic_rot_(6,6)=+2.0d0*Xi
    analytic_rot_(:,6)=analytic_rot_(:,6)/sqrt(6.)
    !
    analytic_rot_=reshape_Z_to_A1(analytic_rot_)
    do ilat=1,Nlat
       analytic_rot(1+(Nspin*Norb)*(ilat-1):ilat*Nspin*Norb,1+(Nspin*Norb)*(ilat-1):ilat*Nspin*Norb)=analytic_rot_(:,:)
    enddo
    !
    G_out=zero
    do i=1,Lreal
       G_out(:,:,i)=matmul(transpose(conjg(analytic_rot)),matmul(G_in(:,:,i),analytic_rot))
    enddo
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   call splot("Glocrot_a_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_realw.ed",wr,-dimag(G_out(io,jo,:))/pi,dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    enddo
    open(unit=106,file='sum_w_Gloc_rot_anlytc.dat',status='unknown',action='write',position='rewind')
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      write(106,*) io,jo,"---",ilat,jlat,ispin,jspin,iorb,jorb,sum(abs(G_out(io,jo,:)))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    endif
  end subroutine rotate_Gloc


  !---------------------------------------------------------------------
  !PURPOSE: 
  !---------------------------------------------------------------------
  subroutine Quantum_operator()
    implicit none
    complex(8),allocatable             :: Gso(:,:,:,:,:,:)
    complex(8),allocatable             :: Stot(:,:,:,:),Ltot(:,:,:,:)
    complex(8),allocatable             :: LdotS(:)
    complex(8)                         :: Sx,Lx,Sy,Ly,Sz,Lz,jz
    integer                            :: ilat,io,jo
    integer                            :: ispin,jspin
    integer                            :: iorb,jorb
    real(8)                            :: wm(Lmats),wr(Lreal),dw
    real(8)                            :: site_mag(Nlat,Norb)
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    allocate( Gso(Nlat,Nspin,Nspin,Norb,Norb,Lmats)); Gso=zero
    call ed_get_gimp_matsubara_lattice(Gso,Nlat)
    !
    !##############################################################
    !
    !                              S
    !
    !##############################################################
    !
    write(*,*) "Computing total Spin operator per site per orbital"
    write(*,*) "Lmats used:",Lmats
    allocate(Stot(Nlat,3,Norb,Norb));Stot=zero
    !
    do ilat=1,Nlat
       do iorb=1,Norb
          do jorb=1,Norb
             !Sx
             Stot(ilat,1,iorb,jorb)=sum(    (Gso(ilat,1,2,iorb,jorb,:)+Gso(ilat,2,1,iorb,jorb,:) ))/beta
             !Sy
             Stot(ilat,2,iorb,jorb)=sum( xi*(Gso(ilat,2,1,iorb,jorb,:)-Gso(ilat,1,2,iorb,jorb,:) ))/beta
             !Sz
             Stot(ilat,3,iorb,jorb)=sum(    (Gso(ilat,1,1,iorb,jorb,:)-Gso(ilat,2,2,iorb,jorb,:) ))/beta
          enddo
       enddo
    enddo
    !
    site_mag=0.d0
    site_mag=ed_get_mag_lattice(Nlat)
    !
    open(unit=105,file='Stot_per_site.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(105,'(a100)') "#diagonal site, diagonal orbital"
    write(105,'(a8,30a20)') "#site","Re{Sx}_1","Re{Sx}_2","Re{Sx}_3" &
                                   ,"Re{Sy}_1","Re{Sy}_2","Re{Sy}_3" &
                                   ,"Re{Sz}_1","Re{Sz}_2","Re{Sz}_3" &
                                   ,"Im{Sx}_1","Im{Sx}_2","Im{Sx}_3" &
                                   ,"Im{Sy}_1","Im{Sy}_2","Im{Sy}_3" &
                                   ,"Im{Sz}_1","Im{Sz}_2","Im{Sz}_3" &
                                   ,"mag" 
    do ilat=1,Nlat
       write(105,'(I8,30F20.12)') ilat,  real(Stot(ilat,1,1,1)), real(Stot(ilat,1,2,2)), real(Stot(ilat,1,3,3)) &
                                      ,  real(Stot(ilat,2,1,1)), real(Stot(ilat,2,2,2)), real(Stot(ilat,2,3,3)) &
                                      ,  real(Stot(ilat,3,1,1)), real(Stot(ilat,3,2,2)), real(Stot(ilat,3,3,3)) &
                                      , aimag(Stot(ilat,1,1,1)),aimag(Stot(ilat,1,2,2)),aimag(Stot(ilat,1,3,3)) &       
                                      , aimag(Stot(ilat,2,1,1)),aimag(Stot(ilat,2,2,2)),aimag(Stot(ilat,2,3,3)) &
                                      , aimag(Stot(ilat,3,1,1)),aimag(Stot(ilat,3,2,2)),aimag(Stot(ilat,3,3,3)) &
                                      , site_mag(ilat,1)/2.,site_mag(ilat,2)/2.,site_mag(ilat,3)/2.
    enddo
    write(105,*)
    write(105,*)
    write(105,'(a100)') "#diagonal site, inter-orbital"
    do ilat=1,Nlat
       write(105,'(a8,I3)') "#site:",ilat
       write(105,'(30a20)') "#Sx(orb_1)","Sx(orb_2)","Sx(orb_3)","Sy(orb_1)","Sy(orb_2)","Sy(orb_3)","Sz(orb_1)","Sz(orb_2)","Sz(orb_3)"
       do iorb=1,Norb
          write(105,'(30F20.12)') (real(Stot(ilat,1,iorb,jorb)),jorb=1,Norb) &
                                 ,(real(Stot(ilat,2,iorb,jorb)),jorb=1,Norb) &
                                 ,(real(Stot(ilat,3,iorb,jorb)),jorb=1,Norb)
       enddo
       write(105,*)
       do iorb=1,Norb
          write(105,'(30F20.12)') (aimag(Stot(ilat,1,iorb,jorb)),jorb=1,Norb) &
                                 ,(aimag(Stot(ilat,2,iorb,jorb)),jorb=1,Norb) &
                                 ,(aimag(Stot(ilat,3,iorb,jorb)),jorb=1,Norb)
       enddo
    enddo
    close(105)
    !
    !##############################################################
    !
    !                              L
    !
    !##############################################################
    !
    write(*,*) "Computing total Orbital operator per site per spin"
    write(*,*) "Lmats used:",Lmats
    allocate(Ltot(Nlat,3,Nspin,Nspin));Ltot=zero
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             !Lx
             Ltot(ilat,1,ispin,jspin)=sum(  xi*(Gso(ilat,ispin,jspin,1,3,:)-Gso(ilat,ispin,jspin,3,1,:))  )/beta
             !Ly
             Ltot(ilat,2,ispin,jspin)=sum(  xi*(Gso(ilat,ispin,jspin,3,2,:)-Gso(ilat,ispin,jspin,2,3,:))  )/beta
             !Lz
             Ltot(ilat,3,ispin,jspin)=sum(  xi*(Gso(ilat,ispin,jspin,1,2,:)-Gso(ilat,ispin,jspin,2,1,:))  )/beta
          enddo
       enddo
    enddo
    !
    open(unit=106,file='Ltot_per_site.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(106,'(a100)') "#diagonal site, diagonal spin"
    write(106,'(a8,30a20)') "#site","Re{Lx}_1","Re{Lx}_2" &
                                   ,"Re{Ly}_1","Re{Ly}_2" &
                                   ,"Re{Lz}_1","Re{Lz}_2" &
                                   ,"Im{Lx}_1","Im{Lx}_2" &
                                   ,"Im{Ly}_1","Im{Ly}_2" &
                                   ,"Im{Lz}_1","Im{Lz}_2"
    do ilat=1,Nlat
       write(106,'(I8,30F20.12)') ilat,  real(Ltot(ilat,1,1,1)), real(Ltot(ilat,1,2,2)) &
                                      ,  real(Ltot(ilat,2,1,1)), real(Ltot(ilat,2,2,2)) &
                                      ,  real(Ltot(ilat,3,1,1)), real(Ltot(ilat,3,2,2)) &
                                      , aimag(Ltot(ilat,1,1,1)),aimag(Ltot(ilat,1,2,2)) &       
                                      , aimag(Ltot(ilat,2,1,1)),aimag(Ltot(ilat,2,2,2)) &
                                      , aimag(Ltot(ilat,3,1,1)),aimag(Ltot(ilat,3,2,2))
    enddo
    write(106,*)
    write(106,*)
    write(106,'(a100)') "#diagonal site, inter-orbital"
    do ilat=1,Nlat
       write(106,'(a8,I3)') "#site:",ilat
       write(106,'(30a20)') "#Lx(spin_1)","Lx(spin_2)","Ly(spin_1)","Ly(spin_2)","Lz(spin_1)","Lz(spin_2)"
       do ispin=1,Nspin
          write(106,'(30F20.12)') (real(Ltot(ilat,1,ispin,jspin)),jspin=1,Nspin) &
                                 ,(real(Ltot(ilat,2,ispin,jspin)),jspin=1,Nspin) &
                                 ,(real(Ltot(ilat,3,ispin,jspin)),jspin=1,Nspin)
       enddo
       write(106,*)
       do ispin=1,Nspin
          write(106,'(30F20.12)') (aimag(Ltot(ilat,1,ispin,jspin)),jspin=1,Nspin) &
                                 ,(aimag(Ltot(ilat,2,ispin,jspin)),jspin=1,Nspin) &
                                 ,(aimag(Ltot(ilat,3,ispin,jspin)),jspin=1,Nspin)
       enddo
    enddo
    close(106)
    !
    !##############################################################
    !
    !                              L.dot.S
    !
    !##############################################################
    !
    write(*,*) "Computing total L dot S operator per site"
    write(*,*) "Lmats used:",Lmats
    allocate(LdotS(Nlat));LdotS=zero
    !
    do ilat=1,Nlat
       LdotS(ilat)=sum(       +xi*Gso(ilat,1,1,1,2,:) &
                              +xi*Gso(ilat,1,2,1,3,:) &  
                              -xi*Gso(ilat,2,2,1,2,:) &  
                              +xi*Gso(ilat,2,1,1,3,:) &  
                              -xi*Gso(ilat,1,1,2,1,:) &  
                              -   Gso(ilat,1,2,2,3,:) &  
                              +xi*Gso(ilat,2,2,2,1,:) &  
                              +   Gso(ilat,2,1,2,3,:) &  
                              -xi*Gso(ilat,1,2,3,1,:) &  
                              +   Gso(ilat,1,2,3,2,:) &  
                              -xi*Gso(ilat,2,1,3,1,:) &  
                              -   Gso(ilat,2,1,3,2,:) &  
                      )/beta
       LdotS(ilat)=LdotS(ilat)/2.d0
    enddo
    !
    open(unit=107,file='Jz_per_site.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(107,'(a8,30a20)') "#site","Re{Sz}_11","Re{Sz}_22","Re{Sz}_33","Im{Sz}_11","Im{Sz}_22","Im{Sz}_33","Re{Tr[Sz]}","Im{Tr[Sz]}" &
                                                           ,"Re{Lz}_uu","Im{Lz}_uu","Re{Lz}_dd","Im{Lz}_dd","Re{Tr[Lz]}","Im{Tr[Lz]}" &
                                                           ,"Re{jz}","Im{jz}","Re{L.S}","Im{L.S}"
    do ilat=1,Nlat
       Sx=trace(Stot(ilat,1,:,:));Sy=trace(Stot(ilat,2,:,:));Sz=trace(Stot(ilat,3,:,:))
       Lx=trace(Ltot(ilat,1,:,:));Ly=trace(Ltot(ilat,2,:,:));Lz=trace(Ltot(ilat,3,:,:))
       jz=Sz+Lz
       write(107,'(I8,30F20.12)') ilat,  real(Stot(ilat,3,1,1)), real(Stot(ilat,3,2,2)), real(Stot(ilat,3,3,3)) &
                                      , aimag(Stot(ilat,3,1,1)),aimag(Stot(ilat,3,2,2)),aimag(Stot(ilat,3,3,3)),real(Sz),aimag(Sz) &
                                      ,  real(Ltot(ilat,3,1,1)), real(Ltot(ilat,3,2,2)) &
                                      , aimag(Ltot(ilat,3,1,1)),aimag(Ltot(ilat,3,2,2)),real(Lz),aimag(Lz) &
                                      ,  real(jz),aimag(jz),real(LdotS(ilat)),aimag(LdotS(ilat))
    enddo
    close(107)
    !
    deallocate(Ltot,Stot,LdotS)
    !
  end subroutine Quantum_operator


  !---------------------------------------------------------------------
  !PURPOSE: old useless
  !---------------------------------------------------------------------
  function Jz_builder(alf) result(A)
    real(8)                                               :: alf
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: A
    complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: B
    integer                                               :: i,j,iorb,jorb,ispin,jspin,io,jo,ilat,jlat
    A = zero;B = zero
    do i=1,3
       B(i,i)=cmplx(1.d0,0.d0)
       B(i+3,i+3)=cmplx(-1.d0,0.d0)
    enddo
    B(1,2)=cmplx(0.d0,-1.d0)
    B(2,1)=cmplx(0.d0,+1.d0)
    B(4,5)=cmplx(0.d0,-1.d0)
    B(5,4)=cmplx(0.d0,+1.d0)
    do ilat=1,Nlat
       do jlat=1,Nlat
          if (ilat.eq.jlat) A(1+(ilat-1)*Nspin*Norb:Nspin*Norb*ilat,1+(jlat-1)*Nspin*Norb:Nspin*Norb*jlat)=B
          if ((ilat.eq.1).and.(jlat.eq.2)) A(1+(ilat-1)*Nspin*Norb:Nspin*Norb*ilat,1+(jlat-1)*Nspin*Norb:Nspin*Norb*jlat)=B*alf
          if ((ilat.eq.3).and.(jlat.eq.4)) A(1+(ilat-1)*Nspin*Norb:Nspin*Norb*ilat,1+(jlat-1)*Nspin*Norb:Nspin*Norb*jlat)=B*alf
          if ((ilat.eq.5).and.(jlat.eq.6)) A(1+(ilat-1)*Nspin*Norb:Nspin*Norb*ilat,1+(jlat-1)*Nspin*Norb:Nspin*Norb*jlat)=B*alf
          if ((ilat.eq.7).and.(jlat.eq.8)) A(1+(ilat-1)*Nspin*Norb:Nspin*Norb*ilat,1+(jlat-1)*Nspin*Norb:Nspin*Norb*jlat)=B*alf
       enddo
    enddo
    do i=1,Nspin*Norb*Nlat
       do j=i+1,Nspin*Norb*Nlat
          A(j,i)=conjg(A(i,j))
       enddo
    enddo
  end function Jz_builder


  !_______________________________________________________________________
  !                            reshape functions
  !_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: reshape functions
  !  Z  = [Nspin,Nspin]*Norb
  !  A1 = [Norb*Norb]*Nspin
  !  A2 = [Nspin,Nspin,Norb,Norb]
  !---------------------------------------------------------------------
  function reshape_Z_to_A1(fg) result(g)
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
  end function reshape_Z_to_A1

  function reshape_A1_to_Z(fg) result(g)
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
                !switchHloc
                g(io1,jo1)  = fg(io2,jo2)
                !
             enddo
          enddo
       enddo
    enddo
  end function reshape_A1_to_Z

  function reshape_A1_to_A2(fg) result(g)
    complex(8),dimension(Nso,Nso)                   :: fg
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
  end function reshape_A1_to_A2

  function reshape_A1_to_A2_L(fg) result(g)
    complex(8),dimension(Nlat*Nso,Nlat*Nso)          :: fg
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: g
    integer                                          :: i,j,iorb,jorb,ispin,jspin,io,jo,ilat
    g = zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   g(ilat,ispin,jspin,iorb,jorb)  = fg(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function reshape_A1_to_A2_L

  function reshape_A1_to_A2_L2(fg) result(g)
    complex(8),dimension(Nlat*Nso,Nlat*Nso)               :: fg
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: g
    integer                                               :: i,j,iorb,jorb,ispin,jspin,io,jo,ilat,jlat
    g = zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      g(ilat,jlat,ispin,jspin,iorb,jorb)  = fg(io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function reshape_A1_to_A2_L2

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
       write(*,*) "Matrices not equal cannot perform inversion test"
       stop
    endif
    dime=maxval(shape(A))
    if (abs(float(dime)-real(sum(matmul(A,B)))).gt.tol) write(*,'(A30)') "inversion test fail"
  end subroutine inversion_test


end program ed_SIO

