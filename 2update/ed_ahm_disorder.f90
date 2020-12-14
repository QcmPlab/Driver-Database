! disorder realization depends on the parameter int idum:
! so that different realizations (statistics) are performed 
! calling this program many times providing a *different* seed 
! +IDUM. 
program ed_ahm_disorder
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  complex(8),allocatable,dimension(:,:,:)     :: Smats,Sreal !self_energies
  complex(8),allocatable,dimension(:,:,:)     :: Gmats,Greal !local green's functions
  complex(8),allocatable,dimension(:,:,:)     :: Delta,errDelta
  real(8),allocatable,dimension(:)            :: erandom
  real(8),allocatable,dimension(:,:)          :: bath,bath_prev,errBath
  logical                                     :: converged,phsym,bool,gfpoles,getkin,bdgflag
  real(8)                                     :: wmixing,Wdis,ts,Eout(2)
  integer                                     :: i,is,iloop,nrandom,idum
  integer                                     :: Nb,Lf
  real(8),dimension(:),allocatable            :: wm,wr
  real(8),dimension(:,:),allocatable          :: nii,dii,pii,eii ![Nlat][Norb]/[4]
  complex(8),dimension(:,:,:),allocatable     :: Hk              ![Nlat][Nlat][Nk=1]
  complex(8),dimension(:,:,:,:,:),allocatable :: Hloc
  character(len=50)                           :: finput

  ! START MPI !
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)


  ! READ INPUT FILES !
  call parse_cmd_variable(finput,"FINPUT",default="inputDSC.conf")
  call parse_input_variable(ts,"TS",finput,default=0.5d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(Wdis,"WDIS",finput,default=0.d0)
  call parse_input_variable(idum,"IDUM",finput,default=1234567)
  call parse_input_variable(phsym,"PHSYM",finput,default=.false.)
  call parse_input_variable(gfpoles,"GFPOLES",finput,default=.false.)
  call parse_input_variable(getkin,"GETKIN",finput,default=.false.)
  call parse_input_variable(bdgflag,"BDGFLAG",finput,default=.false.,comment="Perform Bogoliubov-deGennes calculations rather than DMFT")
  call ed_read_input(trim(finput))
  if(gfpoles.AND.ed_verbose>4)ed_verbose=4

  !CHECK SETUP:
  if(Nspin/=1.OR.Norb/=1)stop "ed_ahm_disorder error: Nspin != 1 OR Norb != 1"


  ! SET THE COMPRESSION THRESHOLD TO 1Mb (1024Kb)
  call set_store_size(1024)


  ! GET THE ACTUAL NUMBER OF LATTICE SITES !
  Nlat = Nside*Nside

  ! ALLOCATE ALL THE REQUIRED VARIABLES !
  allocate(nii(Nlat,Norb))
  allocate(dii(Nlat,Norb))
  allocate(pii(Nlat,Norb))
  allocate(eii(Nlat,4))
  allocate(Smats(2,Nlat,Lmats))
  allocate(Sreal(2,Nlat,Lreal))
  allocate(Gmats(2,Nlat,Lmats))
  allocate(Greal(2,Nlat,Lreal))
  allocate(Delta(2,Nlat,Lmats))


  ! ALLOCATE MATSUBARA AND REAL FREQUENCIES !
  allocate(wm(Lmats),wr(Lreal))
  wini  =wini-Wdis
  wfin  =wfin+Wdis
  wr = linspace(wini,wfin,Lreal)
  wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)


  ! GET RANDOM ENERGIES
  allocate(erandom(Nlat))
  call random_seed(size=nrandom)  
  call random_seed(put=(/(idum,i=1,nrandom)/))
  call random_number(erandom)
  erandom=(2.d0*erandom-1.d0)*Wdis/2.d0
  inquire(file='erandom.restart',exist=bool)
  if(bool)then
     Lf = file_length('erandom.restart')
     if(Lf/=Nlat)stop "ed_ahm_disorder error: found erandom.restart with length different from Nlat"
     call read_data('erandom.restart',erandom)
  endif
  call store_data("erandom.ed",erandom)


  ! ALLOCATE A BATH FOR EACH IMPURITY 
  Nb=get_bath_size()
  allocate(bath(Nlat,Nb), bath_prev(Nlat,Nb))



  ! GET THE TIGHT BINDING HAMILTONIAN FOR THE SQUARE LATTICE 
  allocate(Hk(Nlat,Nlat,1))
  Hk(:,:,1) = one*Htb_square_lattice(Nrow=Nside,Ncol=Nside,ts=ts)
  Hk(:,:,1) = Hk(:,:,1) + diag(Erandom)


  ! SET THE LOCAL PART WITH DISORDER
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  Hloc(1:Nlat,1,1,1,1) = one*Erandom(:)


  ! INITIALIZE DMFT SOLVER
  if(bdgflag)then
     nii=0.5d0
     pii=deltasc
     call read_hfb_solution()
  else
     call ed_init_solver_lattice(bath)
  endif

  if(getkin)then
     if(bdgflag)then
        call read_hfb_solution()
        if(mpiID==0)Eout = Bdg_Kinetic(nii(:,1),pii(:,1),dreal(Hk))
     else
        call read_sigma_matsubara()
        Eout = ed_kinetic_energy_lattice(Hk,[1d0],Smats(1,:,:),Smats(2,:,:))
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
     call MPI_FINALIZE(mpiERR)
     stop
  endif

  ! DMFT LOOP
  ! everywhere iprint=0: printing is postponed to the end
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1

     if(bdgflag)then
        if(mpiID==0)call start_loop(iloop,nloop,"BdG-loop",unit=LOGfile)
        if(mpiID==0)call bdg_solve(nii(:,1),pii(:,1),dreal(Hk))
        call MPI_BCAST(nii,size(nii),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
        call MPI_BCAST(pii,size(pii),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
        dii=nii*nii
        eii = uloc(1)*dii
     else
        if(mpiID==0)call start_loop(iloop,nloop,"DMFT-loop",unit=LOGfile)
        bath_prev = bath
        call ed_solve_lattice(bath,Hloc,iprint=0)
        if(gfpoles)then
           call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
           call MPI_FINALIZE(mpiERR)
           stop
        endif
        nii = ed_get_dens_lattice(Nlat)
        dii = ed_get_docc_lattice(Nlat)
        pii = ed_get_phisc_lattice(Nlat)
        eii = ed_get_eimp_lattice(Nlat)
        call ed_get_sigma_matsubara_lattice(Smats(1,:,:),Nlat)
        call ed_get_self_matsubara_lattice(Smats(2,:,:),Nlat)
        call ed_get_sigma_real_lattice(Sreal(1,:,:),Nlat)
        call ed_get_self_real_lattice(Sreal(2,:,:),Nlat)
        call ed_get_gloc_lattice(Hk,[1d0],Gmats,Greal,Smats,Sreal,iprint=0)
        call ed_get_weiss_lattice(Gmats,Smats,Delta,Hloc,iprint=0)
        call ed_chi2_fitgf_lattice(bath,Delta,Hloc)
        bath = wmixing*bath + (1.d0-wmixing)*bath_prev
        if(phsym)then
           do i=1,Nlat
              call ph_symmetrize_bath(bath(i,:))
           enddo
        endif
     endif

     if(mpiID==0)converged = check_convergence_local(pii,dmft_error,Nsuccess,nloop,id=0,file="error.err")
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     if(bdgflag.AND.converged)then
        do i=1,Lmats
           Smats(1,:,i) =  Uloc(1)*(nii(:,1)-0.5d0)
           Smats(2,:,i) =  Uloc(1)*pii(:,1)
        enddo
        do i=1,Lreal
           Sreal(1,:,i) =  Uloc(1)*(nii(:,1)-0.5d0)
           Sreal(2,:,i) =  Uloc(1)*pii(:,1)
        enddo
        call ed_get_gloc_lattice(Hk,[1d0],Gmats,Greal,Smats,Sreal,iprint=0)
     endif
     call print_sc_out(converged)
     if(mpiID==0)call end_loop()
  enddo
  !+-------------------------------------+!

  if(bdgflag)then
     if(mpiID==0)Eout = Bdg_Kinetic(nii(:,1),pii(:,1),dreal(Hk))
  else
     Eout = ed_kinetic_energy_lattice(Hk,[1d0],Smats(1,:,:),Smats(2,:,:))
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)


contains





  ! !******************************************************************
  ! !******************************************************************




  subroutine read_sigma_matsubara()
    if(mpiID==0)then
       call read_data("LSigma_iw.ed",Smats(1,1:Nlat,1:Lmats),wm(1:Lmats))
       call read_data("LSelf_iw.ed",Smats(2,1:Nlat,1:Lmats),wm(1:Lmats))
    endif
    call MPI_BCAST(Smats,size(Smats),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  end subroutine read_sigma_matsubara

  subroutine read_hfb_solution()
    if(mpiID==0)then
       call read_data("nVSisite.ed",nii(:,1))
       call read_data("phiVSisite.ed",pii(:,1))
    endif
    call MPI_BCAST(nii,size(nii),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BCAST(pii,size(pii),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine read_hfb_solution




  ! !******************************************************************
  ! !******************************************************************





  subroutine BdG_Solve(nii,pii,H)
    real(8),dimension(Nlat),intent(inout)   :: nii
    real(8),dimension(Nlat),intent(inout)   :: pii
    real(8),dimension(Nlat,Nlat),intent(in) :: H
    real(8),dimension(Nlat)                 :: Sigma_HFB
    real(8),dimension(Nlat)                 :: Self_HFB
    real(8),dimension(2*Nlat,2*Nlat)        :: Hnambu
    real(8),dimension(2*Nlat)               :: Enambu
    real(8),dimension(2*Nlat)               :: RhoDiag  !diagonal density matrix 
    real(8),dimension(2*Nlat,2*Nlat)        :: RhoNambu !Nambu density matrix
    integer                                 :: ilat
    write(*,*)"Solving B-dG problem:"
    Sigma_HFB(:) =  Uloc(1)*(nii-0.5d0)
    Self_HFB(:)  =  Uloc(1)*pii(:)
    Hnambu(1:Nlat,1:Nlat)               =  H + diag(Sigma_HFB)
    Hnambu(1:Nlat,Nlat+1:2*Nlat)        =    + diag(Self_HFB)
    Hnambu(Nlat+1:2*Nlat,1:Nlat)        =    + diag(Self_HFB)
    Hnambu(Nlat+1:2*Nlat,Nlat+1:2*Nlat) = -H - diag(Sigma_HFB)
    call eigh(HNambu,ENambu)
    RhoDiag  = fermi(ENambu,beta)
    RhoNambu = matmul(HNambu, matmul(diag(RhoDiag),transpose(HNambu)) ) 
    forall(ilat=1:Nlat)
       nii(ilat) = RhoNambu(ilat,ilat)
       pii(ilat) = RhoNambu(ilat,ilat+Nlat)
    end forall
  end subroutine BdG_Solve


  function BdG_Kinetic(nii,pii,H) result(Ekin)
    real(8),dimension(Nlat),intent(inout)   :: nii
    real(8),dimension(Nlat),intent(inout)   :: pii
    real(8),dimension(Nlat,Nlat),intent(in) :: H
    real(8),dimension(Nlat)                 :: Sigma_HFB
    real(8),dimension(Nlat)                 :: Self_HFB
    real(8),dimension(2*Nlat,2*Nlat)        :: Hnambu,HHNambu
    real(8),dimension(2*Nlat)               :: Enambu
    real(8),dimension(2*Nlat)               :: RhoDiag  !diagonal density matrix 
    real(8),dimension(2*Nlat,2*Nlat)        :: RhoNambu !Nambu density matrix
    integer                                 :: ilat,unit
    real(8),dimension(2)                    :: Ekin
    write(*,*)"Get B-dG Kinetic Energy:"
    Ekin=0d0
    Sigma_HFB(:) =  Uloc(1)*(nii-0.5d0)
    Self_HFB(:)  =  Uloc(1)*pii
    HHNambu=zero
    HHnambu(1:Nlat,1:Nlat)               =  H
    HHnambu(Nlat+1:2*Nlat,Nlat+1:2*Nlat) = -H
    !
    Hnambu(1:Nlat,1:Nlat)               =  H + diag(Sigma_HFB)
    Hnambu(1:Nlat,Nlat+1:2*Nlat)        =    + diag(Self_HFB)
    Hnambu(Nlat+1:2*Nlat,1:Nlat)        =    + diag(Self_HFB)
    Hnambu(Nlat+1:2*Nlat,Nlat+1:2*Nlat) = -H - diag(Sigma_HFB)
    call eigh(HNambu,ENambu)
    RhoDiag   = fermi(ENambu,beta)
    RhoNambu  = matmul(HNambu, matmul(diag(RhoDiag),transpose(HNambu)) ) 
    HNambu    = matmul(HHNambu,RhoNambu)
    Ekin(1) = (3-Nspin)*trace(HNambu(1:Nlat,1:Nlat))/Nlat
    !
    unit = free_unit()
    open(unit,file="kinetic_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",reg(txtfy(1))//"<K>",reg(txtfy(2))//"<Eloc>"
    close(unit)
    !
    unit = free_unit()
    open(unit,file="kinetic_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")Ekin(1),Ekin(2)
    close(unit)
    !
  end function BdG_Kinetic






  ! !******************************************************************
  ! !******************************************************************




  subroutine print_sc_out(converged)
    integer                        :: i,j,is,row,col,unit
    real(8)                        :: nimp,phi,ccdw,docc
    real(8),dimension(Nlat)        :: cdwii,rii,sii,zii
    real(8),dimension(Nside,Nside) :: dij,nij,cij,pij,eij
    real(8),dimension(Nside)       :: grid_x,grid_y
    real(8)                        :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Nlat)      :: data_covariance
    real(8),dimension(2,2)         :: covariance_nd
    real(8),dimension(2)           :: data_mean,data_sdev,Eout
    logical                        :: converged
    complex(8),dimension(2,Lmats)  :: aGmats,aSmats
    complex(8),dimension(2,Lreal)  :: aGreal,aSreal
    character(len=50)              :: suffix,cloop,cfoo

    if(bdgflag)nii=2*nii
    if(mpiID==0)then
       suffix=".ed"
       write(cfoo,'(I4.4)')iloop
       cloop = "_loop"//trim(cfoo)//".ed"
       !Get CDW "order parameter"
       do is=1,Nlat
          row=irow(is)
          col=icol(is)
          cdwii(is) = (-1.d0)**(row+col)*(nii(is,1)-1.d0)
       enddo
       nimp = sum(nii(:,1))/dble(Nlat)
       phi  = sum(pii(:,1))/dble(Nlat)
       docc = sum(dii(:,1))/dble(Nlat)
       ccdw = sum(cdwii)/dble(Nlat)
       print*,"<nimp>  =",nimp
       print*,"<phi>   =",phi
       print*,"<docc>  =",docc
       print*,"<ccdw>  =",ccdw
       call splot("nVSiloop.ed",iloop,nimp,append=.true.)
       call splot("phiVSiloop.ed",iloop,phi,append=.true.)
       call splot("doccVSiloop.ed",iloop,docc,append=.true.)
       call splot("ccdwVSiloop.ed",iloop,ccdw,append=.true.)
       call store_data("nVSisite"//trim(suffix),nii(:,1))
       call store_data("phiVSisite"//trim(suffix),pii(:,1))
       call store_data("doccVSisite"//trim(suffix),dii(:,1))
       call store_data("epotVSisite"//trim(suffix),eii(:,1))
       call store_data("cdwVSisite"//trim(suffix),cdwii)
       !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
       if(converged)then
          unit=free_unit()
          open(unit,file='internal_energy.ed')
          write(unit,'(10(F18.10))')Eout(1),sum(eii(:,1))/Nlat,sum(eii(:,2))/Nlat,sum(eii(:,3))/Nlat,sum(eii(:,4))/Nlat,Eout(2)
          close(unit)
          open(unit,file='internal_energy_info.ed')
          write(unit,"(A,10A18)")"#","1<K>","2<Hi>","3<V=Hi-Ehf>","4<Eloc>","5<Ehf>","6<E0>"
          close(unit)
          call store_data("LG_iw"//trim(suffix),Gmats(1,1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LF_iw"//trim(suffix),Gmats(2,1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LG_realw"//trim(suffix),Greal(1,1:Nlat,1:Lreal),wr(1:Lreal))
          call store_data("LF_realw"//trim(suffix),Greal(2,1:Nlat,1:Lreal),wr(1:Lreal))
          call store_data("LSigma_iw"//trim(suffix),Smats(1,1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LSelf_iw"//trim(suffix),Smats(2,1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LSigma_realw"//trim(suffix),Sreal(1,1:Nlat,1:Lreal),wr(1:Lreal))
          call store_data("LSelf_realw"//trim(suffix),Sreal(2,1:Nlat,1:Lreal),wr(1:Lreal))
          do row=1,Nside
             grid_x(row)=row
             grid_y(row)=row
             do col=1,Nside
                i            = ij2site(row,col)
                nij(row,col) = nii(i,1)
                dij(row,col) = dii(i,1)
                pij(row,col) = pii(i,1)
                eij(row,col) = eii(i,1)
             enddo
          enddo
          call splot3d("3d_nVSij"//trim(suffix),grid_x,grid_y,nij)
          call splot3d("3d_doccVSij"//trim(suffix),grid_x,grid_y,dij)
          call splot3d("3d_phiVSij"//trim(suffix),grid_x,grid_y,pij)
          call splot3d("3d_epotVSij"//trim(suffix),grid_x,grid_y,eij)
          call get_moments(nii(:,1),mean,sdev,var,skew,kurt)
          call splot("statistics.n"//trim(suffix),mean,sdev,var,skew,kurt)
          data_mean(1)=mean
          data_sdev(1)=sdev
          !
          call get_moments(dii(:,1),mean,sdev,var,skew,kurt)
          call splot("statistics.docc"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(pii(:,1),mean,sdev,var,skew,kurt)
          call splot("statistics.phi"//trim(suffix),mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          !
          call get_moments(cdwii,mean,sdev,var,skew,kurt)
          call splot("statistics.cdwn"//trim(suffix),mean,sdev,var,skew,kurt)
          data_covariance(1,:)=nii(:,1)
          data_covariance(2,:)=pii(:,1)
          covariance_nd = get_covariance(data_covariance,data_mean)
          open(10,file="covariance_n.phi"//trim(suffix))
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)
          covariance_nd=0d0
          do i=1,2
             do j=1,2
                if(data_sdev(i)/=0d0.AND.data_sdev(j)/=0d0)then
                   covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
                endif
             enddo
          enddo
          open(10,file="correlation_n.phi"//trim(suffix))
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)
       end if
    end if

    if(bdgflag)nii=nii/2

  end subroutine print_sc_out



  ! !******************************************************************
  ! !******************************************************************



  function Htb_square_lattice(Nrow,Ncol,pbc_row,pbc_col,ts) result(H0)
    integer                                :: Nrow
    integer                                :: Ncol
    logical,optional                       :: pbc_row,pbc_col
    logical                                :: pbc_row_,pbc_col_
    real(8),optional                       :: ts
    real(8)                                :: ts_
    real(8),dimension(Nrow*Ncol,Nrow*Ncol) :: H0
    integer                                :: i,jj,row,col,link(4),j
    integer                                :: unit
    !
    pbc_row_=.true. ; if(present(pbc_row)) pbc_row_=pbc_row
    pbc_col_=.true. ; if(present(pbc_col)) pbc_col_=pbc_col
    ts_=0.5d0;if(present(ts))ts_=ts
    !
    H0 = 0.d0
    unit=free_unit()
    !+- 2D LATTICE (NROW x NCOL) -+!
    if(Nlat /= Nrow*Ncol) stop "get_lattice_hamiltonian error: Nlat != Nrow*Ncol"
    !THESE ARE STILL GLOBAL VARIABLES...
    allocate(icol(Nlat),irow(Nlat))
    allocate(ij2site(Nrow,Ncol))
    do row=0,Nrow-1
       do col=0,Ncol-1
          i=col+ 1 + row*Ncol
          !
          irow(i)=row+1
          icol(i)=col+1
          ij2site(row+1,col+1)=i
          !
          !right hop
          link(1)= i + 1     
          if((col+1)==Ncol) then
             link(1)=0  
             if(pbc_col_)link(1)=1+row*Ncol  
          end if
          !left  hop
          link(3)= i - 1    
          if((col-1)<0)     then
             link(3)=0  
             if(pbc_col_)link(3)=Ncol+row*Ncol
          end if
          !up    hop
          link(2)= i + Ncol 
          if((row+1)==Nrow) then
             link(2)=0  
             if(pbc_row_)link(2)=col+1
          end if
          !down  hop
          link(4)= i - Ncol 
          if((row-1)<0)     then
             link(4)=0  
             if(pbc_row_)link(4)=col+1+(Nrow-1)*Ncol
          end if
          !
          do jj=1,4
             if(link(jj)>0)H0(i,link(jj))=-ts_ !! ts must be negative.
          enddo
          !
       enddo
    enddo
    if(mpiID==0) then
       open(unit,file='Htb_square_lattice.ed')
       do i=1,Nlat
          write(unit,"(5000(F5.2,1x))")(H0(i,j),j=1,Nlat)
       enddo
       close(unit)
    endif
  end function Htb_square_lattice



end program ed_ahm_disorder
