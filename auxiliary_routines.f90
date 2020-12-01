!This is a collection of auxiliary routines which are often used in 
!the drivers, but do not find space in the ED code:





!##################################################################
!##################################################################
!                   TRASFORMATION ROUTINES
!##################################################################
!##################################################################

! !+-----------------------------------------------------------------------------+!
! !PURPOSE:
! ! bcast the Blocks vector [Nlat][Nspin*Norb][Nspin*Norb]
! ! into a large matrix [Nlat*Nspin*Norb][Nlat*Nspin*Norb]
! !+-----------------------------------------------------------------------------+!
! function blocks_to_matrix(Vblocks) result(Matrix)
!   complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb)      :: Vblocks
!   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
!   integer                                               :: i,j,ip
!   Matrix=zero
!   do ip=1,Nlat
!      i = 1 + (ip-1)*Nspin*Norb
!      j = ip*Nspin*Norb
!      Matrix(i:j,i:j) =  Vblocks(ip,:,:)
!   enddo
! end function blocks_to_matrix


!+-----------------------------------------------------------------------------+!
!PURPOSE:
! bcast the diagonal part of a large matrix [Nlat*Nspin*Norb][Nlat*Nspin*Norb]
! into Blocks vector [Nlat][Nspin*Norb][Nspin*Norb]
!+-----------------------------------------------------------------------------+!
function matrix_to_blocks(Matrix) result(Vblocks)
  complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
  complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb)      :: Vblocks
  integer                                               :: i,j,ip
  Vblocks=zero
  do ip=1,Nlat
     i = 1 + (ip-1)*Nspin*Norb
     j = ip*Nspin*Norb
     Vblocks(ip,:,:) = Matrix(i:j,i:j)
  enddo
end function matrix_to_blocks


!+-----------------------------------------------------------------------------+!
!PURPOSE:
! select a single block of the diagonal from a large matrix.
!   + _Nlso = matrix has dimensions Nlso*Nlso
!   + _NNN  = matrix has dimensions Nlat,Nspin,Nspin,Norb,Norb
!+-----------------------------------------------------------------------------+!
function select_block_Nlso(ip,Matrix) result(Vblock)
  integer                                               :: ip
  complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
  complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: Vblock
  integer                                               :: i,j
  Vblock=zero
  i = 1+(ip-1)*Nspin*Norb
  j =       ip*Nspin*Norb
  Vblock(:,:) = Matrix(i:j,i:j)
end function select_block_nlso
!
function select_block_nnn(ip,Matrix) result(Vblock)
  integer                                          :: ip
  complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
  complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
  integer                                          :: is,js,ispin,jspin,iorb,jorb
  Vblock=zero
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              is = iorb + (ispin-1)*Norb !spin-orbit stride
              js = jorb + (jspin-1)*Norb !spin-orbit stride
              Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
           enddo
        enddo
     enddo
  enddo
end function select_block_nnn
















!+-----------------------------------------------------------------------------+!
!PURPOSE:  
! find all inequivalent sites with respect the user defined symmetry operations 
! Build and check maps from the full(independent) lattice to the independent 
! (full) lattice                        !
!+-----------------------------------------------------------------------------+!
subroutine get_independent_sites(symmetry_operations,Nsymm,Nindep)
  integer,intent(in)                 :: Nsymm
  integer,intent(inout)              :: Nindep
  integer                            :: i,unit,isymm
  integer,dimension(Nlat)            :: tmp_search
  integer,dimension(Nlat,Nsymm)      :: tmp_map
  integer                            :: i_ind,check_maps
  character(len=5)                   :: tmp_suffix
  ! integer,dimension(:),allocatable   :: map_lat2ind
  ! integer,dimension(:,:),allocatable :: map_ind2lat
  !integer,allocatable,dimension(:)   :: indep_list
  interface
     function symmetry_operations(site_in) result(sites_out)
       implicit none
       integer                  :: site_in
       integer,allocatable      :: sites_out(:)
     end function symmetry_operations
  end interface
  !+- search and get number of independent sites -+!
  tmp_search=0
  i_ind=0
  do i=1,Nlat
     tmp_map(i,:)=symmetry_operations(i)
     if(tmp_search(i).ge.0) then
        i_ind=i_ind+1
        tmp_search(i)=i
        do isymm=1,Nsymm
           tmp_search(tmp_map(i,isymm))=-1
        end do
     end if
  end do
  write(LOGfile,*) Nlat
  !
  Nindep=i_ind
  ! (remember: each site is connected with Nsymm sites (+ 1 = identity)) !
  allocate(indep_list(Nindep),map_lat2ind(Nlat),map_ind2lat(Nindep,Nsymm+1))  
  !
  !+- get list of independent sites -+!
  i_ind=0
  unit=free_unit()    
  open(unit,file='independent_sites.lattice')
  do i=1,Nlat
     if(tmp_search(i).ge.0) then
        i_ind=i_ind+1
        indep_list(i_ind) = tmp_search(i)
        write(unit,*) dble(icol(indep_list(i_ind))),dble(irow(indep_list(i_ind))),dble(indep_list(i_ind)) 
     end if
  end do
  close(unit)
  !+-  build maps -+!
  !
  write(LOGfile,*) "NINDEP",Nindep
  write(LOGfile,*) indep_list
  do i_ind=1,Nindep
     map_lat2ind(indep_list(i_ind))=i_ind
     do isymm=1,Nsymm
        map_lat2ind(tmp_map(indep_list(i_ind),isymm))=i_ind
     end do
  end do
  ! 
  do i_ind=1,Nindep
     unit=free_unit()
     !write(tmp_suffix,'(I4.4)') i_ind
     ed_file_suffix="_site"//reg(txtfy(i_ind,Npad=4))!trim(tmp_suffix)
     open(unit,file='equivalents'//trim(tmp_suffix)//'.lattice')
     map_ind2lat(i_ind,1) = indep_list(i_ind)
     write(unit,*) icol(indep_list(i_ind)),irow(indep_list(i_ind))
     do isymm=1,Nsymm
        map_ind2lat(i_ind,isymm+1) = tmp_map(indep_list(i_ind),isymm)
        write(unit,*) icol(tmp_map(indep_list(i_ind),isymm)),irow(tmp_map(indep_list(i_ind),isymm))
     end do
     close(unit)
  end do
  !+- check maps +-!
  do i_ind=1,Nindep
     do isymm=1,Nsymm+1
        check_maps=map_ind2lat(i_ind,isymm)
        if(i_ind /= map_lat2ind(check_maps)) stop "WRONG MAPS"
     end do
  end do
  ed_file_suffix=""
  ! allocate(Ineq_sites_list(Nindep))
  ! Ineq_sites_list = indep_list
end subroutine get_independent_sites



  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  
  ! get the block diagonal local part of a given Hamiltonian H(k_\perp;Ri,Rj)
  !+-----------------------------------------------------------------------------+!
  function extract_Hloc_Nlso(Hk,Nlat,Nspin,Norb) result(Hloc)
    complex(8),dimension(:,:,:)                 :: Hk
    integer                                     :: Nlat,Nspin,Norb
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Hloc
    !
    integer                                     :: iorb,ispin,ilat,is
    integer                                     :: jorb,jspin,js
    Hloc = zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hloc(is,js) = sum(Hk(is,js,:))/size(Hk,3)
                enddo
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
  end function extract_Hloc_Nlso

  function extract_Hloc_Nso(Hk,Nspin,Norb) result(Hloc)
    complex(8),dimension(:,:,:)                 :: Hk
    integer                                     :: Nspin,Norb
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Hloc
    !
    integer                                     :: i,iorb,ispin,is
    integer                                     :: j,jorb,jspin,js
    Hloc = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i = stride_index([iorb,ispin],[Norb,Nspin])
                j = stride_index([jorb,jspin],[Norb,Nspin])
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Hloc(is,js) = sum(Hk(is,js,:))/size(Hk,3)
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
  end function extract_Hloc_Nso

