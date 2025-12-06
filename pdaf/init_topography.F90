SUBROUTINE init_topography(dim_state_p, dim_state)

  use mpi
  use mod_parallel_pdaf, &
       only: MPIerr, mype_world
  Use fesom_pdaf, &
       ONLY: volo_full_glob, cellvol, topography_p, topography3D_g, topography3D, &
       inv_volo_full_glob, area_surf_glob, inv_area_surf_glob, mesh_fesom, nlmax, zbar_n_srf, zbar_n_bot, &
       myDim_nod2D, myDim_elem2D, MPI_COMM_FESOM, gather_nod
  USE statevector_pdaf, &
       ONLY: sfields, nfields

    IMPLICIT NONE

! *** Arguments
    INTEGER, INTENT(in) :: dim_state_p
    INTEGER, INTENT(in) :: dim_state

! *** Local variables
    INTEGER :: i, n, b, s, k
    INTEGER :: nz
    REAL, allocatable :: aux(:)  ! Temporary array


    ! *****************************
    ! ***    mesh topography    ***
    ! *****************************
    
    ! create mask:
    ! 1 -- valid model value
    ! 0 -- invalid because of topography
    
    ! array dimensions as pe-local model tracer fields
    allocate(topography3D(nlmax,myDim_nod2D))
    DO i=1, nlmax
      DO n=1, myDim_nod2D
        IF (mesh_fesom%areasvol(i,n)>0) THEN
          topography3D(i,n) = 1.0
        ELSE
          topography3D(i,n) = 0.0
        ENDIF
      ENDDO
    ENDDO
    
    ! array dimensions as global model tracer fields
    allocate(topography3D_g(nlmax,mesh_fesom% nod2D))
    call gather_nod(topography3D,topography3D_g)
    
    ! array dimensions as pe-local state vector
    allocate(topography_p(dim_state_p))
    DO b=1,nfields
      ! surface fields
      IF (sfields(b)%ndims==1) THEN
        topography_p(sfields(b)%off+1:sfields(b)%off+sfields(b)%dim) = 1.0
      ! 3D fields
      ELSEIF (sfields(b)%ndims==2) THEN
        DO i = 1, myDim_nod2D
        DO k = 1, nlmax
           s = (i-1) * nlmax + k
           topography_p(s + sfields(b)%off) = topography3D(k,i)
        ENDDO ! (k=1,nlmax)
        ENDDO ! (i=1,myDim_nod2D)
      ENDIF
    ENDDO ! (b=1,nfields)
    
    ! standard cell volume
    allocate(cellvol(nlmax,myDim_nod2D))
    cellvol=0.0
    do n=1, myDim_nod2D
       do nz=mesh_fesom%ulevels_nod2D(n), mesh_fesom%nlevels_nod2D(n)-1
          if     (nz==mesh_fesom%ulevels_nod2D(n)) then
                 ! surface
                 cellvol(nz,n) = mesh_fesom%areasvol(nz,n) * abs(mesh_fesom%zbar(nz+1) - zbar_n_srf(n))
          elseif (nz==mesh_fesom%nlevels_nod2D(n)-1) then
                 ! bottom
                 cellvol(nz,n) = mesh_fesom%areasvol(nz,n) * abs(zbar_n_bot(n) - mesh_fesom%zbar(nz))
          else
                 ! default
                 cellvol(nz,n) = mesh_fesom%areasvol(nz,n) * abs(mesh_fesom%zbar(nz+1) - mesh_fesom%zbar(nz))
          endif
       enddo
    enddo
    
    ! standard ocean volume
    allocate(aux(1))
    aux=0.0
    aux = sum(cellvol)
    call MPI_AllREDUCE(aux, volo_full_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    inv_volo_full_glob = 1.0/volo_full_glob
    deallocate(aux)
    
    ! ocean area
    allocate(aux(nlmax))
    do n=1, myDim_nod2D
       do nz=1,nlmax
          aux(nz)=aux(nz)+mesh_fesom%areasvol(nz,n)
       enddo
    enddo
    call MPI_AllREDUCE(aux, area_surf_glob, nlmax, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    inv_area_surf_glob = 1.0/area_surf_glob
    deallocate(aux)
    
  
! *** Screen output dimensions ***
    if (mype_world==0) THEN
       write(*,'(/a,2x,a)') 'FESOM-PDAF','*** Overview of model mesh dimensions ***'
       write(*,'(a,2x,a,1x,es13.6)') 'FESOM-PDAF','Standard ocean volume:', volo_full_glob
       write(*,'(a,2x,a,1x,es13.6)') 'FESOM-PDAF','Ocean area surface:   ', area_surf_glob(1)
       WRITE(*,'(a,2x,a,i12)') 'FESOM-PDAF','myDim_elem2D:      ', myDim_elem2D
       WRITE(*,'(a,2x,a,i12)') 'FESOM-PDAF','mesh_fesom%elem2D: ', mesh_fesom%elem2D
       WRITE(*,'(a,2x,a,i12)') 'FESOM-PDAF','mesh_fesom%nl:     ', mesh_fesom%nl
       WRITE(*,'(a,2x,a,i12)') 'FESOM-PDAF','nlmax:             ', nlmax
  endif

END SUBROUTINE init_topography

