MODULE obs_readmdt
   !!======================================================================
   !!                       ***  MODULE obs_readmdt  ***
   !! Observation diagnostics: Read the MDT for SLA data (skeleton for now)
   !!======================================================================
   !! History :      ! 2007-03 (K. Mogensen) Initial skeleton version
   !!                ! 2007-04 (E. Remy) migration and improvement from OPAVAR
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   obs_rea_mdt    : Driver for reading MDT
   !!   obs_offset_mdt : Remove the offset between the model MDT and the used one
   !!----------------------------------------------------------------------
   USE wrk_nemo         ! Memory Allocation
   USE par_kind         ! Precision variables
   USE par_oce          ! Domain parameters
   USE in_out_manager   ! I/O manager
   USE obs_surf_def     ! Surface observation definitions
   USE obs_inter_sup    ! Interpolation support routines
   USE obs_inter_h2d    ! 2D interpolation
   USE obs_utils        ! Various observation tools
   USE iom_nf90         ! IOM NetCDF
   USE netcdf           ! NetCDF library
   USE lib_mpp          ! MPP library
   USE dom_oce, ONLY : &                  ! Domain variables
      &                    tmask, tmask_i, e1t, e2t, gphit, glamt
   USE obs_const, ONLY :   obfillflt      ! Fillvalue
   USE oce      , ONLY :   sshn           ! Model variables

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   obs_rea_mdt     ! called by ?
   PUBLIC   obs_offset_mdt  ! called by ?

   INTEGER , PUBLIC ::   nmsshc    = 1         ! MDT correction scheme
   REAL(wp), PUBLIC ::   rmdtcorr   = 1.61_wp  ! User specified MDT correction
   REAL(wp), PUBLIC ::   rmdtcutoff = 65.0_wp  ! MDT cutoff for computed correction

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_readmdt.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE obs_rea_mdt( kslano, sladata, k2dint )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_rea_mdt ***
      !!
      !! ** Purpose : Read from file the MDT data (skeleton)
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!----------------------------------------------------------------------
      USE iom
      !
      INTEGER                          , INTENT(IN)    ::   kslano    ! Number of SLA Products
      TYPE(obs_surf), DIMENSION(kslano), INTENT(inout) ::   sladata   ! SLA data
      INTEGER                          , INTENT(in)    ::   k2dint    ! ?
      !
      CHARACTER(LEN=12), PARAMETER ::   cpname  = 'obs_rea_mdt'
      CHARACTER(LEN=20), PARAMETER ::   cpmdtname = 'slaReferenceLevel.nc'

      INTEGER ::   jslano              ! Data set loop variable
      INTEGER ::   jobs                ! Obs loop variable
      INTEGER ::   i_var_id, i_stat
      INTEGER ::   nummdt
      !
      REAL(wp), DIMENSION(1)     ::   zext, zobsmask
      REAL(wp), DIMENSION(2,2,1) ::   zweig
      !
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   zmask, zmdtl, zglam, zgphi
      INTEGER , DIMENSION(:,:,:), ALLOCATABLE ::   igrdi, igrdj
      !
      REAL(wp), POINTER, DIMENSION(:,:) ::  z_mdt, zmdtmask
         
      REAL(wp) :: zlam, zphi, zfill, zinfill    ! local scalar
      !!----------------------------------------------------------------------

      CALL wrk_alloc(jpi,jpj,z_mdt,zmdtmask) 

      IF(lwp)WRITE(numout,*) 
      IF(lwp)WRITE(numout,*) ' obs_rea_mdt : Read MDT for referencing altimeter anomalies'
      IF(lwp)WRITE(numout,*) ' ------------- '
      IF(lwp)WRITE(numout,*) '      REFERENCE LEVEL SET TO 0 ( HARD WIRED IN OCCIPUT ) '

!     CALL iom_open( cpmdtname, nummdt )       ! Open the file
      !                                      ! Get the MDT data
!     CALL iom_get ( nummdt, jpdom_data, 'sossheig', z_mdt(:,:), 1 )
!     CALL iom_close(nummdt)                 ! Close the file
      z_mdt(:,:) = 0._wp           ! OCCIPUT
      
!     ! Read in the fill value
!     zinfill = 0.0
!     i_stat = nf90_open( cpmdtname, nf90_nowrite, nummdt )
!     i_stat = nf90_inq_varid( nummdt, 'sossheig', i_var_id )
!     i_stat = nf90_get_att( nummdt, i_var_id, "_FillValue",zinfill)
!     zfill = zinfill
!     i_stat = nf90_close( nummdt )

      ! setup mask based on tmask and MDT mask
      ! set mask to 0 where the MDT is set to fillvalue
!     WHERE(z_mdt(:,:) /= zfill)   ;   zmdtmask(:,:) = tmask(:,:,1)
!     ELSE WHERE                   ;   zmdtmask(:,:) = 0
!     END WHERE
      zmdtmask(:,:) = tmask(:,:,1)   ! OCCIPUT

      ! Remove the offset between the MDT used with the sla and the model MDT
      IF( nmsshc == 1 .OR. nmsshc == 2 )   CALL obs_offset_mdt( z_mdt, zfill )

      ! Intepolate the MDT already on the model grid at the observation point
  
      DO jslano = 1, kslano
         ALLOCATE( &
            & igrdi(2,2,sladata(jslano)%nsurf), &
            & igrdj(2,2,sladata(jslano)%nsurf), &
            & zglam(2,2,sladata(jslano)%nsurf), &
            & zgphi(2,2,sladata(jslano)%nsurf), &
            & zmask(2,2,sladata(jslano)%nsurf), &
            & zmdtl(2,2,sladata(jslano)%nsurf)  &
            & )
         
         DO jobs = 1, sladata(jslano)%nsurf

            igrdi(1,1,jobs) = sladata(jslano)%mi(jobs)-1
            igrdj(1,1,jobs) = sladata(jslano)%mj(jobs)-1
            igrdi(1,2,jobs) = sladata(jslano)%mi(jobs)-1
            igrdj(1,2,jobs) = sladata(jslano)%mj(jobs)
            igrdi(2,1,jobs) = sladata(jslano)%mi(jobs)
            igrdj(2,1,jobs) = sladata(jslano)%mj(jobs)-1
            igrdi(2,2,jobs) = sladata(jslano)%mi(jobs)
            igrdj(2,2,jobs) = sladata(jslano)%mj(jobs)

         END DO

         CALL obs_int_comm_2d( 2, 2, sladata(jslano)%nsurf, igrdi, igrdj, glamt   , zglam )
         CALL obs_int_comm_2d( 2, 2, sladata(jslano)%nsurf, igrdi, igrdj, gphit   , zgphi )
         CALL obs_int_comm_2d( 2, 2, sladata(jslano)%nsurf, igrdi, igrdj, zmdtmask, zmask )
         CALL obs_int_comm_2d( 2, 2, sladata(jslano)%nsurf, igrdi, igrdj, z_mdt   , zmdtl )

         DO jobs = 1, sladata(jslano)%nsurf
            
            zlam = sladata(jslano)%rlam(jobs)
            zphi = sladata(jslano)%rphi(jobs)

            CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
               &                   zglam(:,:,jobs), zgphi(:,:,jobs), &
               &                   zmask(:,:,jobs), zweig, zobsmask )
            
            CALL obs_int_h2d( 1, 1, zweig, zmdtl(:,:,jobs),  zext )
 
            sladata(jslano)%rext(jobs,2) = zext(1)

! mark any masked data with a QC flag
            IF( zobsmask(1) == 0 )   sladata(jslano)%nqc(jobs) = 11

         END DO
         
         DEALLOCATE( &
            & igrdi, &
            & igrdj, &
            & zglam, &
            & zgphi, &
            & zmask, &
            & zmdtl  &
            & )

      END DO

      CALL wrk_dealloc(jpi,jpj,z_mdt,zmdtmask) 
      !
   END SUBROUTINE obs_rea_mdt


   SUBROUTINE obs_offset_mdt( pmdt, pfill )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_offset_mdt ***
      !!
      !! ** Purpose : Compute a correction term for the MDT on the model grid
      !!             !!!!! IF it is on the model grid
      !!
      !! ** Method  : Compute the mean difference between the model and the 
      !!              used MDT and remove the offset. 
      !!
      !! ** Action  : 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pmdt     ! MDT used on the model grid
      REAL(wp)                    , INTENT(in   ) ::   pfill 
      ! 
      INTEGER  :: ji, jj
      REAL(wp) :: zdxdy, zarea, zeta1, zeta2, zcorr_mdt, zcorr_bcketa, zcorr     ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:) :: zpromsk
      CHARACTER(LEN=14), PARAMETER ::   cpname = 'obs_offset_mdt'
      !!----------------------------------------------------------------------

      CALL wrk_alloc( jpi,jpj, zpromsk )

      !  Initialize the local mask, for domain projection 
      !  Also exclude pmdt points which are set to missing data

      DO ji = 1, jpi
        DO jj = 1, jpj
           zpromsk(ji,jj) = tmask_i(ji,jj)
           IF (    ( gphit(ji,jj) .GT.  rmdtcutoff ) &
              &.OR.( gphit(ji,jj) .LT. -rmdtcutoff ) &
              &.OR.( pmdt(ji,jj)  .EQ.  pfill ) ) &
              &        zpromsk(ji,jj) = 0.0
        END DO
      END DO 

      ! Compute MSSH mean over [0,360] x [-rmdtcutoff,rmdtcutoff]

      zarea = 0.0
      zeta1 = 0.0
      zeta2 = 0.0

      DO jj = 1, jpj
         DO ji = 1, jpi
          zdxdy = e1t(ji,jj) * e2t(ji,jj) * zpromsk(ji,jj)
          zarea = zarea + zdxdy
          zeta1 = zeta1 + pmdt(ji,jj) * zdxdy
          zeta2 = zeta2 + sshn (ji,jj) * zdxdy
        END DO      
      END DO

      IF( lk_mpp)   CALL mpp_sum( zeta1 )
      IF( lk_mpp)   CALL mpp_sum( zeta2 )
      IF( lk_mpp)   CALL mpp_sum( zarea )
      
      zcorr_mdt    = zeta1 / zarea
      zcorr_bcketa = zeta2 / zarea

      !  Define correction term

      zcorr = zcorr_mdt - zcorr_bcketa

      !  Correct spatial mean of the MSSH

      IF( nmsshc == 1 )   pmdt(:,:) = pmdt(:,:) - zcorr  

      ! User defined value : 1.6 m for the Rio MDT compared to ORCA2 MDT

      IF( nmsshc == 2 )   pmdt(:,:) = pmdt(:,:) - rmdtcorr

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' obs_readmdt : rmdtcutoff     = ', rmdtcutoff
         WRITE(numout,*) ' -----------   zcorr_mdt     = ', zcorr_mdt
         WRITE(numout,*) '               zcorr_bcketa  = ', zcorr_bcketa
         WRITE(numout,*) '               zcorr         = ', zcorr
         WRITE(numout,*) '               nmsshc        = ', nmsshc
      ENDIF

      IF ( nmsshc == 0 ) WRITE(numout,*) '           MSSH correction is not applied'
      IF ( nmsshc == 1 ) WRITE(numout,*) '           MSSH correction is applied'
      IF ( nmsshc == 2 ) WRITE(numout,*) '           User defined MSSH correction' 

      CALL wrk_dealloc( jpi,jpj, zpromsk )
      !
   END SUBROUTINE obs_offset_mdt
 
   !!======================================================================
END MODULE obs_readmdt
