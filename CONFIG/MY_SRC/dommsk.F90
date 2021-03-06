MODULE dommsk
   !!======================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask 
   !!======================================================================
   !! History :  OPA  ! 1987-07  (G. Madec)  Original code
   !!            6.0  ! 1993-03  (M. Guyon)  symetrical conditions (M. Guyon)
   !!            7.0  ! 1996-01  (G. Madec)  suppression of common work arrays
   !!             -   ! 1996-05  (G. Madec)  mask computed from tmask and sup-
   !!                 !                      pression of the double computation of bmask
   !!            8.0  ! 1997-02  (G. Madec)  mesh information put in domhgr.F
   !!            8.1  ! 1997-07  (G. Madec)  modification of mbathy and fmask
   !!             -   ! 1998-05  (G. Roullet)  free surface
   !!            8.2  ! 2000-03  (G. Madec)  no slip accurate
   !!             -   ! 2001-09  (J.-M. Molines)  Open boundaries
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_msk        : compute land/ocean mask
   !!   dom_msk_nsa    : update land/ocean mask when no-slip accurate option is used.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp
   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing
!{ DRAKKAR
   USE iom             ! For shlat2d
   USE fldread         ! for sn_shlat2d
!}

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_msk         ! routine called by inidom.F90
   PUBLIC   dom_msk_alloc   ! routine called by nemogcm.F90

   !                            !!* Namelist namlbc : lateral boundary condition *
   REAL(wp)        :: rn_shlat   = 2.   ! type of lateral boundary condition on velocity
   LOGICAL, PUBLIC :: ln_vorlat  = .false.   !  consistency of vorticity boundary condition 
   !                                            with analytical eqs.


   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  icoord ! Workspace for dom_msk_nsa()

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id: dommsk.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   INTEGER FUNCTION dom_msk_alloc()
      !!---------------------------------------------------------------------
      !!                 ***  FUNCTION dom_msk_alloc  ***
      !!---------------------------------------------------------------------
      dom_msk_alloc = 0
#if defined key_noslip_accurate
      ALLOCATE(icoord(jpi*jpj*jpk,3), STAT=dom_msk_alloc)
#endif
      IF( dom_msk_alloc /= 0 )   CALL ctl_warn('dom_msk_alloc: failed to allocate icoord array')
      !
   END FUNCTION dom_msk_alloc


   SUBROUTINE dom_msk
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zontal velocity points (u & v), vorticity points (f) and baro-
      !!      tropic stream function  points (b).
      !!
      !! ** Method  :   The ocean/land mask is computed from the basin bathy-
      !!      metry in level (mbathy) which is defined or read in dommba.
      !!      mbathy equals 0 over continental T-point 
      !!      and the number of ocean level over the ocean.
      !!
      !!      At a given position (ji,jj,jk) the ocean/land mask is given by:
      !!      t-point : 0. IF mbathy( ji ,jj) =< 0
      !!                1. IF mbathy( ji ,jj) >= jk
      !!      u-point : 0. IF mbathy( ji ,jj)  or mbathy(ji+1, jj ) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy(ji+1, jj ) >= jk.
      !!      v-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1) >= jk.
      !!      f-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1)
      !!                   or mbathy(ji+1,jj)  or mbathy(ji+1,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1)
      !!                  and mbathy(ji+1,jj) and mbathy(ji+1,jj+1) >= jk.
      !!      b-point : the same definition as for f-point of the first ocean
      !!                level (surface level) but with 0 along coastlines.
      !!      tmask_i : interior ocean mask at t-point, i.e. excluding duplicated
      !!                rows/lines due to cyclic or North Fold boundaries as well
      !!                as MPP halos.
      !!
      !!        The lateral friction is set through the value of fmask along
      !!      the coast and topography. This value is defined by rn_shlat, a
      !!      namelist parameter:
      !!         rn_shlat = 0, free slip  (no shear along the coast)
      !!         rn_shlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < rn_shlat < 2, partial slip   | non-linear velocity profile
      !!         2 < rn_shlat, strong slip        | in the lateral boundary layer
      !!
      !!      N.B. If nperio not equal to 0, the land/ocean mask arrays
      !!      are defined with the proper value at lateral domain boundaries,
      !!      but bmask. indeed, bmask defined the domain over which the
      !!      barotropic stream function is computed. this domain cannot
      !!      contain identical columns because the matrix associated with
      !!      the barotropic stream function equation is then no more inverti-
      !!      ble. therefore bmask is set to 0 along lateral domain boundaries
      !!      even IF nperio is not zero.
      !!
      !!      In case of open boundaries (lk_obc=T):
      !!        - tmask is set to 1 on the points to be computed bay the open
      !!          boundaries routines.
      !!        - bmask is  set to 0 on the open boundaries.
      !!
      !! ** Action :   tmask    : land/ocean mask at t-point (=0. or 1.)
      !!               umask    : land/ocean mask at u-point (=0. or 1.)
      !!               vmask    : land/ocean mask at v-point (=0. or 1.)
      !!               fmask    : land/ocean mask at f-point (=0. or 1.)
      !!                          =rn_shlat along lateral boundaries
      !!               bmask    : land/ocean mask at barotropic stream
      !!                          function point (=0. or 1.) and set to 0 along lateral boundaries
      !!               tmask_i  : interior ocean mask
      !!----------------------------------------------------------------------
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      INTEGER  ::   iif, iil, ii0, ii1, ii   ! local integers
      INTEGER  ::   ijf, ijl, ij0, ij1       !   -       -
      INTEGER , POINTER, DIMENSION(:,:) ::  imsk
      REAL(wp), POINTER, DIMENSION(:,:) ::  zwf

!{ DRAKKAR
      INTEGER  :: inum             !  logical unit for shlat2d
      REAL(wp) :: zshlat           !: locally modified shlat for some strait
      REAL(wp), POINTER, DIMENSION(:,:) :: zshlat2d
      LOGICAL  :: ln_shlat2d = .false.
      TYPE(FLD_N) :: sn_shlat2d
      !!
      NAMELIST/namlbc/ rn_shlat, ln_vorlat, ln_shlat2d, sn_shlat2d
!}
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dom_msk')
      !
      CALL wrk_alloc( jpi, jpj, imsk )
      CALL wrk_alloc( jpi, jpj, zwf  )
      !
      !  use fldread capability for shlat2d 
         !            !    file     ! frequency !  variable  ! time intep !  clim   ! 'yearly' or ! weights  ! rotation   !
         !            !    name     !  (hours)  !   name     !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs      !
      sn_shlat2d = FLD_N( 'shlat_coef',  -12    , 'shlatcoef', .false.    , .true.  ,  'yearly'   , ''       , ''         )

      REWIND( numnam )              ! Namelist namlbc : lateral momentum boundary condition
      READ  ( numnam, namlbc )
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dommsk : ocean mask '
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*) '   Namelist namlbc'
         WRITE(numout,*) '      lateral momentum boundary cond.    rn_shlat  = ',rn_shlat
         WRITE(numout,*) '      consistency with analytical form   ln_vorlat = ',ln_vorlat 
      ENDIF

      IF     (      rn_shlat == 0.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  free-slip '
      ELSEIF (      rn_shlat == 2.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  no-slip '
      ELSEIF ( 0. < rn_shlat .AND. rn_shlat < 2. ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  partial-slip '
      ELSEIF ( 2. < rn_shlat                     ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  strong-slip '
      ELSE
         WRITE(ctmp1,*) ' rn_shlat is negative = ', rn_shlat
         CALL ctl_stop( ctmp1 )
      ENDIF

!{ DRAKKAR
      IF ( ln_shlat2d ) THEN
         IF(lwp) WRITE(numout,*) '         READ shlat as a 2D coefficient in a file '
         CALL wrk_alloc( jpi, jpj, zshlat2d  )
         CALL iom_open(sn_shlat2d%clname, inum)
         CALL iom_get (inum, jpdom_data, sn_shlat2d%clvar, zshlat2d, 1) !
         CALL iom_close(inum)
      ENDIF
!}

      ! 1. Ocean/land mask at t-point (computed from mbathy)
      ! -----------------------------
      ! N.B. tmask has already the right boundary conditions since mbathy is ok
      !
      tmask(:,:,:) = 0._wp
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( REAL( mbathy(ji,jj) - jk, wp ) + 0.1_wp >= 0._wp )   tmask(ji,jj,jk) = 1._wp
            END DO  
         END DO  
      END DO  

!!gm  ????
#if defined key_zdfkpp
      IF( cp_cfg == 'orca' ) THEN
         IF( jp_cfg == 2 )   THEN       ! land point on Bab el Mandeb zonal section
            ij0 =  87   ;   ij1 =  88
            ii0 = 160   ;   ii1 = 161
            tmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 0._wp
         ELSE
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,cform_war)
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*)'          A mask must be applied on Bab el Mandeb strait'
            IF(lwp) WRITE(numout,*)'          in case of ORCAs configurations'
            IF(lwp) WRITE(numout,*)'          This is a problem which is not yet solved'
            IF(lwp) WRITE(numout,*)
         ENDIF
      ENDIF
#endif
!!gm end

      ! Interior domain mask (used for global sum)
      ! --------------------
      tmask_i(:,:) = tmask(:,:,1)
      iif = jpreci                         ! ???
      iil = nlci - jpreci + 1
      ijf = jprecj                         ! ???
      ijl = nlcj - jprecj + 1

      tmask_i( 1 :iif,   :   ) = 0._wp      ! first columns
      tmask_i(iil:jpi,   :   ) = 0._wp      ! last  columns (including mpp extra columns)
      tmask_i(   :   , 1 :ijf) = 0._wp      ! first rows
      tmask_i(   :   ,ijl:jpj) = 0._wp      ! last  rows (including mpp extra rows)

      ! north fold mask
      ! ---------------
      tpol(1:jpiglo) = 1._wp 
      fpol(1:jpiglo) = 1._wp
      IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
         tpol(jpiglo/2+1:jpiglo) = 0._wp
         fpol(     1    :jpiglo) = 0._wp
         IF( mjg(nlej) == jpjglo ) THEN                  ! only half of the nlcj-1 row
            DO ji = iif+1, iil-1
               tmask_i(ji,nlej-1) = tmask_i(ji,nlej-1) * tpol(mig(ji))
            END DO
         ENDIF
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
         tpol(     1    :jpiglo) = 0._wp
         fpol(jpiglo/2+1:jpiglo) = 0._wp
      ENDIF

      ! 2. Ocean/land mask at u-,  v-, and z-points (computed from tmask)
      ! -------------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector loop
               umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
               vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
            END DO
            DO ji = 1, jpim1      ! NO vector opt.
               fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
                  &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
            END DO
         END DO
      END DO
!{ DRAKKAR
      IF( cp_cfg == "orca" .AND. jp_cfg == 25 ) THEN  ! ORCA R025 configuration
          !                                            !  =======================
          ii0 = 212    ;   ii1 = 212        ! East of Ombai strait
          ij0 = 464    ;   ij1 = 465   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the East Ombai Strait'
          ii0 = 210    ;   ii1 = 211        ! West of Ombai strait
          ij0 = 466    ;   ij1 = 466   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the West Ombai Strait '
          ii0 = 210    ;   ii1 = 210        ! exit of Ombai strait
          ij0 = 464    ;   ij1 = 465   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the exit of Ombai Strait '
          ii0 = 172    ;   ii1 = 175        ! Lombok strait
          ij0 = 463    ;   ij1 = 463   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the Lombok Strait'
          ii0 = 278    ;   ii1 = 279        ! Torres Strait
          ij0 = 456    ;   ij1 = 461   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  4.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 4 at the Torres Strait'
 
          !
          ENDIF
!}
      CALL lbc_lnk( umask, 'U', 1._wp )      ! Lateral boundary conditions
      CALL lbc_lnk( vmask, 'V', 1._wp )
      CALL lbc_lnk( fmask, 'F', 1._wp )


      ! 4. ocean/land mask for the elliptic equation
      ! --------------------------------------------
      bmask(:,:) = tmask(:,:,1)       ! elliptic equation is written at t-point
      !
      !                               ! Boundary conditions
      !                                    ! cyclic east-west : bmask must be set to 0. on rows 1 and jpi
      IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
         bmask( 1 ,:) = 0._wp
         bmask(jpi,:) = 0._wp
      ENDIF
      IF( nperio == 2 ) THEN               ! south symmetric :  bmask must be set to 0. on row 1
         bmask(:, 1 ) = 0._wp
      ENDIF
      !                                    ! north fold : 
      IF( nperio == 3 .OR. nperio == 4 ) THEN   ! T-pt pivot : bmask set to 0. on row jpj and on half jpjglo-1 row
         DO ji = 1, jpi                      
            ii = ji + nimpp - 1
            bmask(ji,jpj-1) = bmask(ji,jpj-1) * tpol(ii)
            bmask(ji,jpj  ) = 0._wp
         END DO
      ENDIF
      IF( nperio == 5 .OR. nperio == 6 ) THEN   ! F-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj
         bmask(:,jpj) = 0._wp
      ENDIF
      !
      IF( lk_mpp ) THEN                    ! mpp specificities
         !                                      ! bmask is set to zero on the overlap region
         IF( nbondi /= -1 .AND. nbondi /= 2 )   bmask(  1 :jpreci,:) = 0._wp
         IF( nbondi /=  1 .AND. nbondi /= 2 )   bmask(nlci:jpi   ,:) = 0._wp
         IF( nbondj /= -1 .AND. nbondj /= 2 )   bmask(:,  1 :jprecj) = 0._wp
         IF( nbondj /=  1 .AND. nbondj /= 2 )   bmask(:,nlcj:jpj   ) = 0._wp
         !
         IF( npolj == 3 .OR. npolj == 4 ) THEN  ! north fold : bmask must be set to 0. on rows jpj-1 and jpj
            DO ji = 1, nlci
               ii = ji + nimpp - 1
               bmask(ji,nlcj-1) = bmask(ji,nlcj-1) * tpol(ii)
               bmask(ji,nlcj  ) = 0._wp
            END DO
         ENDIF
         IF( npolj == 5 .OR. npolj == 6 ) THEN  ! F-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj
            DO ji = 1, nlci
               bmask(ji,nlcj  ) = 0._wp
            END DO
         ENDIF
      ENDIF


      ! mask for second order calculation of vorticity
      ! ----------------------------------------------
      CALL dom_msk_nsa

      
      ! Lateral boundary conditions on velocity (modify fmask)
      ! ---------------------------------------     
      DO jk = 1, jpk
         zwf(:,:) = fmask(:,:,jk)         
!{ DRAKKAR
         IF (  ln_shlat2d ) THEN
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( fmask(ji,jj,jk) == 0. ) THEN
                     fmask(ji,jj,jk) = zshlat2d(ji,jj) * MIN( 1._wp, MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                        &                                              zwf(ji-1,jj), zwf(ji,jj-1)  )  )
                  ENDIF
               END DO
            END DO

         ELSE
!}

            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( fmask(ji,jj,jk) == 0._wp ) THEN
                     fmask(ji,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                        &                                           zwf(ji-1,jj), zwf(ji,jj-1)  )  )
                  ENDIF
               END DO
            END DO
!{ DRAKKAR
         ENDIF
!}

!{ DRAKKAR
       ! Locally modify shlat :
       IF( cp_cfg == "orca" .AND. jp_cfg == 025 ) THEN
         !                                           ! =======================
         ! Increased lateral friction in             !  ORCA_R025 configuration
         ! the vicinity of some straits              ! =======================
         !
            !! Gibraltar strait and Gulf of Cadiz
            ij0 = 652    ;   ij1 = 654
            ii0 = 1125   ;   ii1 = 1127
            zshlat=3
         DO jj = mj0(ij0),mj1(ij1)
            DO ji = mi0(ii0), mi1(ii1)
               IF( fmask(ji,jj,jk) == 0._wp ) THEN
                  fmask(ji,jj,jk) = zshlat * MIN( 1., MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                     &                                     zwf(ji-1,jj), zwf(ji,jj-1)  )  )
               ENDIF
            END DO
         END DO

         !
       ENDIF
!}
         DO jj = 2, jpjm1
            IF( fmask(1,jj,jk) == 0._wp ) THEN
               fmask(1  ,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(2,jj), zwf(1,jj+1), zwf(1,jj-1) ) )
            ENDIF
            IF( fmask(jpi,jj,jk) == 0._wp ) THEN
               fmask(jpi,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(jpi,jj+1), zwf(jpim1,jj), zwf(jpi,jj-1) ) )
            ENDIF
         END DO         
         DO ji = 2, jpim1
            IF( fmask(ji,1,jk) == 0._wp ) THEN
               fmask(ji, 1 ,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,1), zwf(ji,2), zwf(ji-1,1) ) )
            ENDIF
            IF( fmask(ji,jpj,jk) == 0._wp ) THEN
               fmask(ji,jpj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jpj), zwf(ji-1,jpj), zwf(ji,jpjm1) ) )
            ENDIF
         END DO
      END DO
      !
      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN   ! ORCA_R2 configuration
         !                                                 ! Increased lateral friction near of some straits
         IF( nn_cla == 0 ) THEN
            !                                ! Gibraltar strait  : partial slip (fmask=0.5)
            ij0 = 101   ;   ij1 = 101
            ii0 = 139   ;   ii1 = 140   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
            ij0 = 102   ;   ij1 = 102
            ii0 = 139   ;   ii1 = 140   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
            !
            !                                ! Bab el Mandeb : partial slip (fmask=1)
            ij0 =  87   ;   ij1 =  88
            ii0 = 160   ;   ii1 = 160   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
            ij0 =  88   ;   ij1 =  88
            ii0 = 159   ;   ii1 = 159   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
            !
         ENDIF
         !                                ! Danish straits  : strong slip (fmask > 2)
! We keep this as an example but it is instable in this case 
!         ij0 = 115   ;   ij1 = 115
!         ii0 = 145   ;   ii1 = 146   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
!         ij0 = 116   ;   ij1 = 116
!         ii0 = 145   ;   ii1 = 146   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
         !
      ENDIF
      !
      IF( cp_cfg == "orca" .AND. jp_cfg == 1 ) THEN   ! ORCA R1 configuration
         !                                                 ! Increased lateral friction near of some straits
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   orca_r1: increase friction near the following straits : '
         IF(lwp) WRITE(numout,*) '      Gibraltar '
         ii0 = 283   ;   ii1 = 284        ! Gibraltar Strait 
         ij0 = 200   ;   ij1 = 200   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2._wp  

         IF(lwp) WRITE(numout,*) '      Bhosporus '
         ii0 = 314   ;   ii1 = 315        ! Bhosporus Strait 
         ij0 = 208   ;   ij1 = 208   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2._wp  

         IF(lwp) WRITE(numout,*) '      Makassar (Top) '
         ii0 =  48   ;   ii1 =  48        ! Makassar Strait (Top) 
         ij0 = 149   ;   ij1 = 150   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  3._wp  

         IF(lwp) WRITE(numout,*) '      Lombok '
         ii0 =  44   ;   ii1 =  44        ! Lombok Strait 
         ij0 = 124   ;   ij1 = 125   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2._wp  

         IF(lwp) WRITE(numout,*) '      Ombai '
         ii0 =  53   ;   ii1 =  53        ! Ombai Strait 
         ij0 = 124   ;   ij1 = 125   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Timor Passage '
         ii0 =  56   ;   ii1 =  56        ! Timor Passage 
         ij0 = 124   ;   ij1 = 125   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      West Halmahera '
         ii0 =  58   ;   ii1 =  58        ! West Halmahera Strait 
         ij0 = 141   ;   ij1 = 142   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) = 3._wp  

         IF(lwp) WRITE(numout,*) '      East Halmahera '
         ii0 =  55   ;   ii1 =  55        ! East Halmahera Strait 
         ij0 = 141   ;   ij1 = 142   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) = 3._wp  
         !
      ENDIF
      !
      CALL lbc_lnk( fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask

      ! CAUTION : The fmask may be further modified in dyn_vor_init ( dynvor.F90 )
            
      IF( nprint == 1 .AND. lwp ) THEN      ! Control print
         imsk(:,:) = INT( tmask_i(:,:) )
         WRITE(numout,*) ' tmask_i : '
         CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                           1, jpj, 1, 1, numout)
         WRITE (numout,*)
         WRITE (numout,*) ' dommsk: tmask for each level'
         WRITE (numout,*) ' ----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( tmask(:,:,jk) )

            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout)
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: vmask for each level'
         WRITE(numout,*) ' -----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( vmask(:,:,jk) )
            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout)
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: fmask for each level'
         WRITE(numout,*) ' -----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( fmask(:,:,jk) )
            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout )
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: bmask '
         WRITE(numout,*) ' ---------------'
         WRITE(numout,*)
         imsk(:,:) = INT( bmask(:,:) )
         CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
            &                              1, jpj, 1, 1, numout )
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, imsk )
      CALL wrk_dealloc( jpi, jpj, zwf  )
!{ DRAKKAR
      IF ( ln_shlat2d ) THEN
        CALL wrk_dealloc( jpi, jpj, zshlat2d  )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dom_msk')
      !
   END SUBROUTINE dom_msk

#if defined key_noslip_accurate
   !!----------------------------------------------------------------------
   !!   'key_noslip_accurate' :         accurate no-slip boundary condition
   !!----------------------------------------------------------------------
   
   SUBROUTINE dom_msk_nsa
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk_nsa  ***
      !! 
      !! ** Purpose :
      !!
      !! ** Method  :
      !!
      !! ** Action :
      !!
      !! History :
      !!        !  00-03  (G. Madec)  no slip accurate
      !!        !  07-07  (G. Hervieux : debug + J.M. Molines : mpp )
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl      ! dummy loop indices
      INTEGER  ::   ii                  ! 
      INTEGER  ::   ine, inw, ins, inn, itest, ierror, iind, ijnd
      INTEGER  ::   ines, inws, inen, inwn
      REAL(wp) ::   zaa
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dom_msk_nsa')
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) 'dom_msk_nsa : noslip accurate boundary condition'
      IF(lwp)WRITE(numout,*) '~~~~~~~~~~~   using Schchepetkin and O Brian scheme'
      IF(lwp)WRITE(numout,*) '~~~~~~~~~~~   Gaelle Hervieux version '

      ! mask for second order calculation of vorticity
      ! ----------------------------------------------
      ! noslip boundary condition: fmask=1  at convex corner, store
      ! index of straight coast meshes ( 'west', refering to a coast,
      ! means west of the ocean, aso)

      DO jk = 1, jpk
         DO jl = 1, 8
            npcoa(jl,jk) = 0
            DO ji = 1, 2*(jpi+jpj)
               nicoa(ji,jl,jk) = 0
               njcoa(ji,jl,jk) = 0
            END DO
         END DO
      END DO

      IF( jperio == 2 ) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' symetric boundary conditions need special'
         WRITE(numout,*) ' treatment not implemented. we stop.'
         STOP
      ENDIF
      ! convex corners
! west-east straight coast
      DO jk = 1, jpkm1
         inw = 0
         ine = 0
         inws = 0
         ines = 0
         inwn = 0
         inen = 0
        DO jj = 2, jpj-1
          DO ji = 2, jpi-1
          ! west straight coast
          IF ( tmask(ji,jj,jk) == 0 .AND. tmask(ji+1,jj,jk) == 1 ) THEN
          ! west-south straight coast
             IF ( tmask(ji,jj+1,jk) == 1 .AND. tmask(ji+1,jj+1,jk) == 1  ) THEN
               inws = inws + 1
                   nicoa(inws,5,jk) = ji
                   njcoa(inws,5,jk) = jj
          ! west-north straight coast
             ELSEIF( tmask(ji,jj-1,jk) == 1 .AND. tmask(ji+1,jj-1,jk) == 1) THEN
               inwn = inwn + 1
                   nicoa(inwn,6,jk) = ji
                   njcoa(inwn,6,jk) = jj-1
             ELSE
               inw = inw + 1
                   nicoa(inw,1,jk) = ji
                   njcoa(inw,1,jk) = jj
             ENDIF
               ! east straight coast
          ELSEIF(tmask(ji,jj,jk)== 1 .AND.tmask(ji+1,jj,jk) == 0 ) THEN
          ! east-south straight coast
             IF    (tmask(ji,jj+1,jk) == 1 .AND. tmask(ji+1,jj+1,jk) == 1 ) THEN
               ines = ines + 1
                   nicoa(ines,7,jk) = ji
                   njcoa(ines,7,jk) = jj
          ! east-north straight coast
             ELSEIF(tmask(ji,jj-1,jk) == 1 .AND. tmask(ji+1,jj-1,jk) == 1 ) THEN
               inen = inen + 1
                   nicoa(inen,8,jk) = ji
                   njcoa(inen,8,jk) = jj-1
             ELSE
               ine = ine + 1
                   nicoa(ine,2,jk) = ji
                   njcoa(ine,2,jk) = jj
             ENDIF
          ENDIF
          END DO
         END DO
         npcoa(1,jk) = inw
         npcoa(2,jk) = ine
         npcoa(5,jk) = inws
         npcoa(6,jk) = inwn
         npcoa(7,jk) = ines
         npcoa(8,jk) = inen
      END DO

! south-north straight coast
      DO jk = 1, jpkm1
         ins = 0
         inn = 0
         DO ji = 2, jpi-1
          DO jj = 2, jpj-1
          ! south straight coast
          IF (tmask(ji  ,jj,jk) == 0 .AND. tmask(ji,jj+1,jk) == 1 &
           .AND.  tmask(ji+1,jj+1,jk) == 1 .AND.  tmask(ji-1,jj+1,jk) == 1) THEN
             ins = ins + 1
                 nicoa(ins,3,jk) = ji
                 njcoa(ins,3,jk) = jj
           ! north straight coast
          ELSEIF(tmask(ji  ,jj  ,jk)== 1 .AND. tmask(ji,jj+1,jk) == 0 &
            .AND.  tmask(ji+1,jj,jk) == 1 .AND.  tmask(ji-1,jj,jk) == 1 ) THEN
              inn = inn + 1
                 nicoa(inn,4,jk) = ji
                 njcoa(inn,4,jk) = jj
          ENDIF
          END DO
         END DO
         npcoa(3,jk) = ins
         npcoa(4,jk) = inn
      END DO
      ! Control print
      ! -------------
      itest = 2 * ( jpi + jpj )
      DO jk = 1, jpk
         IF( npcoa(1,jk) > itest .OR. npcoa(2,jk) > itest .OR.   &
             npcoa(3,jk) > itest .OR. npcoa(4,jk) > itest ) THEN
            WRITE(numout,*) ' nicoa(1,4,:) , njcoa =, ', nicoa(1,4,jk), njcoa(1,4,jk)
            WRITE(numout,*) ' nicoa(1,3,:) , njcoa =, ', nicoa(1,3,jk), njcoa(1,3,jk)
            WRITE(numout,*) ' nicoa(1,1,:) , njcoa =, ', nicoa(1,1,jk), njcoa(1,1,jk)
            WRITE(numout,*) ' nicoa(1,2,:) , njcoa =, ', nicoa(1,2,jk), njcoa(1,2,jk)
            WRITE(numout,*) ' '
            CALL ctl_stop( 'We stop...' )
        ENDIF
      END DO
      !
      IF( nn_timing == 1 )  CALL timing_stop('dom_msk_nsa')

   END SUBROUTINE dom_msk_nsa



#else
   !!----------------------------------------------------------------------
   !!   Default option :                                      Empty routine
   !!----------------------------------------------------------------------
   SUBROUTINE dom_msk_nsa       
   END SUBROUTINE dom_msk_nsa
#endif
   
   !!======================================================================
END MODULE dommsk
