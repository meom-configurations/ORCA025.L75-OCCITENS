MODULE stopar
   !!======================================================================
   !!                       ***  MODULE  stopar  ***
   !! Stochastic parameters : definition and time stepping
   !!=====================================================================
   !! History :  3.3  ! 2011-10 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sto_par       : update the stochastic parameters
   !!   sto_par_init  : define the stochastic parameterization
   !!   sto_rst_read  : read restart file for stochastic parameters
   !!   sto_rst_write : write restart file for stochastic parameters
   !!   sto_par_white : fill input array with white Gaussian noise
   !!   sto_par_flt   : apply horizontal Laplacian filter to input array
   !!----------------------------------------------------------------------
   USE storng          ! random number generator (external module)
   USE par_oce         ! ocean parameters
   USE dom_oce         ! ocean space and time domain variables
   USE lbclnk          ! lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module


   IMPLICIT NONE
   PRIVATE

   PUBLIC   sto_par_init    ! called by nemogcm.F90 
   PUBLIC   sto_par         ! called by step.F90 
   PUBLIC   sto_rst_write   ! called by step.F90 

   LOGICAL            :: ln_rststo = .FALSE.  ! restart stochastic parameters from restart file
   CHARACTER(len=255) :: cn_storst_in = "restart_sto"     ! suffix of sto restart name (input)
   CHARACTER(len=255) :: cn_storst_out = "restart_sto"    ! suffix of sto restart name (output)
   INTEGER            :: numstor, numstow     ! logical unit for restart (read and write)

   INTEGER           :: jpsto2d = 0          ! number of 2D stochastic parameters
   INTEGER           :: jpsto3d = 0          ! number of 3D stochastic parameters

   REAL(wp), PUBLIC, DIMENSION(:,:,:),   ALLOCATABLE :: sto2d      ! 2D stochastic parameters
   REAL(wp), PUBLIC, DIMENSION(:,:,:,:), ALLOCATABLE :: sto3d      ! 3D stochastic parameters
   REAL(wp),         DIMENSION(:,:),     ALLOCATABLE :: sto_tmp    ! temporary workspace
   REAL(wp),         DIMENSION(:,:),     ALLOCATABLE :: sto2d_abc  ! a, b, c parameters (for 2D arrays)
   REAL(wp),         DIMENSION(:,:),     ALLOCATABLE :: sto3d_abc  ! a, b, c parameters (for 3D arrays)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto2d_ave  ! mean value (for 2D arrays)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto3d_ave  ! mean value (for 3D arrays)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto2d_std  ! standard deviation (for 2D arrays)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto3d_std  ! standard deviation (for 3D arrays)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto2d_lim  ! limitation factor (for 2D arrays)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto3d_lim  ! limitation factor (for 3D arrays)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto2d_tcor ! time correlation (for 2D arrays)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto3d_tcor ! time correlation (for 3D arrays)
   INTEGER,          DIMENSION(:),       ALLOCATABLE :: sto2d_ord  ! order of autoregressive process
   INTEGER,          DIMENSION(:),       ALLOCATABLE :: sto3d_ord  ! order of autoregressive process

   CHARACTER(len=1), DIMENSION(:),       ALLOCATABLE :: sto2d_typ  ! nature of grid point (T, U, V, W, F, I)
   CHARACTER(len=1), DIMENSION(:),       ALLOCATABLE :: sto3d_typ  ! nature of grid point (T, U, V, W, F, I)
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto2d_sgn  ! control of the sign accross the north fold
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto3d_sgn  ! control of the sign accross the north fold
   INTEGER,          DIMENSION(:),       ALLOCATABLE :: sto2d_flt  ! number of passes of Laplacian filter
   INTEGER,          DIMENSION(:),       ALLOCATABLE :: sto3d_flt  ! number of passes of Laplacian filter
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto2d_fac  ! factor to restore std after filtering
   REAL(wp),         DIMENSION(:),       ALLOCATABLE :: sto3d_fac  ! factor to restore std after filtering

   LOGICAL, PUBLIC :: ln_sto_ldf = .FALSE.    ! stochastic lateral diffusion
   INTEGER, PUBLIC :: jsto_ldf                ! index of lateral diffusion stochastic parameter
   REAL(wp)        :: rn_ldf_std              ! lateral diffusion standard deviation (in percent)
   REAL(wp)        :: rn_ldf_tcor             ! lateral diffusion correlation timescale (in timesteps)

   LOGICAL, PUBLIC :: ln_sto_hpg = .FALSE.    ! stochastic horizontal pressure gradient
   INTEGER, PUBLIC :: jsto_hpgi               ! index of stochastic hpg parameter (i direction)
   INTEGER, PUBLIC :: jsto_hpgj               ! index of stochastic hpg parameter (j direction)
   REAL(wp)        :: rn_hpg_std              ! density gradient standard deviation (in percent)
   REAL(wp)        :: rn_hpg_tcor             ! density gradient correlation timescale (in timesteps)

   LOGICAL, PUBLIC :: ln_sto_eos = .FALSE.    ! stochastic equation of state
   INTEGER, PUBLIC :: nn_sto_eos = 1          ! number of degrees of freedom in stochastic equation of state
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosi ! index of stochastic eos parameter (i direction)
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosj ! index of stochastic eos parameter (j direction)
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: jsto_eosk ! index of stochastic eos parameter (k direction)
   REAL(wp)        :: rn_eos_stdxy            ! random walk horz. standard deviation (in grid points)
   REAL(wp)        :: rn_eos_stdz             ! random walk vert. standard deviation (in grid points)
   REAL(wp)        :: rn_eos_tcor             ! random walk correlation timescale (in timesteps)
   REAL(wp)        :: rn_eos_lim = 20.0_wp    ! limitation factor
   INTEGER         :: nn_eos_flt = 0          ! number of passes of Laplacian filter
   INTEGER         :: nn_eos_ord = 1          ! order of autoregressive processes

   ! Public array with density correction
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: drho_ran

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynhpg.F90 2528 2010-12-27 17:33:53Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sto_par( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_par  ***
      !! 
      !! ** Purpose :   update the stochastic parameters
      !!
      !! ** Method  :   model basic stochastic parameters
      !!                as a first order autoregressive process AR(1),
      !!                governed by the equation:
      !!                   X(t) = a * X(t-1) + b * w + c
      !!                where the parameters a, b and c are related
      !!                to expected value, standard deviation
      !!                and time correlation (all stationary in time) by:
      !!                   E   [X(t)]        = c / ( 1 - a )
      !!                   STD [X(t)]        = b / SQRT( 1 - a * a )
      !!                   COR [X(t),X(t-k)] = a ** k
      !!                and w is a Gaussian white noise.
      !!
      !!                Higher order autoregressive proces can be optionally generated
      !!                by replacing the white noise by a lower order process.
      !!
      !!                1) The statistics of the stochastic parameters (X) are assumed
      !!                constant in space (homogeneous) and time (stationary).
      !!                This could be generalized by replacing the constant
      !!                a, b, c parameters by functions of space and time.
      !!                
      !!                2) The computation is performed independently for every model
      !!                grid point, which corresponds to assume that the stochastic
      !!                parameters are uncorrelated in space.
      !!                This could be generalized by including a spatial filter: Y = Filt[ X ]
      !!                (possibly non-homgeneous and non-stationary) in the computation,
      !!                or by solving an elliptic equation: L[ Y ] = X.
      !!
      !!                3) The stochastic model for the parameters could also
      !!                be generalized to depend on the current state of the ocean (not done here).
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  :: ji, jj, jk, jsto, jflt
      REAL(wp) :: stomax

      !
      ! Update 2D stochastic arrays
      !
      DO jsto = 1, jpsto2d
        ! Store array from previous time step
        sto_tmp(:,:) = sto2d(:,:,jsto)

        IF ( sto2d_ord(jsto) == 1 ) THEN
          ! Draw new random numbers from N(0,1) --> w
          CALL sto_par_white( sto2d(:,:,jsto) )
          ! Apply horizontal Laplacian filter to w
          DO jflt = 1, sto2d_flt(jsto)
            CALL lbc_lnk( sto2d(:,:,jsto), sto2d_typ(jsto), sto2d_sgn(jsto) )
            CALL sto_par_flt( sto2d(:,:,jsto) )
          END DO
          ! Factor to restore standard deviation after filtering
          sto2d(:,:,jsto) = sto2d(:,:,jsto) * sto2d_fac(jsto)
        ELSE
          ! Use previous process (one order lower) instead of white noise
          sto2d(:,:,jsto) = sto2d(:,:,jsto-1)
        ENDIF

        ! Multiply white noise (or lower order process) by b --> b * w
        sto2d(:,:,jsto) = sto2d(:,:,jsto) * sto2d_abc(jsto,2)
        ! Update autoregressive processes --> a * X(t-1) + b * w
        sto2d(:,:,jsto) = sto2d(:,:,jsto) + sto_tmp(:,:) * sto2d_abc(jsto,1)
        ! Add parameter c --> a * X(t-1) + b * w + c
        sto2d(:,:,jsto) = sto2d(:,:,jsto) + sto2d_abc(jsto,3)
        ! Limit random parameters to std times the limitation factor
        stomax = sto2d_std(jsto) * sto2d_lim(jsto)
        sto2d(:,:,jsto) = SIGN(MIN(stomax,ABS(sto2d(:,:,jsto))),sto2d(:,:,jsto))

        ! Lateral boundary conditions on sto2d
        CALL lbc_lnk( sto2d(:,:,jsto), sto2d_typ(jsto), sto2d_sgn(jsto) )
      END DO
      !
      ! Update 3D stochastic arrays
      !
      DO jsto = 1, jpsto3d
         DO jk = 1, jpk
           ! Store array from previous time step
           sto_tmp(:,:) = sto3d(:,:,jk,jsto)

           IF ( sto3d_ord(jsto) == 1 ) THEN
             ! Draw new random numbers from N(0,1) --> w
             CALL sto_par_white( sto3d(:,:,jk,jsto) )
             ! Apply horizontal Laplacian filter to w
             DO jflt = 1, sto3d_flt(jsto)
               CALL lbc_lnk( sto3d(:,:,jk,jsto), sto3d_typ(jsto), sto3d_sgn(jsto) )
               CALL sto_par_flt( sto3d(:,:,jk,jsto) )
             END DO
             ! Factor to restore standard deviation after filtering
             sto3d(:,:,jk,jsto) = sto3d(:,:,jk,jsto) * sto3d_fac(jsto)
           ELSE
             ! Use previous process (one order lower) instead of white noise
             sto3d(:,:,jk,jsto) = sto3d(:,:,jk,jsto-1)
           ENDIF

           ! Multiply white noise by b --> b * w
           sto3d(:,:,jk,jsto) = sto3d(:,:,jk,jsto) * sto3d_abc(jsto,2)
           ! Update autoregressive processes --> a * X(t-1) + b * w
           sto3d(:,:,jk,jsto) = sto3d(:,:,jk,jsto) + sto_tmp(:,:) * sto3d_abc(jsto,1)
           ! Add parameter c --> a * X(t-1) + b * w + c
           sto3d(:,:,jk,jsto) = sto3d(:,:,jk,jsto) + sto3d_abc(jsto,3)
           ! Limit random parameters to std times the limitation factor
           stomax = sto3d_std(jsto) * sto3d_lim(jsto)
           sto3d(:,:,jk,jsto) = SIGN(MIN(stomax,ABS(sto3d(:,:,jk,jsto))),sto3d(:,:,jk,jsto))
         END DO
         ! Lateral boundary conditions on sto3d
         CALL lbc_lnk( sto3d(:,:,:,jsto), sto3d_typ(jsto), sto3d_sgn(jsto) )
      END DO

   END SUBROUTINE sto_par


   SUBROUTINE sto_par_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_par_init  ***
      !! 
      !! ** Purpose :   define the stochastic parameterization
      !!----------------------------------------------------------------------
      NAMELIST/namsto/ ln_sto_ldf, rn_ldf_std, rn_ldf_tcor, &
        &              ln_sto_hpg, rn_hpg_std, rn_hpg_tcor, &
        &              ln_sto_eos, nn_sto_eos, rn_eos_stdxy, rn_eos_stdz, &
        &              rn_eos_tcor, nn_eos_ord, nn_eos_flt, rn_eos_lim, &
        &              ln_rststo, cn_storst_in, cn_storst_out
       CHARACTER(LEN=40) :: cl_no
      !!----------------------------------------------------------------------
      INTEGER :: jsto, jmem, jarea, jdof, jord, jordm1, jk, jflt
      INTEGER(KIND=8) :: zseed1, zseed2, zseed3, zseed4
      REAL(wp) :: rinflate

      ! Read namsto namelist : stochastic parameterization
      REWIND( numnam )
      READ  ( numnam, namsto )
!{ JMM
      WRITE(cl_no,*) nn_no-1 ; cl_no = TRIM(ADJUSTL(cl_no) )
      cn_storst_in= TRIM(cn_storst_in)//'-'//TRIM(cl_no)

      IF (ln_ens_rst_in) THEN 
         cn_storst_in = TRIM(crstdir_in)//"/"//TRIM(cn_storst_in)//"."//TRIM(c_nmem)
      ELSE
         cn_storst_in = TRIM(crstdir_in)//"/"//TRIM(cn_storst_in)
      ENDIF
!}

      ! Parameter print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sto_par_init : stochastic parameterization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsto : stochastic parameterization'
         WRITE(numout,*) '      restart stochastic parameters           ln_rststo     = ', ln_rststo
         WRITE(numout,*) '      suffix of sto restart name (input)      cn_storst_in  = ', cn_storst_in
         WRITE(numout,*) '      suffix of sto restart name (output)     cn_storst_out = ', cn_storst_out
         WRITE(numout,*) '      stochastic lateral diffusion            ln_sto_ldf    = ', ln_sto_ldf
         WRITE(numout,*) '      lateral diffusion std (in percent)      rn_ldf_std    = ', rn_ldf_std
         WRITE(numout,*) '      lateral diffusion tcor (in timesteps)   rn_ldf_tcor   = ', rn_ldf_tcor
         WRITE(numout,*) '      stochastic horizontal pressure gradient ln_sto_hpg    = ', ln_sto_hpg
         WRITE(numout,*) '      density gradient std (in percent)       rn_hpg_std    = ', rn_hpg_std
         WRITE(numout,*) '      density gradient tcor (in timesteps)    rn_hpg_tcor   = ', rn_hpg_tcor
         WRITE(numout,*) '      stochastic equation of state            ln_sto_eos    = ', ln_sto_eos
         WRITE(numout,*) '      number of degrees of freedom            nn_sto_eos    = ', nn_sto_eos
         WRITE(numout,*) '      random walk horz. std (in grid points)  rn_eos_stdxy  = ', rn_eos_stdxy
         WRITE(numout,*) '      random walk vert. std (in grid points)  rn_eos_stdz   = ', rn_eos_stdz
         WRITE(numout,*) '      random walk tcor (in timesteps)         rn_eos_tcor   = ', rn_eos_tcor
         WRITE(numout,*) '      order of autoregressive  processes      nn_eos_ord    = ', nn_eos_ord
         WRITE(numout,*) '      passes of Laplacian filter              nn_eos_flt    = ', nn_eos_flt
         WRITE(numout,*) '      limitation factor                       rn_eos_lim    = ', rn_eos_lim
      ENDIF

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   stochastic parameterization :'

      ! Set number of 2D stochastic arrays
      jpsto2d = 0
      IF( ln_sto_ldf ) THEN
         IF(lwp) WRITE(numout,*) '       - stochastic lateral diffusion'
         jpsto2d   = jpsto2d + 1
         jsto_ldf  = jpsto2d
      ENDIF
      IF( ln_sto_eos ) THEN
         IF(lwp) WRITE(numout,*) '       - stochastic equation of state'
         ALLOCATE(jsto_eosi(nn_sto_eos))
         ALLOCATE(jsto_eosj(nn_sto_eos))
         ALLOCATE(jsto_eosk(nn_sto_eos))
         DO jdof = 1, nn_sto_eos
            jpsto2d   = jpsto2d + 3 * nn_eos_ord
            jsto_eosi(jdof) = jpsto2d - 2 * nn_eos_ord
            jsto_eosj(jdof) = jpsto2d - 1 * nn_eos_ord
            jsto_eosk(jdof) = jpsto2d
         END DO
      ENDIF

      ! Set number of 3D stochastic arrays
      jpsto3d = 0
      IF( ln_sto_hpg ) THEN
         IF(lwp) WRITE(numout,*) '       - stochastic horizontal pressure gradient'
         jpsto3d   = jpsto3d + 2
         jsto_hpgi = jpsto3d - 1
         jsto_hpgj = jpsto3d
      ENDIF

      ! Allocate 2D stochastic arrays
      IF ( jpsto2d > 0 ) THEN
         ALLOCATE ( sto2d(jpi,jpj,jpsto2d) )
         ALLOCATE ( sto2d_abc(jpsto2d,3) )
         ALLOCATE ( sto2d_ave(jpsto2d) )
         ALLOCATE ( sto2d_std(jpsto2d) )
         ALLOCATE ( sto2d_lim(jpsto2d) )
         ALLOCATE ( sto2d_tcor(jpsto2d) )
         ALLOCATE ( sto2d_ord(jpsto2d) )
         ALLOCATE ( sto2d_typ(jpsto2d) )
         ALLOCATE ( sto2d_sgn(jpsto2d) )
         ALLOCATE ( sto2d_flt(jpsto2d) )
         ALLOCATE ( sto2d_fac(jpsto2d) )
      ENDIF

      ! Allocate 3D stochastic arrays
      IF ( jpsto3d > 0 ) THEN
         ALLOCATE ( sto3d(jpi,jpj,jpk,jpsto3d) )
         ALLOCATE ( sto3d_abc(jpsto3d,3) )
         ALLOCATE ( sto3d_ave(jpsto3d) )
         ALLOCATE ( sto3d_std(jpsto3d) )
         ALLOCATE ( sto3d_lim(jpsto3d) )
         ALLOCATE ( sto3d_tcor(jpsto3d) )
         ALLOCATE ( sto3d_ord(jpsto3d) )
         ALLOCATE ( sto3d_typ(jpsto3d) )
         ALLOCATE ( sto3d_sgn(jpsto3d) )
         ALLOCATE ( sto3d_flt(jpsto3d) )
         ALLOCATE ( sto3d_fac(jpsto3d) )
      ENDIF

      ! Allocate temporary workspace
      IF ( jpsto2d > 0 .OR. jpsto3d > 0 ) THEN
         ALLOCATE ( sto_tmp(jpi,jpj) ) ; sto_tmp(:,:) = 0._wp
      ENDIF

      ! 1) For every stochastic parameter:
      ! ----------------------------------
      ! - set nature of grid point and control of the sign
      !       across the north fold (sto2d_typ, sto2d_sgn)
      ! - set number of passes of Laplacian filter (sto2d_flt)
      ! - set order of every autoregressive process (sto2d_ord)
      DO jsto = 1, jpsto2d
         sto2d_typ(jsto) = 'T'
         sto2d_sgn(jsto) = 1._wp
         sto2d_flt(jsto) = 0
         sto2d_ord(jsto) = 1
         DO jdof = 1, nn_sto_eos
         DO jord = 0, nn_eos_ord-1
            IF ( jsto+jord == jsto_eosi(jdof) ) THEN ! Stochastic equation of state i (ave=0)
               sto2d_ord(jsto) = nn_eos_ord - jord
               sto2d_sgn(jsto) = -1._wp
               sto2d_flt(jsto) = nn_eos_flt
            ENDIF
            IF ( jsto+jord == jsto_eosj(jdof) ) THEN ! Stochastic equation of state j (ave=0)
               sto2d_ord(jsto) = nn_eos_ord - jord
               sto2d_sgn(jsto) = -1._wp
               sto2d_flt(jsto) = nn_eos_flt
            ENDIF
            IF ( jsto+jord == jsto_eosk(jdof) ) THEN ! Stochastic equation of state k (ave=0)
               sto2d_ord(jsto) = nn_eos_ord - jord
               sto2d_flt(jsto) = nn_eos_flt
            ENDIF
         END DO
         END DO
         sto2d_fac(jsto) = sto_par_flt_fac ( sto2d_flt(jsto) )
      END DO
      !
      DO jsto = 1, jpsto3d
         sto3d_typ(jsto) = 'T'
         sto3d_sgn(jsto) = 1._wp
         sto3d_flt(jsto) = 0
         sto3d_ord(jsto) = 1
         IF ( jsto == jsto_hpgi ) THEN ! Stochastic density gradient i (ave=1)
            sto3d_typ(jsto) = 'U'
         ENDIF
         IF ( jsto == jsto_hpgj ) THEN ! Stochastic density gradient j (ave=1)
            sto3d_typ(jsto) = 'V'
         ENDIF
         sto3d_fac(jsto) = sto_par_flt_fac ( sto3d_flt(jsto) )
      END DO

      ! 2) For every stochastic parameter:
      ! ----------------------------------
      ! set average, standard deviation and time correlation
      DO jsto = 1, jpsto2d
         sto2d_ave(jsto)  = 0._wp
         sto2d_std(jsto)  = 1._wp
         sto2d_tcor(jsto) = 1._wp
         sto2d_lim(jsto)  = 3._wp
         IF ( jsto == jsto_ldf  ) THEN ! Stochastic lateral diffusion (ave=1)
            sto2d_ave(jsto)  = 1._wp
            sto2d_std(jsto)  = rn_ldf_std
            sto2d_tcor(jsto) = rn_ldf_tcor
         ENDIF
         DO jdof = 1, nn_sto_eos
         DO jord = 0, nn_eos_ord-1
            IF ( jsto+jord == jsto_eosi(jdof) ) THEN ! Stochastic equation of state i (ave=0)
               sto2d_std(jsto)  = rn_eos_stdxy
               sto2d_tcor(jsto) = rn_eos_tcor
               sto2d_lim(jsto)  = rn_eos_lim
            ENDIF
            IF ( jsto+jord == jsto_eosj(jdof) ) THEN ! Stochastic equation of state j (ave=0)
               sto2d_std(jsto)  = rn_eos_stdxy
               sto2d_tcor(jsto) = rn_eos_tcor
               sto2d_lim(jsto)  = rn_eos_lim
            ENDIF
            IF ( jsto+jord == jsto_eosk(jdof) ) THEN ! Stochastic equation of state k (ave=0)
               sto2d_std(jsto)  = rn_eos_stdz
               sto2d_tcor(jsto) = rn_eos_tcor
               sto2d_lim(jsto)  = rn_eos_lim
            ENDIF
         END DO
         END DO
      END DO
      !
      DO jsto = 1, jpsto3d
         sto3d_ave(jsto)  = 0._wp
         sto3d_std(jsto)  = 1._wp
         sto3d_tcor(jsto) = 1._wp
         sto3d_lim(jsto)  = 3._wp
         IF ( jsto == jsto_hpgi ) THEN ! Stochastic density gradient i (ave=1)
            sto3d_ave(jsto)  = 1._wp
            sto3d_std(jsto)  = rn_hpg_std
            sto3d_tcor(jsto) = rn_hpg_tcor
         ENDIF
         IF ( jsto == jsto_hpgj ) THEN ! Stochastic density gradient j (ave=1)
            sto3d_ave(jsto)  = 1._wp
            sto3d_std(jsto)  = rn_hpg_std
            sto3d_tcor(jsto) = rn_hpg_tcor
         ENDIF
      END DO

      ! 3) For every stochastic parameter:
      ! ----------------------------------
      ! - compute parameters (a, b, c) of the AR1 autoregressive process
      !   from expected value (ave), standard deviation (std)
      !   and time correlation (tcor):
      !     a = EXP ( - 1 / tcor )           --> sto2d_abc(:,1)
      !     b = std * SQRT( 1 - a * a )      --> sto2d_abc(:,2)
      !     c = ave * ( 1 - a )              --> sto2d_abc(:,3)
      ! - for higher order processes (ARn, n>1), use approximate formula
      !   for the b parameter (valid for tcor>>1 time step)
      DO jsto = 1, jpsto2d
         IF ( sto2d_tcor(jsto) == 0._wp ) THEN
            sto2d_abc(jsto,1) = 0._wp
         ELSE
            sto2d_abc(jsto,1) = EXP ( - 1._wp / sto2d_tcor(jsto) )
         ENDIF
         IF ( sto2d_ord(jsto) == 1 ) THEN      ! Exact formula for 1st order process
            rinflate = sto2d_std(jsto)
         ELSE
            ! Approximate formula, valid for tcor >> 1
            jordm1 = sto2d_ord(jsto) - 1
            rinflate = SQRT ( REAL( jordm1 , wp ) / REAL( 2*(2*jordm1-1) , wp ) )
         ENDIF
         sto2d_abc(jsto,2) = rinflate * SQRT ( 1._wp - sto2d_abc(jsto,1) &
                                                     * sto2d_abc(jsto,1) )
         sto2d_abc(jsto,3) = sto2d_ave(jsto) * ( 1._wp - sto2d_abc(jsto,1) )
      END DO
      !
      DO jsto = 1, jpsto3d
         IF ( sto3d_tcor(jsto) == 0._wp ) THEN
            sto3d_abc(jsto,1) = 0._wp
         ELSE
            sto3d_abc(jsto,1) = EXP ( - 1._wp / sto3d_tcor(jsto) )
         ENDIF
         IF ( sto3d_ord(jsto) == 1 ) THEN      ! Exact formula for 1st order process
            rinflate = sto3d_std(jsto)
         ELSE
            ! Approximate formula, valid for tcor >> 1
            jordm1 = sto3d_ord(jsto) - 1
            rinflate = SQRT ( REAL( jordm1 , wp ) / REAL( 2*(2*jordm1-1) , wp ) )
         ENDIF
         sto3d_abc(jsto,2) = rinflate * SQRT ( 1._wp - sto3d_abc(jsto,1) &
                                                     * sto3d_abc(jsto,1) )
         sto3d_abc(jsto,3) = sto2d_ave(jsto) * ( 1._wp - sto3d_abc(jsto,1) )
      END DO

      ! 4) Initialize seeds for random number generator
      ! -----------------------------------------------
      ! using different seeds for different processors (jarea)
      ! and different ensemble members (jmem)
      CALL kiss_reset( )
      DO jarea = 1, narea
         DO jmem = 0, nmember
            zseed1 = kiss() ; zseed2 = kiss() ; zseed3 = kiss() ; zseed4 = kiss()
         END DO
      END DO
      CALL kiss_seed( zseed1, zseed2, zseed3, zseed4 )

      ! 5) Initialize stochastic parameters to: ave + std * w
      ! -----------------------------------------------------
      DO jsto = 1, jpsto2d
         ! Draw random numbers from N(0,1) --> w
         CALL sto_par_white( sto2d(:,:,jsto) )
         ! Apply horizontal Laplacian filter to w
         DO jflt = 1, sto2d_flt(jsto)
            CALL lbc_lnk( sto2d(:,:,jsto), sto2d_typ(jsto), sto2d_sgn(jsto) )
            CALL sto_par_flt( sto2d(:,:,jsto) )
         END DO
         ! Multiply by standard devation and add average value
         sto2d(:,:,jsto) = sto2d(:,:,jsto) * sto2d_std(jsto) + sto2d_ave(jsto)
      END DO
      !
      DO jsto = 1, jpsto3d
         DO jk = 1, jpk
            ! Draw random numbers from N(0,1) --> w
            CALL sto_par_white( sto3d(:,:,jk,jsto) )
            ! Apply horizontal Laplacian filter to w
            DO jflt = 1, sto3d_flt(jsto)
               CALL lbc_lnk( sto3d(:,:,jk,jsto), sto3d_typ(jsto), sto3d_sgn(jsto) )
               CALL sto_par_flt( sto3d(:,:,jk,jsto) )
            END DO
            ! Multiply by standard devation and add average value
            sto3d(:,:,jk,jsto) = sto3d(:,:,jk,jsto) * sto3d_std(jsto) + sto3d_ave(jsto)
         END DO
      END DO

      ! 6) Restart stochastic parameters from file
      ! ------------------------------------------
      IF( ln_rststo ) CALL sto_rst_read

      ! Allocate drho_ran
      ALLOCATE(drho_ran(jpi,jpj,jpk))

   END SUBROUTINE sto_par_init


   SUBROUTINE sto_rst_read
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_rst_read  ***
      !! 

      !! ** Purpose :   read stochastic parameters from restart file
      !!----------------------------------------------------------------------

      INTEGER  :: jsto, jseed
      INTEGER(KIND=8)     ::   ziseed(4)           ! RNG seeds in integer type
      REAL(KIND=8)        ::   zrseed(4)           ! RNG seeds in real type (with same bits to save in restart)
      CHARACTER(LEN=9)    ::   clsto2d='sto2d_000' ! stochastic parameter variable name
      CHARACTER(LEN=9)    ::   clsto3d='sto3d_000' ! stochastic parameter variable name
      CHARACTER(LEN=10)   ::   clseed='seed0_0000' ! seed variable name

      IF ( jpsto2d > 0 .OR. jpsto3d > 0 ) THEN

         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sto_rst_read : read stochastic parameters from restart file'
            WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF

         ! Open the restart file
         CALL iom_open( cn_storst_in, numstor, kiolib = jprstlib )

         ! Get stochastic parameters from restart file:
         ! 2D stochastic parameters
         DO jsto = 1 , jpsto2d
            WRITE(clsto2d(7:9),'(i3.3)') jsto
            CALL iom_get( numstor, jpdom_autoglo, clsto2d , sto2d(:,:,jsto) )
         END DO
         ! 3D stochastic parameters
         DO jsto = 1 , jpsto3d
            WRITE(clsto3d(7:9),'(i3.3)') jsto
            CALL iom_get( numstor, jpdom_autoglo, clsto3d , sto3d(:,:,:,jsto) )
         END DO

         IF (ln_ens_rst_in) THEN
            ! Get saved state of the random number generator
            DO jseed = 1 , 4
               WRITE(clseed(5:5) ,'(i1.1)') jseed
               WRITE(clseed(7:10),'(i4.4)') narea
               CALL iom_get( numstor, clseed , zrseed(jseed) )
            END DO
            ziseed = TRANSFER( zrseed , ziseed)
            CALL kiss_seed( ziseed(1) , ziseed(2) , ziseed(3) , ziseed(4) )
         ENDIF

         ! Close the restart file
         CALL iom_close( numstor )

      ENDIF

   END SUBROUTINE sto_rst_read


   SUBROUTINE sto_rst_write( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_rst_write  ***
      !! 
      !! ** Purpose :   write stochastic parameters in restart file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !!
      INTEGER  :: jsto, jseed
      INTEGER(KIND=8)     ::   ziseed(4)           ! RNG seeds in integer type
      REAL(KIND=8)        ::   zrseed(4)           ! RNG seeds in real type (with same bits to save in restart)
      CHARACTER(LEN=20)   ::   clkt                ! ocean time-step defined as a character
      CHARACTER(LEN=255)  ::   clname              ! restart file name
      CHARACTER(LEN=9)    ::   clsto2d='sto2d_000' ! stochastic parameter variable name
      CHARACTER(LEN=9)    ::   clsto3d='sto3d_000' ! stochastic parameter variable name
      CHARACTER(LEN=10)   ::   clseed='seed0_0000' ! seed variable name

      IF ( jpsto2d > 0 .OR. jpsto3d > 0 ) THEN

         IF( kt == nitrst .OR. kt == nitend ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'sto_rst_write : write stochastic parameters in restart file'
               WRITE(numout,*) '~~~~~~~~~~~~~'
            ENDIF
         ENDIF

         ! Put stochastic parameters in restart files
         ! (as opened at previous timestep, see below)
         IF( kt > nit000) THEN
         IF( kt == nitrst .OR. kt == nitend ) THEN
            ! get and save current state of the random number generator
            CALL kiss_state( ziseed(1) , ziseed(2) , ziseed(3) , ziseed(4) )
            zrseed = TRANSFER( ziseed , zrseed)
            DO jseed = 1 , 4
               WRITE(clseed(5:5) ,'(i1.1)') jseed
               WRITE(clseed(7:10),'(i4.4)') narea
               CALL iom_rstput( kt, nitrst, numstow, clseed , zrseed(jseed) )
            END DO
            ! 2D stochastic parameters
            DO jsto = 1 , jpsto2d
               WRITE(clsto2d(7:9),'(i3.3)') jsto
               CALL iom_rstput( kt, nitrst, numstow, clsto2d , sto2d(:,:,jsto) )
            END DO
            ! 3D stochastic parameters
            DO jsto = 1 , jpsto3d
               WRITE(clsto3d(7:9),'(i3.3)') jsto
               CALL iom_rstput( kt, nitrst, numstow, clsto3d , sto3d(:,:,:,jsto) )
            END DO
            ! Save drho_ran in restart file
            CALL iom_rstput( kt, nitrst, numstow, 'drho' , drho_ran(:,:,:) )
            ! close the restart file
            CALL iom_close( numstow )
         ENDIF
         ENDIF

         ! Open the restart file one timestep before writing restart
         IF( kt < nitend) THEN
         IF( kt == nitrst - 1 .OR. nstock == 1 .OR. kt == nitend-1 ) THEN
            ! create the filename
            IF( nitrst > 999999999 ) THEN   ;   WRITE(clkt, *       ) nitrst
            ELSE                            ;   WRITE(clkt, '(i8.8)') nitrst
            ENDIF
            clname = TRIM(crstdir_out)//"/"//TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_storst_out)
            ! print information
            IF(lwp) THEN
               WRITE(numout,*) '             open stochastic parameters restart file: '//clname
               IF( kt == nitrst - 1 ) THEN
                  WRITE(numout,*) '             kt = nitrst - 1 = ', kt
               ELSE
                  WRITE(numout,*) '             kt = '             , kt
               ENDIF
            ENDIF
            ! open the restart file
            CALL iom_open( clname, numstow, ldwrt = .TRUE., kiolib = jprstlib )
         ENDIF
         ENDIF

      ENDIF

   END SUBROUTINE sto_rst_write


   SUBROUTINE sto_par_white( psto )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_par_white  ***
      !! 
      !! ** Purpose :   fill input array with white Gaussian noise
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out)           ::   psto
      !!
      INTEGER  :: ji, jj
      REAL(KIND=8) :: gran   ! Gaussian random number (forced KIND=8 as in kiss_gaussian)

      DO jj = 1, jpj
         DO ji = 1, jpi
            CALL kiss_gaussian( gran )
            psto(ji,jj) = gran
         END DO
      END DO

   END SUBROUTINE sto_par_white


   SUBROUTINE sto_par_flt( psto )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_par_flt  ***
      !! 
      !! ** Purpose :   apply horizontal Laplacian filter to input array
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out)           ::   psto
      !!
      INTEGER  :: ji, jj

      DO jj = 2, jpj-1
         DO ji = 2, jpi-1
            psto(ji,jj) = 0.5_wp * psto(ji,jj) + 0.125_wp * &
                              &  ( psto(ji-1,jj) + psto(ji+1,jj) +  &
                              &    psto(ji,jj-1) + psto(ji,jj+1) )
         END DO
      END DO

   END SUBROUTINE sto_par_flt


   REAL(wp) FUNCTION sto_par_flt_fac( kpasses )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION sto_par_flt_fac  ***
      !! 
      !! ** Purpose :   compute factor to restore standard deviation
      !!                as a function of the number of passes
      !!                of the Laplacian filter
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kpasses
      !!
      INTEGER :: jpasses, ji, jj, jflti, jfltj
      INTEGER, DIMENSION(-1:1,-1:1) :: pflt0
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pfltb
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pflta
      REAL(wp) :: ratio

      pflt0(-1,-1) = 0 ; pflt0(-1,0) = 1 ; pflt0(-1,1) = 0
      pflt0( 0,-1) = 1 ; pflt0( 0,0) = 4 ; pflt0( 0,1) = 1
      pflt0( 1,-1) = 0 ; pflt0( 1,0) = 1 ; pflt0( 1,1) = 0

      ALLOCATE(pfltb(-kpasses-1:kpasses+1,-kpasses-1:kpasses+1))
      ALLOCATE(pflta(-kpasses-1:kpasses+1,-kpasses-1:kpasses+1))

      pfltb(:,:) = 0
      pfltb(0,0) = 1
      DO jpasses = 1, kpasses
        pflta(:,:) = 0
        DO jflti= -1, 1
        DO jfltj= -1, 1
          DO ji= -kpasses, kpasses
          DO jj= -kpasses, kpasses
            pflta(ji,jj) = pflta(ji,jj) + pfltb(ji+jflti,jj+jfltj) * pflt0(jflti,jfltj)
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        pfltb(:,:) = pflta(:,:)
      ENDDO

      ratio = SUM(pfltb(:,:))
      ratio = ratio * ratio / SUM(pfltb(:,:)*pfltb(:,:))
      ratio = SQRT(ratio)

      DEALLOCATE(pfltb,pflta)

      sto_par_flt_fac = ratio

   END FUNCTION sto_par_flt_fac


END MODULE stopar

