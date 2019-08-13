!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chemistry_mod.F90
!
! !DESCRIPTION: Module CHEMISTRY\_MOD is used to call the proper chemistry 
!  subroutine for the various GEOS-Chem simulations. 
!\\
!\\
! !INTERFACE:
!
MODULE Chemistry_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)
  USE Geos_Timers_Mod  ! For GEOS-Chem timers (optional)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: INIT_CHEMISTRY
  PUBLIC  :: DO_CHEMISTRY
  PUBLIC  :: RECOMPUTE_OD
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CHEM_PASSIVE_SPECIES
!
! !REVISION HISTORY: 
!  (1 ) Bug fix in DO_CHEMISTRY (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Now references "tagged_ox_mod.f"(bmy, 8/18/03)
!  (4 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (5 ) Bug fix: Now also call OPTDEPTH for GEOS-4 (bmy, 1/27/04)
!  (6 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/5/04)
!  (7 ) Now references "seasalt_mod.f" (rjp, bec, bmy, 4/20/04)
!  (8 ) Now references "logical_mod.f", "tracer_mod.f", "diag20_mod.f", and
!        "diag65_mod.f", and "aerosol_mod." (bmy, 7/20/04)
!  (9 ) Now references "mercury_mod.f" (bmy, 12/7/04)
!  (10) Updated for SO4s, NITs chemistry (bec, bmy, 4/13/05)
!  (11) Now call CHEM_HCN_CH3CN from "hcn_ch3cn_mod.f".  Also remove all
!        references to the obsolete CO-OH param simulation. (xyp, bmy, 6/24/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now call MAKE_RH from "main.f" (bmy, 3/16/06)
!  (14) Updated for SOA from isoprene (dkh, bmy, 6/1/06)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (16) For now, replace use RPMARES instead of ISORROPIA. (bmy, 4/2/08)
!  (17) Added KPP chemistry driver subroutine (phs,ks,dhk, 09/15/09)
!  (18) Added public member function recompute_OD (skim, 02/03/11)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  28 Jan 2010 - C. Carouge, R. Yantosca - Modified for ISORROPIA II
!  08 Aug 2012 - R. Yantosca - Now align IF statements better
!  10 Aug 2012 - R. Yantosca - Cosmetic changes
!  25 Mar 2013 - M. Payer    - Now pass State_Chm to several routines
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 May 2014 - C. Keller   - Added INIT_CHEMISTRY
!  15 Dec 2014 - M. Yannetti - KPP code is commented out unless compiling KPP
!  08 Jan 2015 - M. Sulprizio- Now restrict KPP to REAL*8 to allow for KPP code
!                              to compile properly
!  13 Aug 2015 - E. Lundgren - Tracer units are now kg/kg and converted to
!                              kg within DO_CHEMISTRY
!  03 Nov 2016 - C. Keller   - Added wrapper routine for passive tracers.
!  17 Nov 2017 - R. Yantosca - Now in F90 format; added Diag_OH_HO2_O1D_O3P
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

  INTEGER :: id_DST1, id_NK1   ! Species ID flags

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_chemistry
!
! !DESCRIPTION: Subroutine DO\_CHEMISTRY is the driver routine which calls 
!  the appropriate chemistry subroutine for the various GEOS-Chem simulations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_CHEMISTRY( am_I_Root,  Input_Opt,  State_Chm,                &
                           State_Diag, State_Grid, State_Met, RC            )
!
! !USES:
!
    USE AEROSOL_MOD,     ONLY : AEROSOL_CONC
    USE AEROSOL_MOD,     ONLY : RDAER
    USE AEROSOL_MOD,     ONLY : SOILDUST
    USE C2H6_MOD,        ONLY : CHEMC2H6
    USE CARBON_MOD,      ONLY : CHEMCARBON
#if defined( BPCH_DIAG )
    USE CMN_DIAG_MOD  
#endif
    USE CMN_SIZE_MOD
    USE Diagnostics_Mod, ONLY : Compute_Column_Mass
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
    USE DUST_MOD,        ONLY : CHEMDUST
    USE DUST_MOD,        ONLY : RDUST_ONLINE
    USE ErrCode_Mod      
    USE ERROR_MOD        
    USE FlexChem_Mod,    ONLY : Do_FlexChem
    USE GLOBAL_CH4_MOD,  ONLY : CHEMCH4
    USE Input_Opt_Mod,   ONLY : OptInput
    USE ISORROPIAII_MOD, ONLY : DO_ISORROPIAII
    USE MERCURY_MOD,     ONLY : CHEMMERCURY
    USE POPS_MOD,        ONLY : CHEMPOPS
    USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
    USE RPMARES_MOD,     ONLY : DO_RPMARES
    USE SEASALT_MOD,     ONLY : CHEMSEASALT
    USE SULFATE_MOD,     ONLY : CHEMSULFATE
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Chm_Mod,   ONLY : Ind_
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE STRAT_CHEM_MOD,  ONLY : DO_STRAT_CHEM
    USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
    USE TAGGED_O3_MOD,   ONLY : CHEM_TAGGED_O3
    USE TIME_MOD,        ONLY : GET_ELAPSED_SEC
    USE TIME_MOD,        ONLY : GET_TS_CHEM
#if defined( USE_TEND )  
    USE TENDENCIES_MOD   
#endif                   
#if defined( APM )
    USE APM_INIT_MOD,    ONLY : APMIDS
    USE APM_DRIV_MOD,    ONLY : PSO4GAS
    USE APM_DRIV_MOD,    ONLY : AERONUM
    USE APM_DRIV_MOD,    ONLY : APM_DRIV
#endif
#if defined( TOMAS )     
    USE TOMAS_MOD,       ONLY : DO_TOMAS  !(win, 7/14/09)
#endif                   
    USE UCX_MOD,         ONLY : CALC_STRAT_AER
    USE UnitConv_Mod,    ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER            :: N_TROP, N
    INTEGER            :: MONTH
    INTEGER            :: YEAR
    INTEGER            :: WAVELENGTH
    LOGICAL            :: IT_IS_A_C2H6_SIM
    LOGICAL            :: IT_IS_A_CH3I_SIM
    LOGICAL            :: IT_IS_A_CH4_SIM
    LOGICAL            :: IT_IS_A_FULLCHEM_SIM
    LOGICAL            :: IT_IS_A_H2HD_SIM
    LOGICAL            :: IT_IS_A_HCN_SIM
    LOGICAL            :: IT_IS_A_MERCURY_SIM
    LOGICAL            :: IT_IS_A_RnPbBe_SIM
    LOGICAL            :: IT_IS_A_TAGCO_SIM
    LOGICAL            :: IT_IS_A_TAGO3_SIM
    LOGICAL            :: IT_IS_AN_AEROSOL_SIM
    LOGICAL            :: IT_IS_NOT_COPARAM_OR_CH4
    LOGICAL            :: IT_IS_A_POPS_SIM
    LOGICAL            :: LCARB
    LOGICAL            :: LCHEM
    LOGICAL            :: LDUST
    LOGICAL            :: LSCHEM
    LOGICAL            :: LPRT
    LOGICAL            :: LSSALT
    LOGICAL            :: LSULF
    LOGICAL            :: LSOA
    LOGICAL            :: LNLPBL
    LOGICAL            :: LUCX
    REAL(fp)           :: DT_Chem
#if defined( APM )
    INTEGER            :: I,J,L
    REAL*8             :: CONCTMPSO4(State_Grid%NX,                         &
                                     State_Grid%NY,                         &
                                     State_Grid%NZ)
#endif

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc

    !=======================================================================
    ! DO_CHEMISTRY begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Do_Chemistry  (in module GeosCore/chemistry_mod.F90)'

    ! Copy fields from INPUT_OPT to local variables for use below
    LCARB                    = Input_Opt%LCARB                        
    LCHEM                    = Input_Opt%LCHEM
    LDUST                    = Input_Opt%LDUST
    LSCHEM                   = Input_Opt%LSCHEM
    LPRT                     = Input_Opt%LPRT
    LSSALT                   = Input_Opt%LSSALT
    LSULF                    = Input_Opt%LSULF
    LSOA                     = Input_Opt%LSOA
    LNLPBL                   = Input_Opt%LNLPBL
    LUCX                     = Input_Opt%LUCX
    IT_IS_A_C2H6_SIM         = Input_Opt%ITS_A_C2H6_SIM
    IT_IS_A_CH3I_SIM         = Input_Opt%ITS_A_CH3I_SIM
    IT_IS_A_CH4_SIM          = Input_Opt%ITS_A_CH4_SIM 
    IT_IS_A_FULLCHEM_SIM     = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_A_H2HD_SIM         = Input_Opt%ITS_A_H2HD_SIM
    IT_IS_A_HCN_SIM          = Input_Opt%ITS_A_HCN_SIM
    IT_IS_A_MERCURY_SIM      = Input_Opt%ITS_A_MERCURY_SIM
    IT_IS_A_RnPbBe_SIM       = Input_Opt%ITS_A_RnPbBe_SIM
    IT_IS_A_TAGCO_SIM        = Input_Opt%ITS_A_TAGCO_SIM
    IT_IS_A_TAGO3_SIM        = Input_Opt%ITS_A_TAGO3_SIM
    IT_IS_A_POPS_SIM         = Input_Opt%ITS_A_POPS_SIM
    IT_IS_AN_AEROSOL_SIM     = Input_Opt%ITS_AN_AEROSOL_SIM
    
    ! Save species ID"s on first call
    IF ( FIRST ) THEN
       id_DST1 = Ind_('DST1')
       id_NK1  = Ind_('NK1' )
    ENDIF

    !----------------------------------------------------------
    ! Chemistry budget diagnostics - Part 1 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetChemistry ) THEN
       
       ! Get initial column masses
       CALL Compute_Column_Mass(                                             &
            am_I_Root   = am_I_Root,                                         &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            SpcMapping  = State_Chm%Map_Advect,                              &
            isFull      = State_Diag%Archive_BudgetChemistryFull,            &
            SpcMapFull  = State_Diag%Map_BudgetChemistryFull,                &
            ColMassFull = State_Diag%BudgetMassFull1,                        &
            isTrop      = State_Diag%Archive_BudgetChemistryTrop,            &
            SpcMapTrop  = State_Diag%Map_BudgetChemistryTrop,                &
            ColMassTrop = State_Diag%BudgetMassTrop1,                        &
            isPBL       = State_Diag%Archive_BudgetChemistryPBL,             &
            SpcMapPBL   = State_Diag%Map_BudgetChemistryPBL,                 & 
            ColMassPBL  = State_Diag%BudgetMassPBL1,                         &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Compute_Column_Mass" (initial)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

#if defined( USE_TEND )
    !=======================================================================
    ! Archive species concentrations for tendencies (ckeller,7/15/2015)
    !=======================================================================
    CALL Tend_Stage1( am_I_Root, Input_Opt, State_Chm,                       &
                      State_Met, 'CHEM', RC                                 )
#endif

    !=======================================================================
    ! Convert species units to [kg] for chemistry (ewl, 8/12/15)
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root,  Input_Opt, State_Chm,                &
                            State_Grid, State_Met, 'kg',                     &
                            RC,         OrigUnit=OrigUnit                   )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (kg/kg dry -> kg)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! If LCHEM=T then call the chemistry subroutines
    !=======================================================================
    IF ( LCHEM ) THEN

       !====================================================================
       ! Full-chemistry simulations:
       !
       ! (1) Benchmark; (2) Standard; (3) SimpleSOA; (4) complexSOA, 
       ! (5) complexSOA-SVPOA; (6) aciduptake; (7) marinePOA
       !====================================================================
       IF ( IT_IS_A_FULLCHEM_SIM ) THEN 
             
#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          !----------------------------------------
          ! Dry-run sulfate chem to get cloud pH
          !----------------------------------------
          IF ( LSULF ) THEN

             ! Dry run only
             CALL ChemSulfate( am_I_Root,  Input_Opt,  State_Chm,            &
                               State_Diag, State_Grid, State_Met,            &
                               .FALSE.,    RC                               )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

#if defined( APM )
          ! Save SO4 concentration before chemistry
          N          = APMIDS%id_SO4
          CONCTMPSO4 = State_Chm%Species(:,:,:,N)

          CALL AERONUM( am_I_Root,  Input_Opt,  State_Chm,                   &
                        State_Diag, State_Grid, State_Met, RC               )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemSulfate"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
#endif

          !---------------------------
          ! Call gas-phase chemistry
          !---------------------------
          CALL Do_FlexChem( am_I_Root,  Input_Opt,  State_Chm,               &
                            State_Diag, State_Grid, State_Met, RC           )

          ! Check units (ewl, 10/5/15)
          IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
             ErrMsg = 'Incorrect species units after FLEX_CHEMDR!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Do_FlexChem"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

          !----------------------------------------
          ! Call linearized stratospheric scheme
          !----------------------------------------
          IF ( LSCHEM ) THEN 

#if defined( USE_TIMERS )
             CALL GEOS_Timer_Start( "=> Strat chem", RC )
#endif

             ! Do linearized chemistry for the stratosphere (tropchem)
             ! or the mesosphere (UCX)
             CALL Do_Strat_Chem( am_I_Root,  Input_Opt, State_Chm,           &
                                 State_Grid, State_Met, RC                  )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg = 'Incorrect species units after DO_STRAT_CHEM!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in ""!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

#if defined( USE_TIMERS )
             CALL GEOS_Timer_End( "=> Strat chem", RC )
#endif

          ENDIF

#if defined( APM )
          ! Obtain SO4 production after chemistry
          N = APMIDS%id_SO4
          !$OMP PARALLEL DO         &
          !$OMP DEFAULT( SHARED   ) &
          !$OMP PRIVATE( I, J, L  ) &
          !$OMP SCHEDULE( DYNAMIC )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             IF ( State_Chm%Species(I,J,L,N) > CONCTMPSO4(I,J,L) ) THEN
                PSO4GAS(I,J,L) = State_Chm%Species(I,J,L,N)                  &
                               - CONCTMPSO4(I,J,L)
             ELSE
                PSO4GAS(I,J,L) = 0.D0
             ENDIF
          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
#endif

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

          !--------------------------------
          ! Do seasalt aerosol chemistry
          !--------------------------------
          IF ( LSSALT ) THEN
             CALL ChemSeaSalt( am_I_Root,  Input_Opt,  State_Chm,            &
                               State_Diag, State_Grid, State_Met, RC        )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSeaSalt"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !-------------------------------
          ! Recalculate PSC properties
          !-------------------------------
          IF ( LUCX ) THEN

#if defined( USE_TIMERS )
             CALL GEOS_Timer_End  ( "=> All aerosol chem", RC )
             CALL GEOS_Timer_Start( "=> Strat chem",       RC )
#endif
             
             ! Recalculate PSC
             CALL Calc_Strat_Aer( am_I_Root,  Input_Opt, State_Chm,          &
                                  State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Calc_Strat_Aer"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

#if defined( USE_TIMERS )
             CALL GEOS_Timer_End  ( "=> Strat chem",       RC )
             CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

          ENDIF

          !--------------------------------
          ! Also do sulfate chemistry
          !--------------------------------
          IF ( LSULF ) THEN

             ! Do sulfate chemistry
             CALL ChemSulfate( am_I_Root,  Input_Opt,  State_Chm,             &
                               State_Diag, State_Grid, State_Met,             &
                               .TRUE.,     RC                                )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg =  'Incorrect species units after CHEMSULFATE!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
       
             !-----------------------------------------
             ! Do aerosol thermodynamic equilibrium
             !-----------------------------------------
             IF ( LSSALT ) THEN

#if !defined( NO_ISORROPIA )
#if !defined( APM )
                ! ISORROPIA takes Na+, Cl- into account
                CALL Do_IsorropiaII( am_I_Root,  Input_Opt,  State_Chm,      &
                                     State_Diag, State_Grid, State_Met, RC  )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_ISORROPIAII"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
#endif
#endif

             ELSE

#if defined( APM )
                WRITE(*,*)'Warning: APM does not want to use DO_RPMARES'
                STOP
#endif

                ! RPMARES does not take Na+, Cl- into account
                CALL Do_RPMARES( am_I_Root,  Input_Opt, State_Chm,           &
                                 State_Grid, State_Met, RC                  )

             ENDIF

          ENDIF

          !-----------------------------------
          ! Do carbonaceous aerosol chemistry
          !-----------------------------------
          IF ( LCARB ) THEN
             CALL ChemCarbon( am_I_Root,  Input_Opt,  State_Chm,             &
                              State_Diag, State_Grid, State_Met, RC         )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemCarbon"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------
          ! Do dust aerosol chemistry/removal
          !------------------------------------
          IF ( LDUST .AND. id_DST1 > 0 ) THEN
             CALL ChemDust( am_I_Root,  Input_Opt,  State_Chm,               &
                            State_Diag, State_Grid, State_Met, RC           )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemDust"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
 
#if defined( APM )
          !--------------------------------------------
          ! Do APM aerosol microphysics
          !--------------------------------------------
          CALL APM_DRIV( am_I_Root,  Input_Opt,  State_Chm,                  &
                         State_Diag, State_Grid, State_Met, RC              )
                  
          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in routine "APM_DRIV"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
#endif

#if   defined( TOMAS )
          !--------------------------------------------
          ! Do TOMAS aerosol microphysics and dry dep
          !--------------------------------------------
          IF ( id_NK1 > 0 ) THEN 
             CALL Do_TOMAS( am_I_Root, Input_Opt,  State_Chm,               &
                           State_Diag, State_Grid, State_Met, RC           )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg = 'Incorrect species units after DO_TOMAS!' 
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_TOMAS"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
#endif

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> All aerosol chem", RC )
#endif
          
       !====================================================================
       ! Aerosol-only simulation
       !====================================================================
       ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

          !-------------------------------------------------------
          ! Compute aerosol & dust concentrations [kg/m3]
          ! (NOTE: SOILDUST in "aerosol_mod.f" is computed here)
          !-------------------------------------------------------
          CALL Aerosol_Conc( am_I_Root,  Input_Opt,  State_Chm,              &
                             State_Diag, State_Grid, State_Met, RC          )

          ! Check units (ewl, 10/5/15)
          IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
             ErrMsg = 'Incorrect species units after AEROSOL_CONC!'             
             CALL GC_Error( ErrMsg, RC, ThisLoc )
          ENDIF

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Aerosol_Conc"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !-------------------------------------------
          ! Compute AOD's and surface areas at 999 nm
          !-------------------------------------------
          MONTH      = 0
          YEAR       = 0
          WAVELENGTH = 0
          CALL RdAer( am_I_Root,  Input_Opt,  State_Chm,                     &
                      State_Diag, State_Grid, State_Met, RC,                 &
                      MONTH,      YEAR,       WAVELENGTH                    )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "RdAer"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !--------------------------------------------
          ! Aerosol Thermodynamic Equilibrium
          !--------------------------------------------
          IF ( LSULF ) THEN
             IF ( LSSALT ) THEN

#if !defined( NO_ISORROPIA )
#if !defined( APM )
                ! ISORROPIA takes Na+, Cl- into account
                CALL Do_IsorropiaII( am_I_Root,  Input_Opt,  State_Chm,      &
                                     State_Diag, State_Grid, State_Met, RC  )
#endif
#endif

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_IsorropiaII"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

             ELSE

#if defined( APM )
                WRITE(*,*)'Warning: APM does not want to use DO_RPMARES'
                STOP
#endif

                ! RPMARES does not take Na+, Cl- into account
                ! (skip for crystalline & aqueous offline run)
                CALL Do_RPMARES( am_I_Root,  Input_Opt, State_Chm,           &
                                 State_Grid, State_Met, RC                  )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_RPMARES"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
             ENDIF
          ENDIF

          !-----------------------------
          ! Seasalt Aerosols
          !-----------------------------
          IF ( LSSALT ) THEN
             CALL ChemSeaSalt( am_I_Root,  Input_Opt,  State_Chm,            &
                               State_Diag, State_Grid, State_Met, RC        )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSeaSalt"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !-------------------
          ! Sulfate aerosols
          !-------------------
          IF ( LSULF ) THEN
 
             ! Do sulfate chemistry
             CALL ChemSulfate( am_I_Root,  Input_Opt,  State_Chm,            &
                               State_Diag, State_Grid, State_Met,            &
                               .TRUE.,     RC                               )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
            
          !-----------------------------------------
          ! Carbon and Secondary Organic Aerosols
          !-----------------------------------------
          IF ( LCARB ) THEN
             CALL ChemCarbon( am_I_Root,  Input_Opt,  State_Chm,             &
                              State_Diag, State_Grid, State_Met, RC         )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in ""!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------
          ! Mineral Dust Aerosols
          !------------------------
          IF ( LDUST ) THEN 

             ! Do dust aerosol chemistry
             CALL ChemDust( am_I_Root,  Input_Opt,  State_Chm,               &
                            State_Diag, State_Grid, State_Met, RC           )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemDust"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Compute dust OD's & surface areas
             WAVELENGTH = 0
             CALL Rdust_Online( am_I_Root,  Input_Opt,  State_Chm,           &
                                State_Diag, State_Grid, State_Met,           &
                                SOILDUST,   WAVELENGTH, RC                  )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Rdust_Online"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> All aerosol chem", RC )
#endif

       !====================================================================
       ! Rn-Pb-Be
       !====================================================================
       ELSE IF ( IT_IS_A_RnPbBe_SIM ) THEN

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          ! Do Rn-Pb-Be chemistry
          CALL ChemRnPbBe( am_I_Root,  Input_Opt,  State_Chm,                &
                           State_Diag, State_Grid, State_Met, RC            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemRnPbBe"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

       !====================================================================
       ! Tagged O3
       !====================================================================
       ELSE IF ( IT_IS_A_TAGO3_SIM ) THEN 

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          !-----------------------------------------------
          ! Do Tagged O3 chemistry
          !-----------------------------------------------
          CALL Chem_Tagged_O3( am_I_Root,  Input_Opt,  State_Chm,            &
                               State_Diag, State_Grid, State_Met, RC        )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Tagged_O3"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          
#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

          !-----------------------------------------------
          ! Call linearized stratospheric scheme (LINOZ)
          !-----------------------------------------------
          IF ( LSCHEM ) THEN 

#if defined( USE_TIMERS )
             CALL GEOS_Timer_Start( "=> Strat chem", RC )
#endif

             ! Do LINOZ for Ozone
             CALL Do_Strat_Chem( am_I_Root,  Input_Opt, State_Chm,           &
                                 State_Grid, State_Met, RC                  )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_Strat_Chem"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

#if defined( USE_TIMERS )
             CALL GEOS_Timer_End( "=> Strat chem", RC )
#endif

          ENDIF

       !====================================================================
       ! Tagged CO
       !====================================================================
       ELSE IF ( IT_IS_A_TAGCO_SIM ) THEN

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          ! Do tagged CO chemistry
          CALL Chem_Tagged_CO( am_I_Root,  Input_Opt,  State_Chm,            &
                               State_Diag, State_Grid, State_Met, RC        )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Tagged_CO"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

       !====================================================================
       ! C2H6
       !====================================================================
       ELSE IF ( IT_IS_A_C2H6_SIM ) THEN
          CALL ChemC2H6( am_I_Root, Input_Opt, State_Chm, State_Grid,       &
                         State_Met, RC                                     )
 
          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemC2H6"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !====================================================================
       ! CH4
       !====================================================================
       ELSE IF ( IT_IS_A_CH4_SIM ) THEN

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif 

          CALL ChemCh4( am_I_Root,  Input_Opt,  State_Chm,                &
                        State_Diag, State_Grid, State_Met, RC            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemCh4"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

       !====================================================================
       ! Mercury
       !====================================================================
       ELSE IF ( IT_IS_A_MERCURY_SIM ) THEN
 
#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          ! Do Hg chemistry
          CALL ChemMercury( am_I_Root,  Input_Opt,  State_Chm,               &
                            State_Diag, State_Grid, State_Met, RC           )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemMercury"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

       !====================================================================
       ! POPs
       !====================================================================
       ELSE IF ( IT_IS_A_POPS_SIM ) THEN
 
#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          ! Do POPS chemistry
          CALL ChemPOPs( am_I_Root,  Input_Opt,  State_Chm,                  &
                         State_Diag, State_Grid, State_Met, RC              )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemPOPs"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif
       ENDIF

       !====================================================================
       ! PASSIVE SPECIES
       !
       ! This performs a simple loss chemistry on passive species.  Call 
       ! this routine for all simulation types since passive species can 
       ! be defined for various simulations (as additional species to the 
       ! default! ones). ckeller, 09/04/15
       !
       ! NOTE: To speed up execution, only call Chem_Passive_Species
       ! if there is at least one passive species with a finite 
       ! atmospheric lifetime.  There is no reason to apply a loss rate
       ! of unity to those passive species whose lifetime is infinity.  
       ! This will speed up GEOS-Chem simulations. (bmy, 12/13/17)
       !====================================================================
       IF ( Input_Opt%NPassive_Decay > 0 ) THEN

          ! Apply loss rate to passive species with finite lifetimes
          CALL Chem_Passive_Species( am_I_Root,  Input_Opt, State_Chm,       & 
                                     State_Grid, State_Met, RC              )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Passive_Species"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !### Debug
          IF ( LPRT .and. am_I_Root ) THEN
             CALL Debug_Msg( '### MAIN: a CHEMISTRY' )
          ENDIF

       ENDIF

    ENDIF
     
    !=======================================================================
    ! Convert species units back to original unit (ewl, 8/12/15)
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root,  Input_Opt, State_Chm,                &
                            State_Grid, State_Met, OrigUnit,  RC            )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Chemistry timestep [s]
    DT_Chem = Get_Ts_Chem()

#if defined( USE_TEND )
    !=======================================================================
    ! Calculate tendencies and write to diagnostics (ckeller,7/15/2015)
    !=======================================================================

    ! Compute tendencies
    CALL Tend_Stage2( am_I_Root,  Input_Opt, State_Chm,                      &
                      State_Grid, State_Met, 'CHEM',                         &
                      DT_Chem,    RC ) 

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in ""!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

    !----------------------------------------------------------
    ! Chemistry budget diagnostics - Part 2 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetChemistry ) THEN

       ! Get final column masses
       CALL Compute_Column_Mass(                                             &
            am_I_Root   = am_I_Root,                                         &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            SpcMapping  = State_Chm%Map_Advect,                              &
            isFull      = State_Diag%Archive_BudgetChemistryFull,            &
            SpcMapFull  = State_Diag%Map_BudgetChemistryFull,                &
            ColMassFull = State_Diag%BudgetMassFull2,                        &
            isTrop      = State_Diag%Archive_BudgetChemistryTrop,            &
            SpcMapTrop  = State_Diag%Map_BudgetChemistryTrop,                &
            ColMassTrop = State_Diag%BudgetMassTrop2,                        &
            isPBL       = State_Diag%Archive_BudgetChemistryPBL,             &
            SpcMapPBL   = State_Diag%Map_BudgetChemistryPBL,                 & 
            ColMassPBL  = State_Diag%BudgetMassPBL2,                         &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Compute_Column_Mass" (final)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF 
       
       ! Compute chemistry budget diagnostics
       CALL Compute_Budget_Diagnostics(                                      &
            am_I_Root     = am_I_Root,                                       &
            State_Grid    = State_Grid,                                      &
            TS            = DT_Chem,                                         &
            SpcMapping    = State_Chm%Map_Advect,                            &
            isFull        = State_Diag%Archive_BudgetChemistryFull,          &
            SpcMapFull    = State_Diag%Map_BudgetChemistryFull,              &
            diagFull      = State_Diag%BudgetChemistryFull,                  &
            MassInitFull  = State_Diag%BudgetMassFull1,                      &
            MassFinalFull = State_Diag%BudgetMassFull2,                      &
            isTrop        = State_Diag%Archive_BudgetChemistryTrop,          &
            SpcMapTrop    = State_Diag%Map_BudgetChemistryTrop,              &
            diagTrop      = State_Diag%BudgetChemistryTrop,                  &
            MassInitTrop  = State_Diag%BudgetMassTrop1,                      &
            MassFinalTrop = State_Diag%BudgetMassTrop2,                      &
            isPBL         = State_Diag%Archive_BudgetChemistryPBL,           &
            SpcMapPBL     = State_Diag%Map_BudgetChemistryPBL,               &
            diagPBL       = State_Diag%BudgetChemistryPBL,                   &
            MassInitPBL   = State_Diag%BudgetMassPBL1,                       &
            MassFinalPBL  = State_Diag%BudgetMassPBL2,                       &
            RC            = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Compute_Budget_Diagnostics"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE DO_CHEMISTRY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: recompute_od
!
! !DESCRIPTION: Subroutine RECOMPUTE\_OD will update the optical depth values 
!  before accumulating or writing the diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RECOMPUTE_OD( am_I_Root,  Input_Opt,  State_Chm,                &
                           State_Diag, State_Grid, State_Met, RC            )
!
! !USES:
!
    ! References to F90 modules
    USE AEROSOL_MOD,    ONLY : AEROSOL_CONC
    USE AEROSOL_MOD,    ONLY : RDAER
    USE AEROSOL_MOD,    ONLY : SOILDUST
    USE DUST_MOD,       ONLY : RDUST_ONLINE
!    USE DUST_MOD,       ONLY : RDUST_OFFLINE
    USE ErrCode_Mod
    USE ERROR_MOD,      ONLY : Debug_Msg
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : GET_MONTH
    USE TIME_MOD,       ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  03 Fev 2011 - Adapted from chemdr.f by skim
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  13 Nov 2012 - R. Yantosca - Now pass Input_Opt and RC arguments for GIGC
!  15 Nov 2012 - M. Payer    - Now pass all met fields via State_Met
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  12 Aug 2015 - E. Lundgren  - Input tracer units are now [kg/kg] and 
!                               are converted to [kg] for recomputing OD
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: IT_IS_A_FULLCHEM_SIM
    LOGICAL            :: IT_IS_AN_AEROSOL_SIM
    LOGICAL            :: LCARB, LCHEM,  LDUST
    LOGICAL            :: LPRT,  LSSALT, LSULF,      LSOA
    INTEGER            :: MONTH, YEAR,   WAVELENGTH

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! RECOMPUTE_OD begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Recompute_OD  (in module GeosCore/chemistry_mod.F)'

    ! Get month and year
    MONTH                = GET_MONTH()
    YEAR                 = GET_YEAR()

    ! Copy fields from INPUT_OPT to local variables for use below
    LCARB                = Input_Opt%LCARB 
    LCHEM                = Input_Opt%LCHEM
    LDUST                = Input_Opt%LDUST
    LPRT                 = Input_Opt%LPRT
    LSSALT               = Input_Opt%LSSALT
    LSULF                = Input_Opt%LSULF
    LSOA                 = Input_Opt%LSOA
    IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM 

    ! First make sure chemistry is turned on
    IF ( LCHEM ) THEN

       ! Then make sure that the simulations use aerosol species
       IF ( IT_IS_A_FULLCHEM_SIM .or. IT_IS_AN_AEROSOL_SIM ) THEN

          ! And then make sure that the aersol species are defined
          IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

             ! Skip this section if all of these are turned off
             CALL AEROSOL_CONC( am_I_Root,  Input_Opt,  State_Chm,            &
                                State_Diag, State_Grid, State_Met, RC        )

             !==============================================================
             ! Call RDAER -- computes aerosol optical depths
             !==============================================================

             ! Calculate the AOD at the wavelength specified in jv_spec_aod
             WAVELENGTH = 1
             CALL RDAER( am_I_Root,  Input_Opt,  State_Chm,                  &
                         State_Diag, State_Grid, State_Met, RC,              &
                         MONTH,     YEAR,       WAVELENGTH                  )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "RdAer"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             !### Debug
             IF ( LPRT .and. am_I_Root ) THEN 
                CALL Debug_Msg( '### RECOMPUTE_OD: after RDAER' )
             ENDIF

             !==============================================================
             ! If LDUST is turned on, then we have online dust aerosol in
             ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order 
             ! to compute aerosol optical depth for FAST-JX, etc.
             !
             ! If LDUST is turned off, then we don't have online dust 
             ! aerosol in GEOS-CHEM...so read monthly-mean dust files
             ! from disk. (rjp, tdf, bmy, 4/1/04)
             !==============================================================
             IF ( LDUST ) THEN
                CALL RDUST_ONLINE( am_I_Root,  Input_Opt,  State_Chm,       &
                                   State_Diag, State_Grid, State_Met,       &
                                   SOILDUST,   WAVELENGTH, RC              )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Rdust_Online"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

!------------------------------------------------------------------------------
! Prior to 3/3/19:
! Remove RDUST_OFFLINE -- dust should always be on in fullchem and aerosol 
! simulations (mps, 3/3/19)
!#if  !defined( TOMAS )
!             ELSE
!                CALL RDUST_OFFLINE( am_I_Root,  Input_Opt,  State_Chm,      &
!                                    State_Diag, State_Grid, State_Met,      &
!                                    MONTH, YEAR,      WAVELENGTH, RC       )
!
!                ! Trap potential errors
!                IF ( RC /= GC_SUCCESS ) THEN
!                   ErrMsg = 'Error encountered in "Rdust_Offline"!'
!                   CALL GC_Error( ErrMsg, RC, ThisLoc )
!                   RETURN
!                ENDIF
!#endif
!------------------------------------------------------------------------------
             ENDIF

             !### Debug
             IF ( LPRT .and. am_I_Root ) THEN
                CALL DEBUG_MSG( '### RECOMPUTE_OD: after RDUST' )
             ENDIF
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE RECOMPUTE_OD
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_passive_species
!
! !DESCRIPTION: Subroutine CHEM\_PASSIVE\_SPECIES performs loss chemistry 
!  on passive species with finite atmospheric lifetimes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Chem_Passive_Species( am_I_Root,  Input_Opt, State_Chm,         &
                                   State_Grid, State_Met, RC                ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : ind_ 
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Time_Mod,       ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN   ) :: am_I_Root   ! root CPU?
    TYPE(OptInput), INTENT(IN   ) :: Input_Opt   ! Input options object
    TYPE(GrdState), INTENT(IN   ) :: State_Grid  ! Grid state object
    TYPE(MetState), INTENT(IN   ) :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller   - Initial version 
!  03 Nov 2016 - C. Keller   - Moved to chemistry_mod
!  26 Jun 2017 - R. Yantosca - GC_ERROR is now contained in errcode_mod.F90
!  14 Jul 2017 - E. Lundgren - Remove dependency on passive_species_mod.F90
!  02 Aug 2017 - R. Yantosca - Turn off debug print unless ND70 is activated
!  13 Dec 2017 - R. Yantosca - Now apply decay only to those passive species
!                              with finite atmospheric lifetimes
!  04 Jan 2019 - M. Sulprizio- Add capability to specify TAU in half-life;
!                              e-folding time will always be used for now
!  29 Jan 2019 - R. Yantosca - Bug fix: State_Chm should be INTENT(INOUT)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: prtDebug
    LOGICAL             :: Is_HalfLife
    INTEGER             :: I,       J,      L
    INTEGER             :: N,       GCId,   Id
    REAL(fp)            :: DT,      Decay,  Rate

    ! SAVEd scalars
    LOGICAL,  SAVE      :: First = .TRUE.

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg,  ThisLoc
!
! !DEFINED PARAMETERS:
!   
    REAL(fp), PARAMETER :: ln2 = 0.693147181E+00_fp

    !=======================================================================
    ! Chem_Passive_Species begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    prtDebug = ( am_I_Root .and. Input_Opt%LPRT )
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at Chem_Passive_Species (in module GeosCore/chemistry_mod.F)'

    DT       = GET_TS_CHEM() ! timestep in seconds

    ! For now, always compute decay using e-folding time
    Is_HalfLife = .FALSE.

    !=======================================================================
    ! Apply decay loss rate only to those passive species that have a
    ! finite atmospheric lifetime (this speeds up execution)
    !=======================================================================

    ! Loop over all decaying passive species
    DO N = 1, Input_Opt%NPassive_Decay

       !----------------------------------
       ! Find the GEOS-Chem species Id
       !----------------------------------

       ! Get the Id of the species in the passive decay menu
       Id   = Input_Opt%Passive_DecayID(N)

       ! Convert this to a GEOS-Chem species Id number
       GcId = Ind_( TRIM( Input_Opt%PASSIVE_NAME(Id) ) )

       ! Make sure the model ID is valid
       IF ( GcId < 0 ) THEN
          ErrMsg = 'Could not find the GEOS-Chem species ID # '        // &
                   'for passive species : '                            // &
                   TRIM( Input_Opt%PASSIVE_NAME(Id) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !----------------------------------
       ! Compute the decay rate
       !----------------------------------

       ! Compute the decay rate for each passive species
       IF ( Is_HalfLife ) THEN
          Decay = ln2 / Input_Opt%PASSIVE_TAU(Id)
       ELSE
          Decay = 1.0 / Input_Opt%PASSIVE_TAU(Id)
       ENDIF
       Rate  = EXP( - DT * Decay )

       !### Debug output
       IF ( First ) THEN
          IF ( prtDebug ) THEN
             WRITE( 6,100 ) ADJUSTL( Input_Opt%PASSIVE_NAME(Id) ),           &
                            GcId, Rate
 100         FORMAT( '     -  Pass. species name, Id, loss rate: ',&
                      a15, i5, 1x, es13.6 )
          ENDIF
       ENDIF

       !----------------------------------
       ! Apply loss
       !----------------------------------

       !$OMP PARALLEL DO                  &
       !$OMP DEFAULT( SHARED            ) &
       !$OMP PRIVATE( I, J, L           )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Chm%Species(I,J,L,GcId) = State_Chm%Species(I,J,L,GcId)      &
                                        * Rate
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO
 
    ! Reset after the first time
    IF ( First) First = .FALSE.

  END SUBROUTINE Chem_Passive_Species
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_chemistry
!
! !DESCRIPTION: Subroutine INIT\_CHEMISTRY initializes chemistry
! variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Chemistry( am_I_Root,  Input_Opt,  State_Chm,              &
                             State_Diag, State_Grid, RC                     ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE FAST_JX_MOD,    ONLY : Init_FJX
    USE FlexChem_Mod,   ONLY : Init_FlexChem
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)     :: am_I_Root   ! Is this the root CPU?
    TYPE(GrdState), INTENT(IN)     :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(INOUT)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT)  :: State_Diag  ! Diagnostics State object
    INTEGER,        INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  19 May 2014 - C. Keller   - Initial version (stripped from do_chemistry
!                              and chemdr.F)
!  20 Jun 2014 - R. Yantosca - Now pass Input_Opt to INIT_FJX
!  23 Jun 2016 - R. Yantosca - Remove call to SETTRACE, it's obsolete
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! INIT_CHEMISTRY begins here!
    !=======================================================================
    
    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_Chemistry  (in module GeosCore/chemistry_mod.F)'

    ! Skip if we have already done this
    IF ( FIRST ) THEN

       ! Adjust first flag
       FIRST  = .FALSE.

       ! Define species ID's
       id_DST1 = Ind_( 'DST1' )
       id_NK1  = Ind_( 'NK1'  )

       !--------------------------------------------------------------------
       ! Initialize FlexChem
       !--------------------------------------------------------------------
       CALL Init_FlexChem( am_I_Root, Input_Opt, State_Chm, State_Diag, RC  )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_FlexChem"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Initialize Fast-JX photolysis
       !--------------------------------------------------------------------
       CALL Init_FJX( am_I_Root,  Input_Opt, State_Chm, State_Diag, &
                      State_Grid, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_FJX"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE Init_Chemistry
!EOC
END MODULE Chemistry_Mod
