!==============================================================================
!==============================================================================
! Module with named state variables, parameters, ordinary variables, forcings
!==============================================================================
!==============================================================================

MODULE MESOdeclaration

IMPLICIT NONE

!--------------------------- 
! Parameters
!---------------------------
 INTEGER, PARAMETER :: nPars = 71  ! number of parameters

 DOUBLE PRECISION ::                                       &
&     DIC_dC_b , &  !  KARLINE: THIS SHOULD BE A PARAMETER!!!!!!!!!!!!
&     H        , &  !  1.    , # height of the tank [m]
&     Kexc_DIC , &  !  5.0   , # Gas exchange coeff dic [umol C l-1 uatm-1 d-1] VdM04
&     pCO2atm  , &  !400.0   , # partial pressure CO2 in air
&     Kexc_Oxy , &  !  4.8   , # Gas exchange coeff oxygen [±20 cm h-1] Kihm
&     Q10      , &  !  2.    , # temp limitation of process rates [-]
 
&     ON_nit   , &  !  2.0   , # O : N ratio in nitrification
&     nitr_NH4 , &  !  0.03  , # nitrification rate (d-1)
        
&     alpha    , &  !  0.025 , # slope of the P-I curve [umol m-2 s-1 d-1] VdM04
&     K_DIC_P1 , &  !250.0   , # half-saturation constant DIC uptake [mmol C m-3] 
&     K_NH4_P1 , &  !  0.5   , # half-saturation constant NH4 uptake [mmol N m-3] 
&     K_NOX_P1 , &  !  0.5   , # half-saturation constant NOX uptake [mmol N m-3] 
&     psi_P1   , &  !  1.5   , # NH3 inhibition for NOX uptake [m3 mmol N-1]
&     uN_P1    , &  !  2.0   , # net growth rate for NOX and NH4 [mmol N mmol C-1 d-1]
&     uC_P1    , &  !  5.5   , # net growth rate [d-1]
&     minCN_P1 , &  !  6.0   , # minimum CN ratio [mmol C mmol N-1]
&     maxCN_P1 , &  ! 20.0   , # maximum CN ratio [mmol C mmol N-1]
&     gamma1_P1, &  !  0.02  , # leakage fraction [-]
&     gamma3_P1, &  !  0.05  , # DOC produced due to Nlim [-]
&     resb_P1_C, &  !  0.01  , # basal respiration [d-1]
&     resa_P1_C, &  !  0.10  , # activity respiration [-]
&     mort_P1_N, &  !  0.045 , # mortality rate [d-1] 

&     K_DIC_P2 , &  ! 250.0   , # half-saturation constant DIC uptake [mmol C m-3] 
&     K_NH4_P2 , &  !   0.5   , # half-saturation constant NH4 uptake [mmol N m-3] 
&     K_NOX_P2 , &  !   0.5   , # half-saturation constant NOX uptake [mmol N m-3] 
&     psi_P2   , &  !   1.5   , # NH3 inhibition for NOX uptake [m3 mmol N-1]
&     uN_P2    , &  !   2.0   , # net growth rate phyto1 for NOX and NH4 [mmol N mmol C-1 d-1]
&     uC_P2    , &  !   5.5   , # net growth rate phyto1 [d-1]
&     minCN_P2 , &  !   6.0   , # minimum CN ratio [mmol C mmol N-1]
&     maxCN_P2 , &  !  20.0   , # maximum CN ratio [mmol C mmol N-1]
&     gamma1_P2, &  !   0.02  , # leakage fraction [-]
&     gamma3_P2, &  !   0.05  , # DOC produced due to Nlim [-]
&     resb_P2_C, &  !   0.01  , # basal respiration [d-1]
&     resa_P2_C, &  !   0.10  , # activity respiration [-]
&     mort_P2_N, &  !   0.045 , # mortality rate [d-1] 
        
&     delta    , &  !   0.65  , # part of mortality and grazing flux to DMl (1-delta to DMs)
&     epsilon  , &  !   0.34  , # part of mortality and grazing flux to DM (1-epsilon to PM)
        
&     CN_B     , &  !   5.1   , # molar C:N ratio of bacteria
&     uB       , &  !  13.3   , # max uptake of DMl by Bacteria, d-1 {VdM04}
&     K_DMl_B  , &  !  25.0   , # half-saturation uptake of DMl by Bacteria [mmolC m-3], {VdM04}
&     K_NH4_B  , &  !   0.5   , # half-saturation uptake of NH4 by Bacteria [mmolN m-3], {VdM04}
&     omega    , &  !   0.19  , # bacterial growth efficiency [-]
&     mort_B   , &  !   0.05  , # bacterial mortality (quadratic term) [mmol C-1 d-1], {VdM04}
        
&     maxCN_PMw, &  !  15.0   , # max molar CN ratio of PMw [C/N]
&     deg_PMw_C, &  !   0.29  , # degradation rate of C-detritus into DMl + DMs [d-1]
&     deg_PMw_N, &  !   0.33  , # degradation rate of N-detritus into DMl + DMs [d-1]
        
&     uDMs     , &  !   4.0   , # max dissolution rate of DMs_C [d-1]  
&     K_DMs    , &  ! 417.0   , # half-saturation constant for dissolution  of DMs_C [mmolC m-3]
      
&     w        , &  !   0.05  , # settling velocity of PMw [m/d] [default = 0.5]
        
&     maxCN_PMs, &  !  15.0   , # max molar CN ratio of PMs [C/N]
&     deg_PMs_C, &  !   0.29  , # degradation rate of C-detritus into DMl + DMs [d-1]
&     deg_PMs_N, &  !   0.33  , # degradation rate of N-detritus into DMl + DMs [d-1]
      
&     CN_mZP    , &  !  6.25  , # molar CN ratio zooplankton
&     gra_mZP   , &  !  1.0   , # max grazing mZP [d-1]
&     K_mZP     , &  !  1.0   , # half-saturation grazgin by mZP [mmolN m-3]
&     pP1_mZP   , &  !  0.4   , # preference for P1 by mZP [-]
&     pP2_mZP   , &  !  0.3   , # preference for P2 by mZP [-]
&     pB_mZP    , &  !  0.2   , # preference for B by mZP [-] 
&     pPMw_mZP  , &  !  0.1   , # preference for PMw by mZP [-]
&     phi       , &  !  0.23  , # mZP feeding losses to DOM [-]
&     beta_mZP  , &  !  0.75  , # ass eff mZP [-]
&     resb_mZP_C, &  !  0.03  , # basal resp mZP [d-1]
&     resa_mZP_C, &  !  0.25  , # fraction of assimilate food that is respired by mZP [-]
&     mort_mZP  , &  !  0.3   , # zooplankton mortality [d-1]
&     K_mor_mZP , &  !  0.2   , # half-saturation constant for zooplankton mortality [mmolN m-3]

&     CN_Oy     , &  !  6.25  , # CN ratio oysters [C/N}]
&     filt_Oy   , &  !  0.002 , # filtration rate oysters [m^3 mmol C-1 d-1]
&     fae_Oy    , &  !  0.0   , # fraction released as faeces
&     pse_Oy    , &  !  0.65  , # fraction released as pseudeo-faeces
&     resa_Oy_C , &  !  0.25  , # fraction of assimilated food that is respired by Oy [-]
&     resb_Oy_C      ! 0.005    # maint resp Oy [d-1]
  
 COMMON /xcbparms/DIC_dC_b, H, Kexc_DIC, pCO2atm, Kexc_Oxy, Q10, ON_nit, nitr_NH4, &
&   alpha, K_DIC_P1, K_NH4_P1, K_NOX_P1, psi_P1, uN_P1, uC_P1, minCN_P1, maxCN_P1, &
&   gamma1_P1, gamma3_P1, resb_P1_C, resa_P1_C, mort_P1_N, K_DIC_P2, K_NH4_P2,     &
&   K_NOX_P2, psi_P2, uN_P2, uC_P2, minCN_P2, maxCN_P2, gamma1_P2, gamma3_P2,      &
&   resb_P2_C, resa_P2_C, mort_P2_N, delta, epsilon, CN_B, uB, K_DMl_B, K_NH4_B,   &
&   omega, mort_B, maxCN_PMw, deg_PMw_C, deg_PMw_N, uDMs, K_DMs, w, maxCN_PMs,     & 
&   deg_PMs_C, deg_PMs_N, CN_mZP, gra_mZP, K_mZP, pP1_mZP, pP2_mZP, pB_mZP,        &
&   pPMw_mZP, phi, beta_mZP, resb_mZP_C, resa_mZP_C, mort_mZP, K_mor_mZP, CN_Oy,   &
&   filt_Oy, fae_Oy, pse_Oy, resa_Oy_C, resb_Oy_C

!---------------------------
! State variables
!---------------------------

  INTEGER, PARAMETER :: nSvars = 35  ! number of state variables
  
  DOUBLE PRECISION :: Oxy, NOX, NH4, DIC, DIC_13C, P1_N, P1_C, P1_13C,             &
&       P2_N, P2_C, P2_13C, DMl_N, DMl_C, DMl_13C, DMs_N, DMs_C, DMs_13C,          &
&       PMw_N, PMw_C, PMw_13C, PMs_N, PMs_C, PMs_13C, B_N, B_C, B_13C,             &
&       mZP_N, mZP_C, mZP_13C, Oy_N, Oy_C, Oy_13C, MB_N, MB_C, MB_13C

  COMMON/xcbsvars/    Oxy, NOX, NH4, DIC, DIC_13C, P1_N, P1_C, P1_13C,             &
&       P2_N, P2_C, P2_13C, DMl_N, DMl_C, DMl_13C, DMs_N, DMs_C, DMs_13C,          &
&       PMw_N, PMw_C, PMw_13C, PMs_N, PMs_C, PMs_13C, B_N, B_C, B_13C,             &
&       mZP_N, mZP_C, mZP_13C, Oy_N, Oy_C, Oy_13C, MB_N, MB_C, MB_13C

!---------------------------
! Derivatives
!---------------------------

  DOUBLE PRECISION :: dOxy, dNOX, dNH4, dDIC, dDIC_13C, dP1_N, dP1_C, dP1_13C,     &
&      dP2_N, dP2_C, dP2_13C, dDMl_N, dDMl_C, dDMl_13C, dDMs_N, dDMs_C, dDMs_13C,  &
&      dPMw_N, dPMw_C, dPMw_13C, dPMs_N, dPMs_C, dPMs_13C, dB_N, dB_C, dB_13C,     &
&      dmZP_N, dmZP_C, dmZP_13C, dOy_N, dOy_C, dOy_13C, dMB_N, dMB_C, dMB_13C

  COMMON/xcbdsvars/  dOxy, dNOX, dNH4, dDIC, dDIC_13C, dP1_N, dP1_C, dP1_13C,      &
&      dP2_N, dP2_C, dP2_13C, dDMl_N, dDMl_C, dDMl_13C, dDMs_N, dDMs_C, dDMs_13C,  &
&      dPMw_N, dPMw_C, dPMw_13C, dPMs_N, dPMs_C, dPMs_13C, dB_N, dB_C, dB_13C,     &
&      dmZP_N, dmZP_C, dmZP_13C, dOy_N, dOy_C, dOy_13C, dMB_N, dMB_C, dMB_13C

!---------------------------
! Output variables
!---------------------------

  INTEGER, PARAMETER :: noutvars = 27  ! Number of output variables

  DOUBLE PRECISION :: Oxysat, dC_DIC, dC_P1, dC_P2, dC_DMl, dC_DMs,                &
&    dC_PMw, dC_PMs, dC_B, dC_mZP, dC_Oy, DM_C, CO2exc, CO2exc_13C, pCO2,          & 
&    fpar, ftemp, NOXlim_P2, NH4lim_P2, gra_Oy_C, pp_P1_C, pp_P1_N,                &
&     pp_P2_C, pp_P2_N, exc_B_N, exc_mZP_N, NH4_B_N

  COMMON/xcbvars/    Oxysat, dC_DIC, dC_P1, dC_P2, dC_DMl, dC_DMs,                 &
&    dC_PMw, dC_PMs, dC_B, dC_mZP, dC_Oy, DM_C, CO2exc, CO2exc_13C, pCO2,          & 
&    fpar, ftemp, NOXlim_P2, NH4lim_P2, gra_Oy_C, pp_P1_C, pp_P1_N,                &
&     pp_P2_C, pp_P2_N, exc_B_N, exc_mZP_N, NH4_B_N

!---------------------------
! Forcing functions
!---------------------------

  INTEGER, PARAMETER :: nforcs = 7  ! Number of forcing variables

  DOUBLE PRECISION :: par, temp, pH, satOxy, K0, K1, K2
  COMMON/xcbforcs/par, temp, pH, satOxy, K0, K1, K2  

END MODULE MESOdeclaration

!==============================================================================
! Initialisation at start of run (puts pointers to parameters/forcings)
!==============================================================================

SUBROUTINE initpars(odeparms)
!------------------------------------------------------------------------------
! Initialisation of parameter values
!------------------------------------------------------------------------------
USE MESOdeclaration, only : nPars

IMPLICIT NONE
EXTERNAL :: odeparms
DOUBLE PRECISION :: pars(nPars)   ! the same length as all parameters in module
COMMON/xcbparms/ pars

 CALL odeparms(nPars, pars)
 RETURN 

END SUBROUTINE initpars

!==============================================================================

SUBROUTINE initforcs(odeforcs)
!------------------------------------------------------------------------------
! Initialisation of forcing functions
!------------------------------------------------------------------------------
USE MESOdeclaration, only : nForcs

IMPLICIT NONE
EXTERNAL :: odeforcs
DOUBLE PRECISION :: forcs(nForcs)   ! the same length as all forcings in module
COMMON/xcbforcs/ forcs

 CALL odeforcs(nForcs, forcs)
 RETURN 
 
END SUBROUTINE initforcs

!==============================================================================
!==============================================================================
! Dynamic subroutine
!==============================================================================
!==============================================================================

SUBROUTINE MesoModel(neq, t, y, ydot, yout, ip)

!------------------------------------------------------------------------------
! At each time step (t), estimates the derivatives (ydot) and output variables
! (yout) from the state variables (y)
! Uses parameters, state variables, and forcings by their name (in module)
!------------------------------------------------------------------------------

USE MESOdeclaration  ! module with named variables, forcings, parameters
IMPLICIT NONE

INTEGER          :: neq, ip(*)
DOUBLE PRECISION :: t, y(nSvars), ydot(neq)
DOUBLE PRECISION :: yout(nSvars)

DOUBLE PRECISION, EXTERNAL :: FC_dC, DC_FC

! ordinary variables
DOUBLE PRECISION :: F_DIC, F_P1, F_P2, F_DMl, F_DMs, F_PMw, F_PMs, F_B, F_mZP, F_Oy
DOUBLE PRECISION :: Tlim, pro, CO2
DOUBLE PRECISION :: CN_P1, PARlim_P1, NH4lim_P1, NOXlim_P1, DIClim_P1, NH4_P1_N,   &
&  NOX_P1_N, exc_P1_C, reb_P1_C, rea_P1_C, res_P1_C, mor_P1_N, mor_P1_C
  
DOUBLE PRECISION :: CN_P2, PARlim_P2, DIClim_P2, NH4_P2_N, NOX_P2_N, exc_P2_C,     &
&  reb_P2_C, rea_P2_C, res_P2_C, mor_P2_N, mor_P2_C

DOUBLE PRECISION :: gro_B_C, gro_B_N, res_B_C, mor_B_N, mor_B_C, CN_PMw,           &
&  dis_PMw_C, dis_PMw_N, sed_PMw_C, sed_PMw_N, CN_PMs, dis_PMs_C, dis_PMs_N

DOUBLE PRECISION :: gra_mZP_denominator, gP1_mZP_N, gP2_mZP_N, gB_mZP_N,           &       
&  gPMw_mZP_N, gP1_mZP_C, gP2_mZP_C, gB_mZP_C, gPMw_mZP_C, gra_mZP_N, gra_mZP_C,   &
&  ass_mZP_N, ass_mZP_C, reb_mZP_C, rea_mZP_C, cor_resa_mZP, res_mZP_C

DOUBLE PRECISION :: mor_mZP_N, mor_mZP_C, gP1_Oy_C, gP1_Oy_N, gP2_Oy_C, gP2_Oy_N,  &
&  gPMw_Oy_C, gPMw_Oy_N, gmZP_Oy_C, gmZP_Oy_N, gra_Oy_N, ass_Oy_C, ass_Oy_N,       &
&  rea_Oy_C, reb_Oy_C, res_Oy_C, cor_resa_Oy, exc_Oy_N   

DOUBLE PRECISION :: CN_DMs, dis_DMs_C, dis_DMs_N, nit_NH4, Oxyexc, upt_B_C, upt_B_N 

!----------------------------------------------------------
! Check if input from R is consistent with this subroutine
!----------------------------------------------------------

 IF (neq .NE. nSvars) CALL rexit("neq should be = 35")
 IF (ip(1) < nOutvars) CALL rexit("number of output variables too small")
 
!----------------------------------------------------------
 CALL copySvars (y)  ! Copy from long vector y into named state variables

! Temperature limitation function
 Tlim = dexp(dlog(Q10)*((temp-20.d0)/10.d0))
  
! Fraction 13C/C (F_) and delta13C (dC_)
  F_DIC  = DIC_13C / DIC
  dC_DIC = FC_dC(F_DIC)
  F_P1   = P1_13C / P1_C
  dC_P1  = FC_dC(F_P1)
  F_P2   = P2_13C / P2_C
  dC_P2  = FC_dC(F_P2)
  F_DMl  = DMl_13C / DMl_C
  dC_DMl = FC_dC(F_DMl)
  F_DMs  = DMs_13C / DMs_C
  dC_DMs = FC_dC(F_DMs)
  F_PMw  = PMw_13C / PMw_C
  dC_PMw = FC_dC(F_PMw)
  F_PMs  = PMs_13C / PMs_C
  dC_PMs = FC_dC(F_PMs)
  F_B    = B_13C / B_C
  dC_B   = FC_dC(F_B)
  F_mZP  = mZP_13C / mZP_C
  dC_mZP = FC_dC(F_mZP)
  F_Oy   = Oy_13C / Oy_C
  dC_Oy  = FC_dC(F_Oy)

! Oxygen - note: satOxy is a forcing function here
  Oxyexc = Kexc_Oxy * (satOxy - Oxy)

! Saturated percentage oxygen
  Oxysat = oxy/satOxy*100d0
  
! DIN - KARLINE: NEED OXYGEN LIMITATION HERE
  nit_NH4 = nitr_NH4*NH4
  
! DIC, can update with variable Temp, salinity and pH fixed for now
! note: K0, K1, K2 are forcing functions

  pro    = 10.d0**(-pH )
  CO2    = DIC * pro**2 / (pro**2 + k1*pro + k1*k2)
  pCO2   = CO2 / k0
  CO2exc = Kexc_DIC * (pCO2 - pCO2atm)
  
  IF (CO2exc > 0) THEN
    CO2exc_13C = CO2exc*F_DIC
  ELSE
    CO2exc_13C = CO2exc*dC_FC(DIC_dC_b)
  END IF    
  
! Phytoplankton
  CN_P1     = P1_C / P1_N
  PARlim_P1 = (1.d0-exp(-alpha*par/uC_P1))
  NH4lim_P1 = NH4 / (K_NH4_P1 + NH4)
  NOXlim_P1 = (NOX*dexp(-psi_P1*NH4)) / (K_NOX_P1 + NOX)
  DIClim_P1 = DIC / (K_DIC_P1 + DIC)
  NH4_P1_N  = (NH4lim_P1 * uN_P1) * (1.d0 - minCN_P1/CN_P1) * P1_N * Tlim
  NOX_P1_N  = (NOXlim_P1 * uN_P1) * (1.d0 - minCN_P1/CN_P1) * P1_N * Tlim
  pp_P1_N   = NH4_P1_N + NOX_P1_N
  pp_P1_C   = PARlim_P1 * DIClim_P1 * (1.d0 - CN_P1/maxCN_P1) * uC_P1 * P1_C * Tlim
  exc_P1_C  = gamma3_P1 *PARlim_P1 *DIClim_P1 *(CN_P1 - minCN_P1)/maxCN_P1 *uC_P1 *P1_C *Tlim
  reb_P1_C  = resb_P1_C * P1_C * Tlim
  rea_P1_C  = resa_P1_C *  pp_P1_C
  res_P1_C  = reb_P1_C + rea_P1_C
  mor_P1_N  = mort_P1_N * P1_N * Tlim
  mor_P1_C  = mort_P1_N * CN_P1 
  
  CN_P2     = P2_C / P2_N
  PARlim_P2 = (1.d0-exp(-alpha*par/uC_P2))
  NH4lim_P2 = NH4 / (K_NH4_P2 + NH4)
  NOXlim_P2 = (NOX*exp(-psi_P2*NH4)) / (K_NOX_P2 + NOX)
  DIClim_P2 = DIC / (K_DIC_P2 + DIC)
  NH4_P2_N  = (NH4lim_P2 * uN_P2) * (1.d0 - minCN_P2/CN_P2) * P2_N * Tlim
  NOX_P2_N  = (NOXlim_P2 * uN_P2) * (1.d0 - minCN_P2/CN_P2) * P2_N * Tlim
  pp_P2_N   = NH4_P2_N + NOX_P2_N
  pp_P2_C   = PARlim_P2 * DIClim_P2 * (1.d0 - CN_P2/maxCN_P2) * uC_P2 * P2_C * Tlim
  exc_P2_C  = gamma3_P2 *PARlim_P2 *DIClim_P2 *(CN_P2 - minCN_P2)/maxCN_P2 *uC_P2 *P2_C *Tlim
  reb_P2_C  = resb_P2_C * P2_C * Tlim
  rea_P2_C  = resa_P2_C *  pp_P2_C
  res_P2_C  = reb_P2_C + rea_P2_C
  mor_P2_N  = mort_P2_N * P2_N * Tlim
  mor_P2_C  = mort_P2_N * CN_P2 

  ! Dissolved organic matter
  CN_DMs    = DMs_C / DMs_N
  dis_DMs_C = uDMs * B_C * DMs_C / (K_DMs + DMs_C)
  dis_DMs_N = dis_DMs_C / CN_DMs
  
  ! Bacteria
  upt_B_C   = uB * B_C * DMl_C / (K_DMl_B + DMl_C)
  upt_B_N   = upt_B_C * DMl_N / DMl_C
  NH4_B_N   = uB * B_N * NH4 / (K_NH4_B + NH4) ! potential NH4 uptake
  exc_B_N   = upt_B_N - upt_B_C*omega/CN_B 
  gro_B_C   = omega * upt_B_C

  IF (exc_B_N > 0) THEN
    gro_B_N = omega * upt_B_C / CN_B
    res_B_C = (1.d0-omega) * upt_B_C
  ELSE IF (-exc_B_N <= NH4_B_N) THEN
    gro_B_N = omega * upt_B_C / CN_B
    res_B_C = (1.d0-omega) * upt_B_C
  ELSE
    gro_B_N = upt_B_N + NH4_B_N 
    res_B_C = gro_B_C * (1.d0/omega - 1.d0)
    exc_B_N = -NH4_B_N                     ! exc_BPN is negative and == potential uptake
  END IF

  mor_B_N = mort_B  * B_N**2 *Tlim
  mor_B_C = mor_B_N * CN_B
  
! Suspended particulate organic matter
  CN_PMw    = PMw_C / PMw_N
  dis_PMw_C = deg_PMw_C*(1.d0-CN_PMw/maxCN_PMw)*PMw_C*Tlim
  dis_PMw_N = deg_PMw_N*(1.d0-CN_PMw/maxCN_PMw)*PMw_N*Tlim
  sed_PMw_C = w * PMw_C
  sed_PMw_N = w * PMw_N
  
! Sedimented particulate organic matter
  CN_PMs    = PMs_C / PMs_N
  dis_PMs_C = deg_PMs_C*(1.d0-CN_PMs/maxCN_PMs)*PMs_C*Tlim
  dis_PMs_N = deg_PMs_N*(1.d0-CN_PMs/maxCN_PMs)*PMs_N*Tlim
   
 
! Zooplankton
  gra_mZP_denominator = (K_mZP*(pP1_mZP*P1_N + pP2_mZP*P2_N + pB_mZP*B_N + pPMw_mZP*PMw_N) +  &
&                   pP1_mZP*P1_N**2 + pP2_mZP*P2_N**2 + pB_mZP*B_N**2 + pPMw_mZP*PMw_N**2)
  gP1_mZP_N  = gra_mZP*mZP_N*pP1_mZP * P1_N**2 / gra_mZP_denominator * Tlim
  gP2_mZP_N  = gra_mZP*mZP_N*pP2_mZP * P2_N**2 / gra_mZP_denominator * Tlim
  gB_mZP_N   = gra_mZP*mZP_N*pB_mZP  *  B_N**2 / gra_mZP_denominator * Tlim
  gPMw_mZP_N = gra_mZP*mZP_N*pPMw_mZP*PMw_N**2 / gra_mZP_denominator * Tlim
  gP1_mZP_C  = gP1_mZP_N *CN_P1
  gP2_mZP_C  = gP2_mZP_N *CN_P2
  gB_mZP_C   = gB_mZP_N  *CN_B
  gPMw_mZP_C = gPMw_mZP_N*CN_PMw
  gra_mZP_N  = gP1_mZP_N       + gP2_mZP_N       + gB_mZP_N      + gPMw_mZP_N
  gra_mZP_C  = gP1_mZP_C       + gP2_mZP_C       + gB_mZP_C      + gPMw_mZP_C
  ass_mZP_N  = (1-phi)*beta_mZP*gra_mZP_N
  ass_mZP_C  = (1-phi)*beta_mZP*gra_mZP_C
  reb_mZP_C  = resb_mZP_C*mZP_C*Tlim
  rea_mZP_C  = resa_mZP_C*ass_mZP_C
  cor_resa_mZP = 1
  res_mZP_C  = reb_mZP_C + rea_mZP_C
  exc_mZP_N  = ass_mZP_N - (ass_mZP_C - res_mZP_C)/CN_mZP

  IF (exc_mZP_N  <= 0) THEN
    exc_mZP_N = 0 
    cor_resa_mZP = (ass_mZP_C - ass_mZP_N*CN_mZP - reb_mZP_C) / rea_mZP_C 
    rea_mZP_C    = cor_resa_mZP * rea_mZP_C
    
    IF (rea_mZP_C < 0.d0) THEN
        rea_mZP_C = 0.d0 
        reb_mZP_C = ass_mZP_C - ass_mZP_N*CN_mZP
    END IF
    
    res_mZP_C  = reb_mZP_C + rea_mZP_C
  END IF
  
  mor_mZP_N   = mort_mZP*mZP_N**2 / (K_mor_mZP + mZP_N) * Tlim
  mor_mZP_C   = mor_mZP_N*CN_mZP
  
! Oysters
  gP1_Oy_C   = filt_Oy * P1_C  * Oy_C * Tlim
  gP1_Oy_N   = gP1_Oy_C / CN_P1             
  gP2_Oy_C   = filt_Oy * P2_C  * Oy_C * Tlim
  gP2_Oy_N   = gP2_Oy_C / CN_P2             
  gPMw_Oy_C  = filt_Oy * PMw_C * Oy_C * Tlim
  gPMw_Oy_N  = gPMw_Oy_C / CN_PMw           
  gmZP_Oy_C  = filt_Oy * mZP_C * Oy_C * Tlim
  gmZP_Oy_N  = gmZP_Oy_C / CN_mZP           
  gra_Oy_C   = gP1_Oy_C + gP2_Oy_C + gPMw_Oy_C + gmZP_Oy_C
  gra_Oy_N   = gP1_Oy_N + gP2_Oy_N + gPMw_Oy_N + gmZP_Oy_N
  ass_Oy_C   = (1.d0-fae_Oy-pse_Oy) * gra_Oy_C
  ass_Oy_N   = (1.d0-fae_Oy-pse_Oy) * gra_Oy_N
  rea_Oy_C   = resa_Oy_C * ass_Oy_C
  reb_Oy_C   = resb_Oy_C * Oy_C * Tlim      
  res_Oy_C   = rea_Oy_C + reb_Oy_C
  cor_resa_Oy= 1.d0
  exc_Oy_N   = ass_Oy_N - (ass_Oy_C - res_Oy_C)/CN_Oy

  IF (exc_Oy_N  <= 0) THEN
    exc_Oy_N = 0.d0 
    cor_resa_Oy = (ass_Oy_C - ass_Oy_N*CN_Oy - reb_Oy_C) / rea_Oy_C 
    rea_Oy_C    = cor_resa_Oy * rea_Oy_C
    
    IF (rea_Oy_C < 0.d0) THEN
        rea_Oy_C = 0.d0 
        reb_Oy_C = ass_Oy_C - ass_Oy_N*CN_Oy
    END IF
    
    res_Oy_C  = reb_Oy_C + rea_Oy_C
  END IF

!  #####################
!  ## derivatives ##
!  #####################
  dOxy     =  + pp_P1_C        +  pp_P2_C + Oxyexc/H                                   &
&            - ON_nit*nit_NH4 - res_P1_C - res_P2_C - res_B_C - res_mZP_C - res_Oy_C/H

  dNOX     = +nit_NH4  -  NOX_P1_N - NOX_P2_N    
  dNH4     = +exc_B_N  + exc_mZP_N + exc_Oy_N/H                                        &
&                        -NH4_P1_N - NH4_P2_N  - nit_NH4

    dDIC     = +res_P1_C -  pp_P1_C  - exc_P1_C                                        &
&                      +res_P2_C -  pp_P2_C  - exc_P2_C                                &
&                      +res_B_C  + res_mZP_C + res_Oy_C/H                              &
&                    -CO2exc
                
    dDIC_13C = +reb_P1_C *F_P1  + rea_P1_C*F_DIC -  pp_P1_C*F_DIC - exc_P1_C*F_DIC     &
&              +reb_P2_C *F_P2  + rea_P2_C*F_DIC -  pp_P2_C*F_DIC - exc_P2_C*F_DIC     &
&              +res_B_C  *F_DMl                                                        &   
&              +(1.d0-phi)*beta_mZP*resa_mZP_C*cor_resa_mZP*gP1_mZP_C *F_P1            &
&              +(1.d0-phi)*beta_mZP*resa_mZP_C*cor_resa_mZP*gP2_mZP_C *F_P2            &
&              +(1.d0-phi)*beta_mZP*resa_mZP_C*cor_resa_mZP*gB_mZP_C  *F_B             &
&              +(1.d0-phi)*beta_mZP*resa_mZP_C*cor_resa_mZP*gPMw_mZP_C*F_PMw           &
&                +reb_mZP_C*F_mZP                                                      &
&              +(1.d0-fae_Oy-pse_Oy)*resa_Oy_C*cor_resa_Oy*gP1_Oy_C   *F_P1  / H       &
&              +(1.d0-fae_Oy-pse_Oy)*resa_Oy_C*cor_resa_Oy*gP2_Oy_C   *F_P2  / H       &
&              +(1.d0-fae_Oy-pse_Oy)*resa_Oy_C*cor_resa_Oy*gPMw_Oy_C  *F_PMw / H       &
&              +(1.d0-fae_Oy-pse_Oy)*resa_Oy_C*cor_resa_Oy*gmZP_Oy_C  *F_mZP / H       &
&              +reb_Oy_C*F_Oy / H                                                      &
&              -CO2exc_13C
              
    dP1_N    = (1.d0-gamma1_P1)*pp_P1_N       - mor_P1_N             - gP1_mZP_N - gP1_Oy_N/H
    dP1_C    = (1.d0-gamma1_P1)*pp_P1_C       - mor_P1_C  - res_P1_C - gP1_mZP_C - gP1_Oy_C/H
    dP1_13C  = (1.d0-gamma1_P1)*pp_P1_C*F_DIC - mor_P1_C*F_P1                          &
&                 -reb_P1_C    *F_P1          - rea_P1_C*F_DIC                         &
&                 -gP1_mZP_C   *F_P1          - gP1_Oy_C*F_P1/H
                              
    dP2_N    = (1.d0-gamma1_P2)*pp_P2_N       - mor_P2_N             - gP2_mZP_N - gP2_Oy_N/H
    dP2_C    = (1.d0-gamma1_P2)*pp_P2_C       - mor_P2_C  - res_P2_C - gP2_mZP_C - gP2_Oy_C/H
    dP2_13C  = (1.d0-gamma1_P2)*pp_P2_C*F_DIC - mor_P2_C*F_P2                          &
&                             -reb_P2_C    *F_P2          - rea_P2_C*F_DIC             &
&                             -gP2_mZP_C   *F_P2          - gP2_Oy_C*F_P2/H
                   
    dDMl_N   =   gamma1_P1*pp_P1_N           + gamma1_P2*pp_P2_N                       &
&                            +delta*epsilon*mor_P1_N      + delta*epsilon*mor_P2_N     &
&                            -upt_B_N                     + delta*mor_B_N              &
&                            +delta*dis_PMw_N                                          &
&                            +delta*dis_PMs_N   / H                                    &
&                            +dis_DMs_N                                                &
&                            +phi*delta*gra_mZP_N
                        
    dDMl_C   =   gamma1_P1*pp_P1_C           + gamma1_P2*pp_P2_C                       &
&                            +exc_P1_C                    + exc_P2_C                   &
&                            +delta*epsilon*mor_P1_C      + delta*epsilon*mor_P2_C     &
&                            -upt_B_C                       + delta*mor_B_C            &
&                            +delta*dis_PMw_C                                          &
&                            +delta*dis_PMs_C   / H                                    &
&                            +dis_DMs_C                                                &   
&                            +phi*delta*gra_mZP_C   
                               
    dDMl_13C =  gamma1_P1*pp_P1_C*F_DIC     + gamma1_P2*pp_P2_C*F_DIC                  &
&                            +exc_P1_C*F_DIC              + exc_P2_C*F_DIC             &
&                +delta*epsilon*mor_P1_C*F_P1 + delta*epsilon*mor_P2_C*F_P2            &
&                            -upt_B_C*F_DMl               + delta*mor_B_C*F_B          &
&                            +delta*dis_PMw_C*F_PMw                                    &
&                            +delta*dis_PMs_C*F_PMs / H                                &
&                            +dis_DMs_C*F_DMs                                          &
&                            +phi*delta*gP1_mZP_C *F_P1                                &
&                            +phi*delta*gP2_mZP_C *F_P2                                &
&                            +phi*delta*gB_mZP_C  *F_B                                 &
&                            +phi*delta*gPMw_mZP_C*F_PMw        
                            
    
    dDMs_N   =  (1.d0-delta)*epsilon*mor_P1_N      + (1.d0-delta)*epsilon*mor_P2_N     &
&                            +(1.d0-delta)*mor_B_N                                     &
&                            +(1.d0-delta)*dis_PMw_N                                   &
&                            +(1.d0-delta)*dis_PMs_N    / H                            &
&                            -dis_DMs_N                                                &
&                            +phi*(1.d0-delta)*gra_mZP_N    
                
    dDMs_C   =  (1.d0-delta)*epsilon*mor_P1_C      + (1.d0-delta)*epsilon*mor_P2_C     &
&                +(1.d0-delta)*mor_B_C                                                 &
&                            +(1.d0-delta)*dis_PMw_C                                   &
&                            +(1.d0-delta)*dis_PMs_C    / H                            &  
&                            -dis_DMs_C                                                &
&                            +phi*(1.d0-delta)*gra_mZP_C        
                
    dDMs_13C = (1.d0-delta)*epsilon*mor_P1_C*F_P1 +(1.d0-delta)*epsilon*mor_P2_C*F_P2  &  
&                            +(1.d0-delta)*mor_B_C*F_B                                 &
&                            +(1.d0-delta)*dis_PMw_C*F_PMw                             &
&                            +(1.d0-delta)*dis_PMs_C*F_PMs  / H                        &
&                            -dis_DMs_C*F_DMs                                          &
&                            +phi*(1.d0-delta)*gP1_mZP_C *F_P1                         &
&                            +phi*(1.d0-delta)*gP2_mZP_C *F_P2                         &
&                            +phi*(1.d0-delta)*gB_mZP_C  *F_B                          &
&                            +phi*(1.d0-delta)*gPMw_mZP_C*F_PMw     

    dPMw_N   =  (1.d0-epsilon)*mor_P1_N       + (1.d0-epsilon)*mor_P2_N                &
&                            - dis_PMw_N                                               &
&                            - gPMw_mZP_N                  + mor_mZP_N                 &
&                            - gPMw_Oy_N/H                                             &
&                            - sed_PMw_N/H                                             &
&                            +(1.d0-phi)*(1.d0-beta_mZP)*gra_mZP_N  
              
    dPMw_C   =  (1.d0-epsilon)*mor_P1_C       + (1.d0-epsilon)*mor_P2_C                &
&                            - dis_PMw_C                                               &
&                            - gPMw_mZP_C                   + mor_mZP_C                &
&                            - gPMw_Oy_C/H                                             &
&                            - sed_PMw_C/H                                             &
&                            +(1.d0-phi)*(1-beta_mZP)*gra_mZP_C     
                
    dPMw_13C =  (1.d0-epsilon)*mor_P1_C*F_P1  + (1.d0-epsilon)*mor_P2_C*F_P2           &
&                            - dis_PMw_C *F_PMw                                        &  
&                            - gPMw_mZP_C*F_PMw             + mor_mZP_C*F_mZP          &
&              - gPMw_Oy_C *F_PMw/H                                                    &
&              - sed_PMw_C *F_PMw/H                                                    &
&                            +(1.d0-phi)*(1.d0-beta_mZP)*gP1_mZP_C *F_P1               &
&                            +(1.d0-phi)*(1.d0-beta_mZP)*gP2_mZP_C *F_P2               &
&                            +(1.d0-phi)*(1.d0-beta_mZP)*gB_mZP_C  *F_B                &
&                            +(1.d0-phi)*(1.d0-beta_mZP)*gPMw_mZP_C*F_PMw       

    dPMs_N   =   (fae_Oy+pse_Oy) * gra_Oy_N          + sed_PMw_N  - dis_PMs_N
    dPMs_C   =   (fae_Oy+pse_Oy) * gra_Oy_C          + sed_PMw_C  - dis_PMs_C     
    dPMs_13C =   (fae_Oy+pse_Oy) * gP1_Oy_C *F_P1                                      &
&                            +(fae_Oy+pse_Oy) * gP2_Oy_C *F_P2                         &
&                            +(fae_Oy+pse_Oy) * gPMw_Oy_C*F_PMw                        &
&                            +(fae_Oy+pse_Oy) * gmZP_Oy_C*F_mZP                        &
&                               +sed_PMw_C      * F_PMw                                &
&                               -dis_PMs_C      * F_PMs

    dB_N     =   gro_B_N       - mor_B_N     - gB_mZP_N
    dB_C     =   gro_B_C       - mor_B_C     - gB_mZP_C
    dB_13C   =   gro_B_C*F_DMl - mor_B_C*F_B - gB_mZP_C*F_B
    
    dmZP_N   =   ass_mZP_N     - gmZP_Oy_N/H - mor_mZP_N   - exc_mZP_N
    dmZP_C   =   ass_mZP_C     - gmZP_Oy_C/H - mor_mZP_C   - res_mZP_C
    dmZP_13C =  (1.d0-phi)*beta_mZP*(1-resa_mZP_C*cor_resa_mZP)*gP1_mZP_C *F_P1        &
&              +(1.d0-phi)*beta_mZP*(1-resa_mZP_C*cor_resa_mZP)*gP2_mZP_C *F_P2        &
&              +(1.d0-phi)*beta_mZP*(1-resa_mZP_C*cor_resa_mZP)*gB_mZP_C  *F_B         &
&              +(1.d0-phi)*beta_mZP*(1-resa_mZP_C*cor_resa_mZP)*gPMw_mZP_C*F_PMw       &
&              -gmZP_Oy_C*F_mZP/H                                                      & 
&              -mor_mZP_C*F_mZP                                                        &
&              -reb_mZP_C*F_mZP
                
    
    dOy_N    =  ass_Oy_N  - exc_Oy_N
    dOy_C    =  ass_Oy_C  - res_Oy_C
    
    dOy_13C  =  (1.d0-fae_Oy-pse_Oy)*(1.d0-resa_Oy_C*cor_resa_Oy)*gP1_Oy_C   *F_P1     &
&              +(1.d0-fae_Oy-pse_Oy)*(1.d0-resa_Oy_C*cor_resa_Oy)*gP2_Oy_C   *F_P2     &
&              +(1.d0-fae_Oy-pse_Oy)*(1.d0-resa_Oy_C*cor_resa_Oy)*gPMw_Oy_C  *F_PMw    &
&              +(1.d0-fae_Oy-pse_Oy)*(1.d0-resa_Oy_C*cor_resa_Oy)*gmZP_Oy_C  *F_mZP    &
&              -reb_Oy_C*F_Oy
                
    
! Mass balances
    dMB_N    =  dNOX     + dNH4+ dP1_N   + dP2_N    + dDMl_N   + dDMs_N  + dPMw_N      &
&             + dPMs_N  /H + dB_N   + dmZP_N    + dOy_N  /H             
    dMB_C    =  dDIC           + dP1_C   + dP2_C    + dDMl_C   + dDMs_C  + dPMw_C      &
&               + dPMs_C  /H + dB_C   + dmZP_C    + dOy_C  /H + CO2exc     
    dMB_13C  = dDIC_13C        + dP1_13C + dP2_13C  + dDMl_13C + dDMs_13C + dPMw_13C   &
&               + dPMs_13C/H + dB_13C + dmZP_13C  + dOy_13C/H + CO2exc_13C
 
 
! Output variables

  DM_C = DMl_C+DMs_C
  
  fpar = par
  ftemp = temp

  CALL copydSvars (ydot)
  CALL copyOutvars(yout)

END SUBROUTINE MesoModel

!==============================================================================
! Two functions to esimate stable isotope delta values
!==============================================================================

DOUBLE PRECISION FUNCTION FC_dC (FC)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: FC
DOUBLE PRECISION, PARAMETER :: Cstd = 0.0112372d0
DOUBLE PRECISION :: FC_RC, RC_DC

  FC_RC = FC/(1.d0-FC)
  
  FC_dC = (FC_RC / Cstd - 1.d0)*1d3

END FUNCTION FC_dC

!==============================================================================

DOUBLE PRECISION FUNCTION DC_FC (DC)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: DC
DOUBLE PRECISION, PARAMETER :: Cstd = 0.0112372d0
DOUBLE PRECISION :: DC_RC

  DC_RC = ((DC/1000.D0)+1.D0)*Cstd

  DC_FC = DC_RC / (DC_RC + 1.D0)

END FUNCTION DC_FC

!==============================================================================
! Helper functions to copy from module to R-passed variables
!==============================================================================

SUBROUTINE copySvars(y)
!------------------------------------------------------------------------------
! Copy long vector of state variables passed from solver in common block svars
!------------------------------------------------------------------------------
USE MESOdeclaration, only : nSvars
IMPLICIT NONE

DOUBLE PRECISION :: svars(nSvars), y(nSvars)
COMMON/xcbsvars/ svars
INTEGER :: I
 
 DO I = 1, nSvars
  svars(i) = y(i)
 END DO
 
END SUBROUTINE copySvars   

!==============================================================================

SUBROUTINE copydSvars(dy)
!------------------------------------------------------------------------------
! Copy derivatives from common block to long vector used in solver
!------------------------------------------------------------------------------

USE MESOdeclaration, only : nSvars
IMPLICIT NONE
DOUBLE PRECISION :: dsvars(nSvars), dy(nSvars)
COMMON/xcbdsvars/ dsvars
INTEGER :: I
 
 DO I = 1, nSvars
  dy(i) = dsvars(i)
 END DO
 
END SUBROUTINE copydSvars   

!==============================================================================

SUBROUTINE copyOutvars(rpar)

!------------------------------------------------------------------------------
! Saves all output variables in vector rpar
! rpar will be passed back to the solver and returned to R
! The common block with the output variables is defined as a long vector here
!------------------------------------------------------------------------------

USE MESOdeclaration, only : nOutvars
IMPLICIT NONE
DOUBLE PRECISION :: rpar(*), vars(nOutvars)
COMMON/xcbvars/vars
 
INTEGER :: I
 
 DO I = 1, nOutvars
  rpar(i) = vars(i)
 END DO

END SUBROUTINE copyOutvars
