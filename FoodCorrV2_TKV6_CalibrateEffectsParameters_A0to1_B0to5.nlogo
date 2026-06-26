; Script-file -- Copper DEB-IBM Lymnaea stagnalis

; Author: Karel Vlaeminck
; Adaptation of the extended script from Karel Viaene & Andreas Focks
; Adapted from original script-file by Ben Martin (btmarti25@gmail.com) and Elke Zimmer
; Implementation of the standard DEB equations in an IBM



; =========================================================================================================================================
; ========================== DEFINITION OF PARAMETERS AND STATE VARIABLES =================================================================
; =========================================================================================================================================


; - - - - - - - - - - - - - - - - Global parameters: are accessible for patches and turtles - - - - - - - - - - - - - - - - - - - - - - - -


globals[

  L_0                    ; cm, initial structural volume
  day

  ; Embryo related parameters for calc-embryo-reserve-investment

  embryo-timestep

  e_scaled_embryo
  e_ref
  S_C_embryo
  U_H_embryo
  dU_E_embryo
  dU_H_embryo
  dL_embryo

  ; Parameters used to calculate the costs for an egg / initial reserves

  lower-bound ; lower boundary for shooting method
  upper-bound ; upper boundary for shooting method
  sim         ; this keeps track of how many times the calc-egg-size loop is r


  ; General parameters for processes

  T_arrh                 ; arrhenius temperature for species
  T_Fact                 ; arrhenius factor
  K                      ; half saturation food density
  temperature            ; the temperature of the environment [K]

  ; General global parameters

  id                     ; this is used to create a unique id for each simulation
  species_traits_matrix  ; species file
  temperature-matrix     ; temperate file
  species_clear_name     ; global variable, gives the clear (latin) name of the simulated species
  meanDispersal          ; global constant, gives the average step size per day

  starv_fraction         ; fraction of the maximum length attained that an individual can shrink to before dying

  ; DEB parameters

  p_M                    ; volume-specific somatic maintenance
  E_H^b                  ; maturity at birth
  E_H^j                  ; maturity at transformation
  E_H^p                  ; maturity at puberty
  E_G                    ; specific cost for structure
  v_rate_int             ; energy conductance
  kap_int                ; allocation fraction to soma (structural growth + structural maintenance)
  kap_R_int              ; reproduction efficiency
  k_J_rate_int           ; maturity maintenance rate coefficient
  p_Xm                   ; maximum surface specific ingestion rate
  F_m                    ; maximum surface specific searching rate
  time-between-repro     ; time between reproduction
  sG                     ; Gompertz stress coefficient
  h_a                    ; Weibull aging acceleration
  egg_weight             ; estimated energy of an embryo
  L_embryo               ; estimated intial length of an embryo
  U_E_embryo             ; estimated initial energy of an embryo

  zoom

  ; Convert of Add-my-pet parameters to standard (scaled) DEB parameters

  p_Am                   ; maximum surface specific assimilation rate (also equal to y_EX * p_Xm = 0.8 * 49.275 = 39.42)
  U_H^b_int              ; scaled maturity at birth
  U_H^p_int              ; scaled maturity at puberty
  U_H^j_int              ; scaled maturity at transformation

  ; Compound parameters

  k_M_rate_int           ; specific somatic maintenance rate
  g_int                  ; energy investment ration

  ; Copper mortality parameters

  mort_rate_juv          ; mortality rate acute toxicity metal
  mort_rate_adult
  mort_juv               ; chance of dying due to acute toxicity metal
  mort_adult

  ; Copper stress factors

  stress                 ; the parameter s used in the various PMoA stress effects
  ;f_stress               ; stress on the feeding factor f
  ;g_stress               ; stress on the parameter g
  k_M_stress1            ; stress on the maintenance rate > maintenance costs increase
  k_M_stress2            ; stress on the maintenance rate > growth cost increase
  k_J_stress             ; stress on k_J_rate
  repro_stress1          ; stress on the cost for reproduction
  repro_stress2          ; stress on survival of embryos

  ;LC50
  ;slope
  ;time

  A_eq
  H_eq

  iter_N                 ; number of iterations during multiple

  resize                 ; parameter used to resize size distribution
  egg_cost
  cum_repro              ; cummulative reproduction (biomass production)
  age_1strepro           ; age of first reproduction


  ; Food related parameter

  food-factor            ; factor on assimilation (effect fish flakes vs lettuce)

  ; Toxicokinetics parameters- V3 Kristi Weighman
  k_u                    ; uptake rate
  k_e                    ; elminiation rate


]



; - - - - - - - - - - - - - - - - Parameters for the environment: are accessible for patches - - - - - - - - - - - - - - - - - - - - - - - -


patches-own[

  p_food                ; # / cm^2, prey density ; CONSIDER UNITS!!!! > volume/surface you choose reflects the output
  p_d_food              ; change of food density in time

  p_feedingpotential    ; potential feeding by the available grazers
  p_actualfeeding       ; actual feeding taking into account the phyto density
  p_counts              ; amount of bugs at this patch

  p_slice_n             ; used in communication with JAVA for the spatial distribution of patches
  p_area                ; used in communication with JAVA for the spatial distribution of patches
  p_depth

]



; - - - - - - - - - - - - - - - - Parameters for the individuals: are accessible for turtles - - - - - - - - - - - - - - - - - - - - - - - -


turtles-own[

  ; State variables

  t_L                   ; cm, structural length
  t_dL                  ; change of structural length in time
  t_U_H                 ; t L^2, scaled maturity
  t_dU_H                ; change of scaled maturity in time
  t_U_E                 ; t L^2, scaled reserves
  t_dU_E                ; change of scaled reserves in time
  t_e_scaled            ; scaled reserves per unit of structure
  t_U_R                 ; t L^2, scaled energy in reproduction buffer (not standard DEB)
  t_dU_R                ; change of energy in reproduction buffer (reproduction rate)

  ; Fluxes

  t_S_A                 ; assimilation flux
  t_S_C                 ; mobilisation flux

  ; Standard DEB parameters

  t_g                   ; - , energy investment ratio
  t_v_rate              ; cm /t , energy conductance (velocity)
  t_kap                 ; - , allocation fraction to soma
  t_kap_R               ; - , reproduction efficiency
  t_k_M_rate            ; 1/t, somatic maintenance rate coefficient
  t_k_J_rate            ; 1/t, maturity maintenance rate coefficient
  t_U_H^b               ; t L^2, scaled maturity at birth
  t_U_H^p               ; t L^2, scaled maturity at puberty
  t_scatter-multiplier  ; parameter that is used to randomize the input parameters
  t_U_H^j               ; t L^2, scaled maturity at metamorphosis

  ; Prey dynamics

  t_p_Xm_rate           ; J / (cm^2 t), surface-area-specific maximum ingestion rate
  t_p_Am_rate           ; J / (cm^2 t), surface-area-specific maximum assimilation rate (equal to y_EX * t_p_Xm_rate = 0.8 * t_p_Xm_rate)
  t_K                   ; individual half saturation coefficient
  t_functresp           ; - , scaled functional response
  ;t_food-factor

  ; Ageing

  t_q_acceleration      ; - , ageing acceleration
  t_dq_acceleration     ; change of ageing acceleration in time
  t_h_rate              ; - , hazard rate
  t_dh_rate             ; change of hazard rate in time

  ; Animal specific parameters

  t_repro-time          ; time after a reproduction event (delay parameter before new reproduction event)
  t_die?                ; if 1, the individual will die
  t_juvenile            ; if 1, the individual has reached puberty
  t_adult               ; if 1, the individual is an adult
  t_offspring-number    ; number of offsprings
  t_develop-time        ; development time, time after hatching
  t_Lm                  ; individual ultimate structural length

  ; Metabolic acceleration

  t_M                   ; metabolic acceleration coefficient
  t_r                   ; volumetric growth rate (Zimmer 2013)
  t_i
  t_L^j                 ; structural length at transformation
  t_j
  t_L^p                 ; structural length at puberty

  t_XL
  t_X

  t_bmass ;wet weight Kristi

  ; Stress on PMoAs (moved from globals to turtle-specific)- V5 Kristi

  t_stress                 ; the parameter s used in the various PMoA stress effects
  t_f_stress               ; stress on the feeding factor f
  t_g_stress               ; stress on the parameter g
  t_k_M_stress1            ; stress on the maintenance rate > maintenance costs increase
  t_k_M_stress2            ; stress on the maintenance rate > growth cost increase
  t_k_J_stress             ; stress on k_J_rate
  t_repro_stress1          ; stress on the cost for reproduction
  t_repro_stress2          ; stress on survival of embryos


  ; Toxicokinetics

  t_cV                   ; scaled internal concentration
  t_dcV                  ; change in scaled internal concentration

]



; =========================================================================================================================================
; ============================================ SETUP PROCEDURE: SETTING INITIAL CONDITIONS ================================================
; =========================================================================================================================================


; - - - - - - - - - - - - - - - - - - - - - - - - - - Setup button - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to setup
  __clear-all-and-reset-ticks

  set id random-float 3
  set L_0 .00001        ; set initial length to some very small value (embryo start off as nearly all reserves)

  if Temperature-dependency = "from file"
  [
    set temperature-matrix read-in-matrix temperature-file        ; load file with species traits
    set temperature (item 0 temperature-matrix) + 273
  ]

  if Temperature-dependency = "constant" [ set temperature temperature-constant ]

  init_globals
  food-corr
  ;define_TK ; calibrate uptake and elimination Kristi

  ask patches [
    set p_food K_food       ; set initial value of prey to their carrying capacity
    set p_area 1
    set p_depth 1
    set p_slice_n 1
  ]

  calc-embryo-reserve-investment

  ask patches [ sprout initial_number_bugs ]
  ask turtles  [ individual-variability ]  ; first their individual variability in the parameter is set, then the initial energy is calculated for each

  set temperature 293
  set T_fact 1
  T_corrections
  copper_toxicity

  write_Data_output
end



; - - - - - - - - - - - - - - - - - - - - - - - Initialize global parameters - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to init_globals

  set species_clear_name       "Lymnaea"

  ; Load Add-my-pet parameters
  set meanDispersal            1                                          ; m / d, mean dispersal rate
  set p_M                      157.3                                      ; J / d cm³, volume-specific somatic maintenance
  set zoom                     0.1951                                     ; -, zoom factor
  set E_H^b                    0.3417                                     ; J, maturity at birth
  set E_H^j                    217.3                                      ; J, maturity at transformation
  set E_H^p                    721.7                                      ; J, maturity at puberty
  set E_G                      2800                                       ; J / cm³, specific cost for structure
  set v_rate_int               0.02161                                    ; cm / d, energy conductance
  set kap_int                  0.7785                                     ; -, allocation fraction to soma (structural growth + structural maintenance)
  set kap_R_int                0.5                                        ; -, reproduction efficiency
  set k_J_rate_int             0.03804                                    ; d-1, maturity maintenance rate coefficient
  set p_Xm                     49.275                                     ; J / d cm², maximum surface specific ingestion rate
  set F_m                      6.5                                        ; L / d cm², maximum surface specific searching rate
  set time-between-repro       0                                          ; d, time between reproduction
  set sG                       0.0001                                     ; -, Gompertz stress coefficient
  set h_a                      1.287 * (10 ^ -5)                          ; d-2, Weibull aging acceleration
  set egg_weight               0.064941406                                ; J, estimated energy of an embryo
  set L_embryo                 0.073492959                                ; cm, estimated intial length of an embryo
  set U_E_embryo               0.018512204                                ; J, estimat initial energy of an embryo
  set T_arrh                   8000                                       ; K, Arrhenius temperature

  set egg_cost 0.00779 ; cm²d, egg cost in adult-scaled units (= E0_J / p_Am_adult), E0_J = egg_weight × p_Am_base = 0.064941 × 39.42 = 2.56 J, p_Am_adult = t_M × p_Am_base = 8.34 × 39.42 = 328.74 J/d/cm²

  ; Convert of Add-my-pet parameters to standard (scaled) DEB parameters

  set p_Am p_M * zoom / kap_int                                           ; maximum surface specific assimilation rate (also equal to y_EX * p_Xm = 0.8 * 49.275 = 39.42)
  set U_H^b_int E_H^b / p_Am                                              ; scaled maturity at birth
  set U_H^p_int E_H^p / p_Am                                              ; scaled maturity at puberty
  set U_H^j_int E_H^j / p_Am                                              ; scaled maturity at transformation

  ; Compound parameters

  set k_M_rate_int p_M / E_G                                              ; specific somatic maintenance rate
  set g_int (E_G * v_rate_int / p_Am) / kap_int                           ; energy investment ration

  ; Half saturation constant of the population

  set K p_Xm / F_m

  ; Toxicokinetics (elimination rate) - V3 Kristi
   set k_e 0.0118111 ; Dissolved nickel temp-corrected rate for 20 degrees
   set k_u 0.000131304 ; Dissolved nickel temp-corrected rate for 20 degrees
   ;set k_u 0.000673 ; temp-corrected rate for 20 degrees free ion

end



; - - - - - - - - - - - - - - - - - - - Correction for food through the food factor - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; Assimilation rates are corrected through the food factor, which depends on the available food (lettuce, fish flakes, periphyton)


to food-corr
  set p_Xm                     49.275                                     ; J / d cm², maximum surface specific ingestion rate
  set F_m                      6.5                                        ; L / d cm², maximum surface specific searching rate

  if food-type = "calibrate"
  [ set food-factor f_min + random-float (f_max - f_min)] ; randomize food-factor (calibration)

  if food-type = "Fish flakes"
  [ set food-factor 4.025 ]                                               ; value calibrated

  if food-type = "Lettuce"
  [ set food-factor 1 ]
  ;[ set food-factor 1.15 ]
  ;[ set food-factor 1.288 ]                                               ; value calibrated from juvenile control lettuce

  if food-type = "Periphyton"
  [ set food-factor 1 ]

  if food-type = "Brix"
  [ set food-factor 3.753361 ]

  if food-type = "Mattson"
  ;[ set food-factor 2.438 ]
  ;[ set food-factor 2.01 ] ; 18 Degrees
  [ set food-factor 3.14 ] ; 26 Degrees

  if food-type = "Freshly Hatched"
   [set food-factor f_min + random-float (f_max - f_min)]

  if food-type = "Two Week Old"
  ;[set food-factor f_min + random-float (f_max - f_min)]
  [ set food-factor 1.13 ] ;Kristi

  ;if food-type = "Two Week"
  ;[ set food-factor 1.15 ] ;Kristi

  if food-type = "Two Week Old" and day = foodfactor-change
  [set food-factor f_min + random-float (f_max - f_min)]
  ;[ set food-factor 1.62 ] ;Kristi

  if food-type = "Juvenile" and lifestage = "FH" ;freshly hatched cohort
  [set food-factor 2.56] ;Kristi

  if food-type = "Juvenile" and lifestage = "2WFF"
  [set food-factor 1.13 ] ;Kristi

  if food-type = "Juvenile" and lifestage = "2WL"
  [set food-factor 1.13 ] ;Kristi

  if food-type = "Juvenile" and day >= foodfactor-change and lifestage = "2WFF"
  [ set food-factor 1.99 ] ;Kristi


  if food-type = "Juvenile" and day >= foodfactor-change and lifestage = "2WL"
  [set food-factor 2.38 ] ;Kristi

  if food-type ="AdultTest"
  [set food-factor f_min + random-float (f_max - f_min)]

  if food-type = "AdultTest" and day >= foodfactor-change
  [set food-factor f_min + random-float (f_max - f_min)]

set p_Xm ( p_Xm * food-factor )
set F_m  ( F_m  * food-factor )


end



; =========================================================================================================================================
; =============================================== GO PROCEDURE: RUNNING THE MODEL =========================================================
; =========================================================================================================================================
; The go statement below is the order in which all procedures are run each timestep


to go

  ifelse (ticks / timestep <= stopsim or stopsim = 0) ; as long as the provided end time is not exceeded, simulations continue
  [

    T_corrections
    if day = 0 [food-corr ]
    if day = foodfactor-change [food-corr ]
    toxicokinetics
    copper_toxicity

    ask patches [
      calculate_potential_feeding               ; calculate feeding potential per patch internally if necessary
      calc-dU_E ]

    ask turtles with [t_juvenile = 1]           ; all individuals calculate the change in their state variables based on the current conditions
    [
      calc-dU_H
      calc-dU_R
      calc-dL
      calc-dq_acceleration
      calc-dh_rate
    ]

    density-dependent-mortality

    ask turtles
      [ update-individuals                      ; the the state variables of the individuals updated based on the delta value
        metabolic-acceleration ]                ; metabolic-acceleratin by Zimmer (2014)

    update-environment                          ; the the state variables of the environment are updated based on the delta value
    if Die? [ death? ]

    ask turtles with [t_adult = 1]              ; mature individuals check if they have enough energy to reproduce
    [
      set t_repro-time t_repro-time + ((T_fact * 1) / timestep) ; taking into account time between reproduction
      if t_repro-time > time-between-repro
      [ lay-eggs ]
    ]

    if Movement and count patches > 1 [ move ]  ; if enabled, all individuals move

    tick

    if (plotting) [ do-plots ]                  ; then the plots are updated
    set day ticks / timestep
    write_Data_output

    if count turtles = 0 [stop]

  ]

  [ stop ]

end



; - - - - - - - - - - - - - - - - - - - - - - - - - - - - The MC button  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; The MC statement runs the model for multiple iterations, determin by N_iterations


to MC
  let i 0

  if RecordData
  [
    set-current-directory pathName
    if file-exists? OutputFilename [ user-message ("Output file already present, please delete or select another output file" )]
    write-Heading
  ]

  while [i < N_iterations]
  [
    set i i + 1
    print i

    setup
    copper_toxicity
    set iter_N i

    let j 0
    while [ j < stopsim * timestep ]
    [
      set j j + 1
      go
    ]
  ]
end



; - - - - - - - - - - - - - - - - - - - - - - - Effect parameters calibration algorithm - - - - - - - - - - - - - - - - - - - - - - - - -
; Algorithm to asses calibrate the effect parameters A and H
; Algorithm can also be used to let the model run for multiple copper concentrations


to pop_asses

  let i 0

  ;let loc_metal (list 0); controls only for food calibrations
  ;let loc_metal (list 0 50.6); Mattson 18 degrees
  ;let loc_metal (list 0 24.7); Mattson free ion 18 degrees
  ;let loc_metal (list 0 54.6); Mattson 26 degrees
  ;let loc_metal (list 0 26.6); Mattson free ion 26 degrees
  ;let loc_metal (list 0 60) ;Mattson nominal Ni
  ;let loc_metal (list 0 50.6) ;Mattson dissolved Ni measurements 18 degrees
  ;let loc_metal (list 0 1.93 3.38 5.97 12.67 30.56 81.99 198.57) ;FH dissolved concentrations
  ;let loc_metal (list 0 1.07 3.06 5.72 13.01 31.89 77.39 201.16) ;2W FF dissolved concentrations
  ;let loc_metal (list 0 1.06 2.88 5.31 12.24 29.87 73.72 194.13) ;2W L dissolved concentrations
  let loc_metal (list 0 1.93 3.38 5.97 12.67 30.56 81.99 198.57 0 1.07 3.06 5.72 13.01 31.89 77.39 201.16 0 1.06 2.88 5.31 12.24 29.87 73.72 194.13) ; dissolved Ni concentrations
  ;let loc_metal (list 0 0.261 0.456 0.806 1.710 4.126 11.076 26.857 0 0.206 0.589 1.101 2.504 6.140 14.905 38.780 0 0.208 0.564 1.041 2.399 5.855 14.455 38.100) ;free Ni ion concentrations




  ;let loc_exposure (list 0) ;controls only for food calibrations
  ;let loc_exposure (list 0 0 0 0 0 0 0 0) ;For freshly hatched calibrations
  let loc_exposure (list 0 0 0 0 0 0 0 0 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14) ;For combined freshly hatched and two week calibrations

  ;let loc_lifestage(list "2W")
  ;let loc_lifestage(list "AdultTest")
  let loc_lifestage(list "FH" "FH" "FH" "FH" "FH" "FH" "FH" "FH" "2WFF" "2WFF" "2WFF" "2WFF" "2WFF" "2WFF" "2WFF" "2WFF" "2WL" "2WL" "2WL" "2WL" "2WL" "2WL" "2WL" "2WL")

  write-Heading

  while [i < N_iterations]
  [
    set i (i + 1)
    print i

    ; set effect parameters
 ;   if food-type = "Fish flakes"
 ;   [
 ;     set A_min 77.2
 ;     set A_max 77.6
 ;     set H_min 4.68
 ;     set H_max 4.72
 ;   ]

  ;  if food-type = "Lettuce"
  ;  [
  ;    set A_min 63.8
  ;    set A_max 64.4
  ;    set H_min 1.65
  ;    set H_max 1.69
  ;  ]

    let A_temp A_min + random-float (A_max - A_min)              ; randomize A between given min and max
    let H_temp H_min + random-float (H_max - H_min)              ; randomize H between given min and max

    let counter1 0
    let counter2 0
    let counter3 0

    while [counter2 < length loc_metal]
    [
      set counter2 counter2 + 1
      set counter1 counter1 + 1
      set counter3 counter3 + 1
      setup
      set A_eq A_temp
      set H_eq H_temp
      set metal-conc (item (counter2 - 1) loc_metal)
      set t_exposure (item (counter1 - 1) loc_exposure)
      set lifestage (item (counter3 - 1) loc_lifestage)

      ;setup
      ;toxicokinetics
      copper_toxicity
      set iter_N i

      let j 0
      while [ j < stopsim * timestep ]
       [
         set j j + 1
         go
       ]

    ]
  ]
end



; =========================================================================================================================================
; ================================================= INITIAL ENERGY ========================================================================
; =========================================================================================================================================
; Calculate the initial energy of the first individuals using a bisection method


to calc-embryo-reserve-investment
  set embryo-timestep timestep * 50
  set lower-bound 0
  set upper-bound 1
  set sim 0
  loop
  [
    set sim sim + 1
    set egg_weight .5 * (lower-bound + upper-bound)
    set L_embryo  L_0
    set U_E_embryo egg_weight
    set U_H_embryo  0
    set e_scaled_embryo v_rate_int * (U_E_embryo / L_embryo  ^ 3)
    set e_ref 1
    while [U_H_embryo < U_H^b_int and e_scaled_embryo > 1 ]
      [
        set e_scaled_embryo v_rate_int * (U_E_embryo / L_embryo  ^ 3)
        set S_C_embryo   L_embryo  ^ 2 * (g_int * e_scaled_embryo / (g_int + e_scaled_embryo)) * (1 + (L_embryo  / (g_int * (v_rate_int / ( g_int * k_M_rate_int)))))
        set dU_E_embryo  ( -1 * S_C_embryo )
        set dU_H_embryo  ((1 - kap_int) * S_C_embryo - k_J_rate_int * U_H_embryo  )
        set dL_embryo    ((1 / 3) * (((V_rate_int /( g_int * L_embryo  ^ 2 )) * S_C_embryo) - k_M_rate_int * L_embryo ))
        set U_E_embryo   U_E_embryo +  dU_E_embryo    / (embryo-timestep )
        set U_H_embryo   U_H_embryo  +  dU_H_embryo   / (embryo-timestep )
        set L_embryo     L_embryo  +  dL_embryo    / (embryo-timestep )
        set e_scaled_embryo v_rate_int * (U_E_embryo / L_embryo  ^ 3)
      ]
    if e_scaled_embryo <  (.01 +  e_ref) and e_scaled_embryo > (-.01 + e_ref) and U_H_embryo  >= U_H^b_int
      [
        stop
      ]
    ifelse U_H_embryo  > U_H^b_int
      [
        set upper-bound egg_weight
      ]
      [
        set lower-bound egg_weight
      ]
    if sim > 100
    [
      user-message ("Embryo submodel did not converge. Timestep may need to be smaller." )  stop  ;if the timestep is too big relative to the speed of growth of species this will no converge
    ]
  ]
end



; =========================================================================================================================================
; ========================================================== SUBMODELS ====================================================================
; =========================================================================================================================================


; - - - - - - - - - - - - - - - - - - - - - - - - Individual variability - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to individual-variability
  ; Individuals vary in their DEB paramters on a normal distribution with a mean on the input paramater and a coefficent of variation equal to the cv
  ; Set cv to 0 for no variation
  ;set cv random-float 0.25

  set t_scatter-multiplier e ^ (random-normal 0 cv)
  set t_p_Xm_rate p_Xm * t_scatter-multiplier
  set t_p_Am_rate t_p_Xm_rate * 0.8
  set t_g g_int / t_scatter-multiplier
  set t_U_H^b U_H^b_int / t_scatter-multiplier
  set t_U_H^p U_H^p_int / t_scatter-multiplier
  set t_U_H^j U_H^j_int / t_scatter-multiplier

  set t_v_rate v_rate_int
  set t_kap kap_int
  set t_kap_R kap_R_int
  set t_k_M_rate k_M_rate_int
  set t_k_J_rate k_J_rate_int
  ;set t_K t_p_Xm_rate / F_m

  set t_M 1
  set t_i 1
  set t_j 1
  set t_X 1

  ; Energetics

  set t_adult 0
  set t_juvenile 1
  set t_die? 0
  set t_offspring-number 0
  set t_L L_embryo * t_scatter-multiplier
  set t_U_E U_E_embryo * t_scatter-multiplier
  set t_U_R 0
  set t_dU_R  0
  set t_U_H t_U_H^b
  set t_dU_H 0
  set t_h_rate 0
  set t_dh_rate 0
  set t_q_acceleration 0
  set t_dq_acceleration 0

end



; - - - - - - - - - - - - - - - - - - - - - - - - Temperature effects - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; The effect of the temperature is incorporated through the Arrhenius equations
; The different rates are corrected with a temperature correction factor T_fact


to T_corrections
  if Temperature-dependency = "none" [ set T_fact 1 ]
  if Temperature-dependency = "constant" [ set T_fact exp (T_arrh / 293 - T_arrh / temperature) ]
  if Temperature-dependency = "from file" [ set T_fact exp (T_arrh / 293 - T_arrh / temperature) ]
end


;------------------Toxicokinetics-----------------------------------------------------------------------------------------

;to define_TK
 ;set k_e 0.0109 ; calibrated elimination rate Mattson 18 degrees
 ;set k_u 0.000272 ; calibrated uptake rate Mattson 18 degrees
 ;set k_e 0.02511 ; calibrated elimination rate Mattson 26 degrees
 ;set k_u 0.000467 ; calibrated uptake rate Mattson 26 degrees
 ;set k_e 0.0131 ; temp-corrected rate for 20 degrees
 ;set k_u 0.000328 ; temp-corrected rate for 20 degrees


;  if calibrate_uptake? = "Yes"
;  [ set k_u k_u_min + random-float (k_u_max - k_u_min) ;calibrating elimination rate Kristi
    ;set k_e 0.0109 ]; calibrated elimination rate Mattson 18 degrees
;    set k_e 0.02511 ]; calibrated elimination rate Mattson 26 degrees


;end

to toxicokinetics ; V3 Kristi
    ask turtles
   [if day < t_exposure
   [set t_cV 0
   set t_dcV 0]



  ; Calculates the change in scaled internal concentration over tim
   ;ask turtles
   ;set t_cV 0
   ;if day < t_exposure [set t_cV 0 ]

    if day >= t_exposure and metal-conc > 0
    [ ask turtles with [t_L > 0]  ;V5 Kristi Weighman
    [set t_dcV ((t_Lm / t_L) * k_u * metal-conc) - (k_e * t_cV) - (t_cV * ((3 / t_L) * (t_dL / timestep)))

    if day > t_endexposure and metal-conc > 0
    ;[ set t_dcV ( 0 - ( k_e * t_cV)) ]
    [ set t_dcV ( 0 - ( k_e * t_cV) - ( t_cV * (( 3 / t_L ) * (t_dL / timestep))) )]

    set t_cV t_cV + t_dcV

    ]
    ]]



end
; - -
; - - - - - - - - - - - - - - - - - - - - - - - - Copper toxicity effects - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to copper_toxicity
  if metal-present?
  [
    ask turtles
    [
      set t_f_stress 0
      set t_g_stress 0
      set t_k_M_stress1 0
      set t_k_M_stress2 0
      set t_k_J_stress 0
      set t_repro_stress1 0

      set t_stress logistic-dose-response H_eq A_eq t_cV
      if t_stress < 0 [ set t_stress 0 ]
      if t_stress = 1 [ stop ]

      if PMoA = "Maintenance"
      [
        set t_k_J_stress t_stress
        set t_k_M_stress1 t_stress
      ]

      if PMoA = "Assimilation"
      [
        set t_f_stress t_stress
      ]

      if PMoA = "Growth"
      [
        set t_g_stress t_stress
        set t_k_M_stress2 t_stress
      ]

      set t_k_J_rate ( k_J_rate_int / (1 - t_k_J_stress) )
      set t_g t_g / (1 - t_g_stress)
      set t_k_M_rate (k_M_rate_int / (1 - t_k_M_stress1) * (1 - t_k_M_stress2))
      set t_kap_R (kap_R_int * (1 - t_repro_stress1))
    ]
  ]
end



; - - - - - - - - - - - - - - - - - - - - - - - - Reserve dynamics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; Change in reserves: determined by the difference between assimilation (S_A) and mobilization (S_C) fluxes
; When food-dynamics are constant f = the value of f_scaled set in the user interface
; If food is set to  "logistic" f depends on prey density and the half-saturation coefficient (K)
; For embryos f = 0 because they do not feed exogenously


to calculate_potential_feeding ;
    ask turtles-here with [ t_juvenile = 1 ]
    [
      set t_functresp (p_food / (p_food + t_K))
      set t_S_A  t_functresp * t_L ^ 2
    ]

    set p_feedingpotential ( sum [ t_S_A * (T_fact * t_p_Xm_rate) ] of turtles-here ) / timestep
end



to calc-dU_E ; Change in energy reserves
  ; Calculate the functional response, taking into account the amount of available food (avoids that individuals in a patch eat more than availabe)

  if (Food-dynamics = "logistic") or (Food-dynamics = "constant")
  [
    ifelse (p_feedingpotential >= p_food)    ; If the amount of food is too small to support all the individuals in the patch
    [
      set p_actualfeeding p_food             ; Individuals eat exactly the amount of food available
      ask turtles-here with [t_juvenile = 1]
      [
        ifelse (p_food = 0)
        [ set t_functresp 0]
        [ set t_functresp (p_food / (p_food + t_K)) * (p_actualfeeding / p_feedingpotential) ]        ; The functional response function is calculated with a correction for the amount of food available
      ]
    ]
    [                                                 ; If the amount of food is large enough to support all the individuals in the patch
      set p_actualfeeding p_feedingpotential          ; The amount of food that will be eaten is set to the calculated feeding by all individuals
    ]
  ]


  if (Food-dynamics = "scaled")
  [
    ask turtles-here
    [
      set t_functresp f_scaled
    ]
  ]


  ask turtles-here with [t_juvenile = 1] [
    set t_functresp t_functresp * (1 - t_f_stress)


    set t_S_A  t_functresp * t_L ^ 2

    set t_e_scaled (T_fact * t_v_rate) * (t_U_E / t_L ^ 3)
    set t_S_C t_L ^ 2 * (t_g * t_e_scaled / (t_g + t_e_scaled)) * (1 + (t_L / (t_g * ((T_fact * t_v_rate) / (t_g * (T_fact * t_K_M_rate) )))))

    set t_dU_E (t_S_A - t_S_C)
  ]

end



; - - - - - - - - - - - - - - - - - - - - - - - - Maturity and Reproduction - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; Change in maturity is calculated (for immature individuals only)


to calc-dU_H
  ifelse t_U_H < (t_U_H^p / T_Fact) ; they only invest into maturity until they reach puberty
    [set t_dU_H ((1 - t_kap) * t_S_C - (T_Fact * t_k_J_rate) * t_U_H) ]
    [set t_dU_H 0]
end



; Change in reprobuffer (mature individuals only)
to calc-dU_R
  if t_U_H >= (t_U_H^p / T_Fact)
    [set t_dU_R  ((1 - t_kap) * t_S_C - (T_Fact * t_k_J_rate) * (t_U_H^p / T_Fact)) ]
end



; - - - - - - - - - - - - - - - - - - - - - - - - Dynamics of structural length - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; The following procedure calculates change in structural length, if growth is negative the individual does not have enough energy to pay somatic maintenance and the starvation submodel is run
; Where growth is set to 0 and individuals divirt enough energy from development (for juveniles) or reprodution (for adults) to pay maintenance costs


to calc-dL

  set t_dL ((1 / 3) * (((T_Fact * t_v_rate /( t_g * t_L ^ 2 )) * t_S_C) - (T_Fact * t_K_M_rate) * t_L))

  if t_dL <= 0  ; If growth is negative use starvation strategy 3 from the DEB book            t_e_scaled < t_L / ((T_Fact * t_v_rate) / ( t_g * (1 + g_stress) * (T_Fact * t_K_M_rate)))              t_dL <=0
  ; Structural maintenance is dominant over maturity maintenance
    [
      set t_dL 0 ; No structural growth

      ifelse t_U_H < (t_U_H^p / T_Fact)
       [ set t_dU_H t_dU_H - t_kap * t_L ^ 2 * ( t_L / ((T_fact * t_v_rate) / (t_g * (T_fact * t_K_M_rate))) - t_e_scaled)]
       ; Juveniles will divert energy from maturity to pay maintenance costs
       [ set t_dU_R t_dU_R - t_kap * t_L ^ 2 * ( t_L / ((T_fact * t_v_rate) / (t_g * (T_fact * t_K_M_rate))) - t_e_scaled)]
       ; Adults will divert energy for reproduction to pay maintenance costs

      set t_dU_E  t_S_A - t_e_scaled * t_L ^ 2 ; The mobilized energy flux S_C is simplified > the change is energy reserves is recalculated

      ifelse t_U_H < (t_U_H^p / T_Fact)
      [ if t_dU_H < 0 [ set t_die? 1 ] ]
      [ if t_dU_R < 0 [ set t_die? 1 ] ]
      ; The individual cannot pay maintenance costs, even though it diverts energy from maturity/reproduction
      ; Maintenance costs always have to be paid, otherwise the individual dies
    ]

end



; - - - - - - - - - - - - - - - - - - - - - - - - Lay eggs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;The following procedure is run for mature individuals which have enough energy to reproduce
;They create new offsprings and give it the following state variables and DEB parameters


to lay-eggs

  set t_offspring-number  floor ((t_U_R * t_kap_R) / (egg_cost / T_fact))

  if t_offspring-number > 0
  [
    set t_repro-time 0
    set t_U_R t_U_R - ((floor ((t_U_R ) / (egg_cost / T_fact))) * (egg_cost / T_fact))
    set cum_repro cum_repro + ( sum [ t_offspring-number ] of turtles-here ) ; Calculate cumulative reproduction, number of produced and hatched eggs
  ]

  if count turtles > 1 ;Kristi
  [
    hatch t_offspring-number
    [
      set t_adult 0
      set t_juvenile 0
      set t_die? 0
      set t_offspring-number 0
      set t_L 0
      set t_U_E 0
      set t_U_H 0
      set t_U_R 0
      set t_dU_R  0
      set t_h_rate 0
      set t_dh_rate 0
      set t_q_acceleration 0
      set t_dq_acceleration 0
      set t_S_A 0
      set t_S_C 0
      set t_dL 0

      set t_develop-time sim * timestep * T_Fact
    ]

  ]
end



; - - - - - - - - - - - - - - - - - - - - - - - - Ageing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; The following procedure calculates the change in damage enducing compounds of an individual


to calc-dq_acceleration
  set t_dq_acceleration (t_q_acceleration * (t_L ^ 3 / ((T_Fact * t_v_rate) / (t_g * (T_Fact * t_k_M_rate ))) ^ 3) * sG + H_a) * t_e_scaled * ((T_Fact * t_v_rate / t_L) - ((3 / t_L)*  t_dL)) - ((3 / t_L ) * t_dL) * t_q_acceleration
end



; The following procedure calculates the change in damage in the individual
to calc-dh_rate
  set t_dh_rate t_q_acceleration - ((3 / t_L) * t_dL) * t_h_rate
end



; - - - - - - - - - - - - - - - - - - - - - - - - Movement - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; Calculates the chance of movement according to the meanDispersal rate
; Unit of meanDispersal must correlate with the surface area > 1 m/h > area 1 = 1 m2


to move                                                             ; This will calculate the distance and direction of movement
  ask turtles with [t_juvenile = 1] [ if random-float 1 < ( meanDispersal / timestep / sqrt ([p_area] of patch-here) ) [move-to one-of neighbors with [p_slice_n > 0 ]] ]   ; You can take the sqrt of p_area because patches are assumed squares
end



; - - - - - - - - - - - - - - - - - - - - - - - - Metabolic acceleration - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; Metabolic acceleration subprocess for Lymnaea stagnalis, as poposed by Zimmer (2014)


to metabolic-acceleration

  if t_U_H < (t_U_H^b / T_Fact)
  [ set t_M 1
    ;set t_food-factor food-factor
    ]
  ;set t_food-factor 4.025]

  if t_U_H < (t_U_H^j / T_Fact) and t_U_H > (t_U_H^b / T_Fact)
  [
    set t_M t_L / ( L_Embryo * t_scatter-multiplier )
    ;set t_food-factor food-factor
  ;set t_food-factor 4.025]
  ]

  if t_U_H > (t_U_H^j / T_Fact)
  [
    set t_M t_L^j / ( L_Embryo * t_scatter-multiplier)
    ;set t_food-factor 1
  ]

  set t_p_Am_rate t_M * 0.8 * p_Xm * t_scatter-multiplier
  set t_v_rate t_M * v_rate_int

  set t_Lm  kap_int * t_p_Am_rate / p_M

  set t_U_H^b E_H^b / t_p_Am_rate / t_scatter-multiplier
  set t_U_H^p E_H^p / t_p_Am_rate / t_scatter-multiplier
  set t_k_M_rate p_M / E_G * t_scatter-multiplier
  set t_g (E_G * t_v_rate / t_p_Am_rate) / t_kap / t_scatter-multiplier
  set t_U_H^j E_H^j / t_p_Am_rate / t_scatter-multiplier

  set t_p_Xm_rate p_Xm * t_M * t_scatter-multiplier

end



; - - - - - - - - - - - - - - - - - - - - - - - - Update invdividuals - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to update-individuals

  set t_U_E t_U_E + t_dU_E / timestep
  set t_U_H t_U_H + t_dU_H / timestep
  set t_U_R t_U_R + t_dU_R / timestep
  if (t_U_R < 0) [ set t_U_R 0 ]
  set t_L t_L + t_dL / timestep

  if (t_U_H >= (t_U_H^b / T_Fact))
  [
    set t_q_acceleration t_q_acceleration + t_dq_acceleration  / timestep
    set t_h_rate t_h_rate + t_dh_rate  / timestep
  ]

    if (t_juvenile = 0)
  [
    set t_develop-time t_develop-time - (1 * T_fact)
    if (t_develop-time <= 0) [ individual-variability ]
  ]


  if t_i = 1
  [
    if t_U_H >= (t_U_H^j / T_Fact)
    [
      set t_L^j t_L                                  ; Determine size at transformation
      set t_i 0
    ]
  ]

  if t_j = 1
  [
    if t_U_H >= (t_U_H^p / T_Fact)
    [
      set t_L^p t_L                                  ; Determine size at puberty
      set t_j 0
    ]
  ]

  if t_X = 1
  [
    if t_L > 21 / 10 * shape_factor
    [
      set t_XL day
      set t_X 0
    ]
  ]
set t_bmass (t_L ^ 3) * 0.15 + (((t_U_E + t_U_R) * p_Am * 2) / 45795.12);Conversions: 12 g C/mol C, 2 g DW/g C, 550 000 J/mol C, 39.42= p_am, Specific density for dry weight (g/cm^3): 0.15 Kristi

end



to  death?

  ask turtles
  [
    if (t_h_rate > 1) [ set t_h_rate 1 ]             ; This was added to avoid non-numbers being produced (see next line)
    if random-float 1 <  ( 1 - (1 - t_h_rate) ^ ( 1 / timestep))
    [ set t_die? 1 ]
  ]

  if (metal-present? and mortality?)                 ; Effect of the metal on survival
  [
    ask turtles with [ t_juvenile = 1 and t_adult = 0 ]
    [
      if random-float 1 < mort_juv                   ; Generates a random, if it is lower than the chance of dying, the individual dies
      [ set t_die? 1 ]
    ]

    ask turtles with [ t_adult = 1 ]
    [
      if random-float 1 < mort_adult                 ; Generates a random, if it is lower than the chance of dying, the individual dies
      [ set t_die? 1 ]
    ]
  ]

  ask turtles with [t_die? = 1]
  [ die ]

end



; - - - - - - - - - - - - - - - - - - - - - - - - Update environment - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to update-environment
  ask turtles
  [
    if t_U_H >= (t_U_H^b / T_Fact) [set t_juvenile 1]
    if t_U_H >= (t_U_H^p / T_Fact) [set t_adult 1]
  ]

  set day (ticks / timestep)

  if (Food-dynamics = "logistic")                      ; The following procedure calculates change in prey density when prey dynamics are logistic
  [
    ask patches
    [
      set p_d_food ((r_food * T_Fact) * p_food * (1 - (p_food / K_food))) / timestep - p_actualfeeding
      set p_food p_food + p_d_food
      if p_food < (1) [ set p_food (1) ]
    ]
  ]


  if cum_repro = 0
  [
    set age_1strepro day + 1                            ; Keep age of first reproduction
  ]

  if ticks mod timestep = 0 and Temperature-dependency = "from file"
  [
    let index ticks / timestep
    set temperature (item index temperature-matrix) + 273
  ]
end



; - - - - - - - - - - - - - - - - - - - - - - - - Denstiy dependent mortality - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to density-dependent-mortality

  if dd-mortality = "e-juvenile-only"
  [
    if any? turtles with [t_juvenile = 1 and t_adult = 0]
    [
      ask turtles with [t_juvenile = 1 and t_adult = 0 ]
      [
        if random-float 1 < ( 1 - t_e_scaled ) * ( 1 - (1 - mortality-constant) ^ ( 1 / timestep))
        [
          set t_die? 1
        ]
      ]
    ]
  ]

  if dd-mortality = "e-adult-only"
  [
    if any? turtles with [t_adult = 1]
    [
      ask turtles with [t_adult = 1 ]
      [
        if random-float 1 < ( 1 - t_e_scaled ) * ( 1 - (1 - mortality-constant) ^ ( 1 / timestep))
        [
          set t_die? 1
        ]
      ]
    ]
  ]

  if dd-mortality = "e-all"
  [
    if any? turtles with [t_juvenile = 1]
    [
      ask turtles with [t_juvenile = 1 ]
      [
        if random-float 1 < ( 1 - t_e_scaled ) * ( 1 - (1 - mortality-constant) ^ ( 1 / timestep))
        [
          set t_die? 1
        ]
      ]
    ]
  ]

end



; =========================================================================================================================================
; =============================================================== PlOT ====================================================================
; =========================================================================================================================================


to set-patch-vals [number phy cpw]
 ; print (word "bugs " number " " phy " " cpw)
  ask patches with [p_slice_n = number]
  [
    set p_food phy
  ]
end



to do-plots

  set-current-plot "population density"
  set-plot-pen-interval 1 / timestep
  plot count turtles with [t_juvenile = 1]

 ; set-current-plot "length"
 ; set-plot-pen-interval 1 / timestep
 ; plot mean [t_L ] of turtles with [t_juvenile = 1]

  ifelse Food-dynamics = "scaled"
    [set-current-plot "food density"
    plot 0]
    [
    set-current-plot "food density"
    set-plot-pen-interval 1 / timestep
    plot mean [p_food] of patches]

  set-current-plot "age structure"
  set-current-plot-pen "embryo"
    set-plot-pen-interval 1 / timestep
    ifelse any? turtles with [t_juvenile = 0] [plot count turtles with [t_juvenile = 0]]
    [plot 0]

  set-current-plot-pen "juveniles"
    set-plot-pen-interval 1 / timestep
    ifelse any? turtles with [t_juvenile = 1 and t_adult = 0] [plot count turtles with [t_juvenile = 1 and t_adult = 0]]
    [plot 0]

  set-current-plot-pen "adults"
    set-plot-pen-interval 1 / timestep
    ifelse any? turtles with [t_adult = 1] [plot count turtles with [t_adult = 1]]
    [plot 0]

  set-current-plot "cummulative reproduction"
    set-plot-pen-interval 1 / timestep
    plot cum_repro

  set-current-plot "size distribution adults"
  set-current-plot-pen "default"
    histogram [t_L / shape_factor ] of turtles with [t_adult = 1]   ;  shape_factor  = 0.4272 for Lymnaea (Zimmer, 2014) (*10 to get [mm] in stead of [cm])

   set-current-plot "size distribution juveniles"
   set-current-plot-pen "default"
    histogram [t_L / shape_factor ] of turtles with [t_adult = 0 and t_juvenile = 1]

  set-current-plot "size distribution population"
  set-current-plot-pen "default"
    histogram [t_L / shape_factor ] of turtles with [t_juvenile = 1]

end



; =========================================================================================================================================
; ==================================================== OTHER FUNCTIONS ====================================================================
; =========================================================================================================================================


; - - - - - - - - - - - - - - - - - - - - - - - - Read txt-files - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to-report read-in-matrix [filename ]
  let loc-mat []
  carefully [
    file-open filename
    set loc-mat file-read
    file-close
  ]
  [
    print (word "There was a problem with importing file " filename)
    user-message (word "There was a problem with importing file " filename)
    report []
  ]
  report loc-mat
end



; - - - - - - - - - - - - - - - - - - - - - - - - Write data to txt-file - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to write_Data_output


  if (RecordData and ticks mod (timestep) = 0) [
    set-current-directory pathName
    file-open OutputFilename

    file-write iter_N
    ifelse metal-present? [file-write metal-conc ] [ file-write 0]
    file-write day
    file-write t_exposure

    file-write lifestage
    file-write Food-type                                           ; which food type is selected
    file-write Food-factor                                         ; the factor for food
    file-write cv
    file-write mean [p_food] of patches                            ; food density [J]


    file-write count turtles with [t_juvenile = 1 ]                ; we do not include the eggs

    file-write count turtles with [t_juvenile = 0]                 ; print out embryos
    file-write count turtles with [t_juvenile = 1 and t_adult = 0] ; print out juveniles
    file-write count turtles with [t_adult = 1]                    ; print out adults
   ; file-write count turtles with [t_U_H > (t_U_H^j / T_Fact)]


    let size_adult 0
    ifelse count turtles with [t_adult = 1] > 0
    [ set size_adult mean [t_L / shape_factor * 10] of turtles with [t_adult = 1 ]
      file-write size_adult ]
    [ file-write 0 ]                                               ; mean size adults [mm]

    let size_juv 0
    ifelse count turtles with [t_juvenile = 1 and t_adult = 0 ] > 0
    [ set size_juv mean [t_L / shape_factor * 10] of turtles with [t_juvenile = 1 and t_adult = 0 ]
      file-write size_juv ]
    [ file-write 0 ]                                               ; mean size juveniles [mm]

    let size_pop 0
    ifelse count turtles with [t_juvenile = 1] > 0
    [ set size_pop mean [t_L / shape_factor * 10] of turtles with [t_juvenile = 1]
      file-write size_pop ]
    [ file-write 0]

    file-write cum_repro                                           ; cumulative reproduction [#eggs]
    file-write age_1strepro                                        ; age of first reproduction


    file-write A_eq
    file-write H_eq
    file-write k_u
    file-write k_e

    let int_conc 0
    ifelse count turtles with [t_juvenile = 1] > 0
    [ set int_conc mean [t_cV ] of turtles with [t_juvenile = 1]
      file-write int_conc ]
    [ file-write 0]

    ;let change_int_conc 0
    ;ifelse count turtles with [t_juvenile = 1] > 0
    ;[ set change_int_conc mean [t_dcV ] of turtles with [t_juvenile = 1]
    ;  file-write change_int_conc ]
    ;[ file-write 0]


    let popwwt 0
    ifelse count turtles > 0
    [ set popwwt sum[t_bmass] of turtles / 0.15 + sum [t_bmass] of turtles * 0.10
    file-write popwwt]
    [file-write 0]

    let popwwtsoft 0
    ifelse count turtles > 0
    [ set popwwtsoft sum[t_bmass] of turtles / 0.15
    file-write popwwtsoft]
    [file-write 0]

    let juvwwt 0
    ifelse count turtles with [ t_juvenile = 1 and t_adult = 0 ]> 0
    [ set juvwwt mean[t_bmass] of turtles with [ t_juvenile = 1 and t_adult = 0 ] / 0.15 + mean[t_bmass] of turtles with [ t_juvenile = 1 and t_adult = 0 ] * 0.10
    file-write juvwwt]
    [file-write 0]

    let juvwwtsoft 0
    ifelse count turtles with [ t_juvenile = 1 and t_adult = 0 ] > 0
    [ set juvwwtsoft mean[t_bmass] of turtles with [ t_juvenile = 1 and t_adult = 0 ] / 0.15
    file-write juvwwtsoft]
    [file-write 0]

    let adultwwt 0
    ifelse count turtles with [t_adult = 1]  > 0
    [ set adultwwt mean[t_bmass] of turtles with [t_adult = 1] / 0.15 + mean[t_bmass] of turtles with [t_adult = 1 ] * 0.10
    file-write adultwwt]
    [file-write 0]

    let adultwwtsoft 0
    ifelse count turtles with [t_adult = 1] > 0
    [ set adultwwtsoft mean[t_bmass] of turtles with [t_adult = 1 ] / 0.15
    file-write adultwwtsoft]
    [file-write 0]

    ;let volume 0
    ;ifelse count turtles with [t_juvenile = 1] > 0
    ;[ set volume mean [t_L] of turtles with [t_juvenile = 1] ^ 3
    ;  file-write volume ]
    ;[ file-write 0]

    ;let stress_level 0
    ;ifelse count turtles with [t_juvenile = 1] > 0
    ;[ set stress_level mean [t_stress ] of turtles with [t_juvenile = 1]
    ;file-write stress_level ]
    ;[ file-write 0]

    ;file-write Lm

    ;let percentultimate 0
    ;ifelse count turtles with [t_juvenile = 1] > 0
    ;[ set percentultimate mean [t_L] of turtles with [t_juvenile = 1] / Lm
    ;file-write percentultimate ]
    ;[ file-write 0]

    let struct_size_pop 0
    ifelse count turtles with [t_juvenile = 1] > 0
    [ set struct_size_pop mean [t_L] of turtles with [t_juvenile = 1] * 10
    file-write struct_size_pop ]
    [ file-write 0]

    ;let maturity_pop 0
    ;ifelse count turtles with [t_juvenile = 1] > 0
    ;[ set maturity_pop mean [t_U_H] of turtles with [t_juvenile = 1]
    ;file-write maturity_pop ]
    ;[ file-write 0]

    ;let reserves_pop 0
    ;ifelse count turtles with [t_juvenile = 1] > 0
    ;[ set reserves_pop mean [t_U_E] of turtles with [t_juvenile = 1]
    ;file-write reserves_pop ]
    ;[ file-write 0]

    ;let reprobuffer_pop 0
    ;;ifelse count turtles with [t_juvenile = 1] > 0
    ;[ set reprobuffer_pop mean [t_U_R] of turtles with [t_juvenile = 1]
    ;file-write reprobuffer_pop ]
    ;[ file-write 0]

    let reserveperstruc_pop 0
    ifelse count turtles with [t_juvenile = 1] > 0
    [ set reserveperstruc_pop mean [t_e_scaled] of turtles with [t_juvenile = 1]
    file-write reserveperstruc_pop ]
    [ file-write 0]

    file-print ""
    file-close


  ]
end



; - - - - - - - - - - - - - - - - - - - - - - - - Write heading (names) of data to txt-file - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to write-Heading
  if RecordData [
    set-current-directory pathName
    file-open OutputFilename
    file-write "Iteration"
    file-write "Copper"
    file-write "Day"
    file-write "ExposureDay"
    file-write "LifestageTest"
    file-write "FoodType"
    file-write "FoodFactor"
    file-write "cv"
    file-write "FoodDensity"
    file-write "CountIndividuals"
    file-write "CountEmbryo"
    file-write "CountJuvenile"
    file-write "CountAdult"
    ;file-write "CountAcceleration"
    file-write "MeanLengthAdults"
    file-write "MeanLengthJuveniles"
    file-write "PopulationLength"
    file-write "Cumulativeoffspring"
    file-write "Age-1stReproduction"
    file-write "A_eq"
    file-write "B_eq"
    file-write "UptakeRate"
    file-write "EliminationRate"
    file-write "InternalConcentration"
    ;file-write "ChangeInternalConcentration"
    file-write "PopulationWetWeight"
    file-write "PopulationSoftTissueWetWeight"
    file-write "MeanJuvenileWetWeight"
    file-write "MeanJuvenileSoftTissueWetWeight"
    file-write "MeanAdultWetWeight"
    file-write "MeanAdultSoftTissueWetWeight"
    ;file-write "Volume"
    ;file-write "Stress"
    ;file-write "Lm"
    ;file-write "PortionUltimateLength"
    file-write "PopulationStructuralLength"
    ;file-write "PopulationMaturity"
    ;file-write "PopulationReserves"
    ;file-write "PopulationReproBuffer"
    file-write "PopulationReservesPerStructure"
    file-print ""
    file-close
  ]
end


; - - - - - - - - - - - - - - - - - - - - - - - - Logistic dose response function - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


to-report logistic-dose-response [H EC50 x]                                                ; determines the survival chance using a log-logistic curve with given slope, LC50 and exposure concentration
  ifelse (x > 0) [
    report 1 / (1 + (exp (- H * (ln x - ln EC50 ))))                                       ; log-logistic formula
  ]
  [
    report 0
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
1379
580
1624
776
0
0
165.0
1
10
1
1
1
0
1
1
1
0
0
0
0
1
1
1
ticks
30.0

BUTTON
32
10
117
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
32
43
117
76
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
32
76
117
109
go-once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
134
54
197
114
timestep
24
1
0
Number

INPUTBOX
196
54
259
114
stopsim
42
1
0
Number

MONITOR
134
10
197
55
turtles
count turtles
17
1
11

MONITOR
196
10
259
55
days
day
17
1
11

INPUTBOX
258
54
390
114
N_iterations
10000
1
0
Number

BUTTON
32
109
117
142
NIL
MC
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
134
114
259
174
initial_number_bugs
1
1
0
Number

SWITCH
444
10
572
43
RecordData
RecordData
0
1
-1000

INPUTBOX
444
43
759
103
OutputFileName
Juveniles_CalibrateEP_DissolvedNi_Growth_A0to1_B0to5
1
0
String

INPUTBOX
444
103
759
163
pathName
.
1
0
String

TEXTBOX
962
43
1112
61
Intraspecific variation
11
0.0
1

INPUTBOX
982
67
1059
127
cv
0
1
0
Number

TEXTBOX
1387
22
1589
50
Temperature related parameters
11
0.0
1

CHOOSER
1386
46
1570
91
Temperature-dependency
Temperature-dependency
"none" "from file" "constant"
0

INPUTBOX
1570
31
1699
91
Temperature-constant
293
1
0
Number

INPUTBOX
1386
91
1699
151
temperature-file
./Temperature_Ref.txt
1
0
String

TEXTBOX
1388
203
1538
221
Food dynamics
11
0.0
1

CHOOSER
1386
226
1524
271
Food-dynamics
Food-dynamics
"constant" "logistic" "scaled" "intermittent"
2

TEXTBOX
1699
219
1849
237
Metal toxicity parameters
11
0.0
1

TEXTBOX
1391
335
1577
365
if constant or if logistic
11
0.0
1

TEXTBOX
1392
427
1542
445
if scaled
11
0.0
1

INPUTBOX
1388
353
1465
413
K_food
200
1
0
Number

INPUTBOX
1465
353
1542
413
r_food
1.4
1
0
Number

SLIDER
1390
444
1507
477
f_scaled
f_scaled
0
1
0.99
0.01
1
NIL
HORIZONTAL

SWITCH
1700
239
1832
272
metal-present?
metal-present?
0
1
-1000

INPUTBOX
1700
272
1768
332
metal-conc
3.06
1
0
Number

SWITCH
30
226
135
259
plotting
plotting
1
1
-1000

PLOT
135
226
807
587
population density
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

SWITCH
809
43
929
76
Movement
Movement
1
1
-1000

INPUTBOX
30
259
135
319
shape_factor
0.4272
1
0
Number

PLOT
807
425
1344
587
food density
NIL
NIL
0.0
10.0
0.0
1.1
true
false
"" ""
PENS
"food_density" 1.0 0 -16777216 true "" ""

PLOT
807
226
1344
425
age structure
NIL
NIL
0.0
300.0
0.0
10.0
true
true
"" ""
PENS
"adults" 1.0 0 -13791810 true "" ""
"juveniles" 1.0 0 -11085214 true "" ""
"embryo" 1.0 0 -2674135 true "" ""

PLOT
460
587
807
836
size distribution adults
NIL
NIL
0.0
4.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.25 1 -16777216 true "" ""

PLOT
135
587
460
836
size distribution juveniles
NIL
NIL
0.0
1.5
0.0
10.0
true
false
"" ""
PENS
"default" 0.1 1 -16777216 true "" ""

PLOT
807
587
1081
836
size distribution population
NIL
NIL
0.0
4.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.25 1 -16777216 true "" ""

PLOT
1081
587
1344
836
cummulative reproduction
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

SWITCH
809
109
929
142
mortality?
mortality?
1
1
-1000

BUTTON
1752
73
1849
106
NIL
pop_asses
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
1740
52
1873
80
Population assessment
11
0.0
1

MONITOR
32
142
117
187
Species
species_clear_name
17
1
11

CHOOSER
1135
45
1248
90
dd-mortality
dd-mortality
"e-juvenile-only" "e-adult-only" "e-all" "off"
0

INPUTBOX
1135
90
1248
150
mortality-constant
0.1
1
0
Number

TEXTBOX
1123
23
1286
41
Density dependent mortality
11
0.0
1

SWITCH
809
76
929
109
Die?
Die?
1
1
-1000

CHOOSER
1386
271
1525
316
Food-type
Food-type
"Lettuce" "Fish flakes" "Periphyton" "Brix" "Freshly Hatched" "Mattson" "calibrate" "Juvenile" "AdultTest" "Two Week Old"
7

MONITOR
1524
271
1600
316
NIL
food-factor
17
1
11

MONITOR
258
10
328
55
NIL
temperature
17
1
11

MONITOR
327
10
390
55
NIL
iter_N
17
1
11

TEXTBOX
818
23
929
41
Process controllers
11
0.0
1

MONITOR
1700
365
1764
410
NIL
A_eq
17
1
11

MONITOR
1764
365
1832
410
NIL
H_eq
17
1
11

CHOOSER
1402
495
1540
540
PMoA
PMoA
"Assimilation" "Growth" "Maintenance"
1

INPUTBOX
1708
421
1923
481
A_min
0
1
0
Number

INPUTBOX
1722
496
1937
556
A_max
1
1
0
Number

INPUTBOX
1715
587
1930
647
H_min
0
1
0
Number

INPUTBOX
1727
652
1942
712
H_max
5
1
0
Number

INPUTBOX
1854
276
2009
336
t_endexposure
43
1
0
Number

INPUTBOX
1733
119
1793
179
f_min
0.5
1
0
Number

INPUTBOX
1735
183
1797
243
f_max
10
1
0
Number

CHOOSER
1858
51
1996
96
calibrate_uptake?
calibrate_uptake?
"Yes" "No"
1

INPUTBOX
1864
116
2019
176
k_u_min
1.0E-5
1
0
Number

INPUTBOX
1868
188
2099
248
k_u_max
0.001
1
0
Number

INPUTBOX
2005
372
2160
432
t_exposure
14
1
0
Number

INPUTBOX
2034
479
2281
539
lifestage
2WFF
1
0
String

INPUTBOX
2056
285
2211
345
foodfactor-change
14
1
0
Number

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.3.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
