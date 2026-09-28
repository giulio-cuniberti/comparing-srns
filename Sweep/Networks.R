# Source and product complexes of the 205 benchmark networks.
# Rows represent reactions; columns represent species in the same order in S and P.
# BioModels identifiers and model names are retained. See README.md for provenance.

networks <- list(
  list(
    source = "HARV", id = "BIOMD0000000138",
    name = "Tabak2007_dopamine",
    d = 1, n = 1,
    S = matrix(c(
      0
    ), nrow = 1, byrow = TRUE),
    P = matrix(c(
      1
    ), nrow = 1, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000119",
    name = "Golomb2006_SomaticBursting_nonzero[Ca]",
    d = 1, n = 1,
    S = matrix(c(
      0
    ), nrow = 1, byrow = TRUE),
    P = matrix(c(
      1
    ), nrow = 1, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001026",
    name = "Kurlovics2021 - Metformin partitioning between plasma and RBC with independent Kin and Kout coefficients",
    d = 1, n = 2,
    S = matrix(c(
      1,
      0
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      0,
      1
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001040",
    name = "Kurlovics2021 - Metformin partitioning from plasma to RBC,  single coefficient",
    d = 1, n = 2,
    S = matrix(c(
      1,
      0
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      0,
      1
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000309",
    name = "Tyson2003_NegFB_Homeostasis",
    d = 1, n = 2,
    S = matrix(c(
      0,
      1
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      1,
      0
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000417",
    name = "Ratushny2012_NF",
    d = 1, n = 2,
    S = matrix(c(
      0,
      1
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      1,
      0
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000418",
    name = "Ratushny2012_SPF",
    d = 1, n = 2,
    S = matrix(c(
      0,
      1
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      1,
      0
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000311",
    name = "Tyson2003_Mutual_Activation",
    d = 1, n = 3,
    S = matrix(c(
      0,
      0,
      1
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      1,
      1,
      0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000414",
    name = "Band2012_DII-Venus_ReducedModel",
    d = 1, n = 3,
    S = matrix(c(
      0,
      1,
      1
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      1,
      0,
      0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000617",
    name = "Walsh2014 - Inhibition kinetics of DAPT on APP Cleavage",
    d = 1, n = 4,
    S = matrix(c(
      0,
      0,
      1,
      1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1,
      1,
      0,
      0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001079",
    name = "DeBoeck2021 - Modular approach to modeling the cell cycle, simple cell cycle model",
    d = 1, n = 4,
    S = matrix(c(
      0,
      1,
      1,
      0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1,
      0,
      0,
      1
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000076",
    name = "Cronwright2002_Glycerol_Synthesis",
    d = 1, n = 4,
    S = matrix(c(
      0,
      1,
      1,
      0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1,
      0,
      0,
      1
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000310",
    name = "Tyson2003_Mutual_Inhibition",
    d = 1, n = 4,
    S = matrix(c(
      0,
      0,
      1,
      1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1,
      1,
      0,
      0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000425",
    name = "Tan2012 - Antibiotic Treatment, Inoculum Effect",
    d = 1, n = 4,
    S = matrix(c(
      0,
      2,
      1,
      2
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1,
      1,
      2,
      1
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000068",
    name = "Curien2003_MetThr_synthesis",
    d = 1, n = 6,
    S = matrix(c(
      0,
      1,
      1,
      0,
      1,
      0
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      1,
      0,
      0,
      1,
      0,
      1
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000124",
    name = "Wu2006_K+Channel",
    d = 2, n = 2,
    S = matrix(c(
      0, 0,
      0, 0
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      1, 0,
      0, 1
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000047",
    name = "Oxhamre2005_Ca_oscillation",
    d = 2, n = 3,
    S = matrix(c(
      1, 0,
      1, 0,
      0, 1
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      0, 1,
      0, 1,
      1, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000307",
    name = "Tyson2003_Substrate_Depletion_Osc",
    d = 2, n = 3,
    S = matrix(c(
      0, 1,
      0, 0,
      1, 0
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      1, 0,
      0, 1,
      0, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000225",
    name = "Westermark2003_Pancreatic_GlycOsc_basic",
    d = 2, n = 3,
    S = matrix(c(
      0, 0,
      0, 1,
      1, 0
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      0, 1,
      1, 0,
      0, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000601",
    name = "Rosas2015 - Caffeine-induced luminal SR calcium changes",
    d = 2, n = 3,
    S = matrix(c(
      1, 0,
      0, 1,
      1, 0
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      0, 1,
      1, 0,
      0, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000385",
    name = "Arnold2011_Schultz2003_RuBisCO-CalvinCycle",
    d = 2, n = 4,
    S = matrix(c(
      2, 0,
      2, 0,
      1, 1,
      0, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 2,
      0, 2,
      1, 0,
      1, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000386",
    name = "Arnold2011_Sharkey2007_RuBisCO-CalvinCycle",
    d = 2, n = 4,
    S = matrix(c(
      2, 0,
      2, 0,
      1, 1,
      0, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 2,
      0, 2,
      1, 0,
      1, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000387",
    name = "Arnold2011_Damour2007_RuBisCO-CalvinCycle",
    d = 2, n = 4,
    S = matrix(c(
      2, 0,
      2, 0,
      1, 1,
      0, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 2,
      0, 2,
      1, 0,
      1, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000787",
    name = "Frascoli2014 - A dynamical model of tumour immunotherapy",
    d = 2, n = 4,
    S = matrix(c(
      0, 0,
      1, 1,
      0, 1,
      1, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 1,
      1, 0,
      1, 1,
      0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000825",
    name = "Greene2019 - Differentiate Spontaneous and Induced Evolution to Drug Resistance During Cancer Treatment",
    d = 2, n = 4,
    S = matrix(c(
      0, 1,
      1, 0,
      0, 1,
      0, 1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1, 0,
      1, 1,
      1, 1,
      0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001078",
    name = "Hammaren-Geissen2022_PPToP_Model12",
    d = 2, n = 4,
    S = matrix(c(
      0, 1,
      0, 1,
      1, 0,
      1, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 0,
      1, 0,
      0, 1,
      0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000024",
    name = "Scheper1999_CircClock",
    d = 2, n = 4,
    S = matrix(c(
      0, 1,
      1, 0,
      1, 0,
      0, 1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1, 1,
      1, 1,
      0, 0,
      0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000155",
    name = "Zatorsky2006_p53_Model6",
    d = 2, n = 4,
    S = matrix(c(
      0, 0,
      1, 1,
      1, 0,
      0, 1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1, 0,
      0, 1,
      1, 1,
      0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000420",
    name = "Ratushny2012_ASSURE_I",
    d = 2, n = 4,
    S = matrix(c(
      0, 0,
      1, 0,
      0, 0,
      0, 1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1, 0,
      0, 0,
      0, 1,
      0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000421",
    name = "Ratushny2012_ASSURE_II",
    d = 2, n = 4,
    S = matrix(c(
      0, 0,
      1, 0,
      0, 0,
      0, 1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1, 0,
      0, 0,
      0, 1,
      0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001057",
    name = "Nikolov2020 - p53-miR34 model",
    d = 2, n = 5,
    S = matrix(c(
      0, 0,
      0, 1,
      1, 0,
      1, 0,
      0, 1
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      1, 0,
      1, 1,
      0, 0,
      1, 1,
      0, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000191",
    name = "Montañez2008_Arginine_catabolism",
    d = 2, n = 5,
    S = matrix(c(
      1, 1,
      0, 1,
      1, 1,
      1, 0,
      0, 1
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      0, 2,
      1, 1,
      1, 0,
      0, 0,
      0, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000973",
    name = "Dasgupta2020 - Reduced model of receptor clusturing and aggregation",
    d = 2, n = 6,
    S = matrix(c(
      1, 1,
      0, 1,
      1, 0,
      0, 1,
      1, 0,
      1, 0
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      1, 0,
      0, 0,
      1, 1,
      1, 1,
      0, 0,
      0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000306",
    name = "Tyson2003_Activator_Inhibitor",
    d = 2, n = 6,
    S = matrix(c(
      0, 0,
      0, 0,
      1, 0,
      1, 1,
      1, 0,
      0, 1
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      1, 0,
      1, 0,
      0, 0,
      0, 1,
      1, 1,
      0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000458",
    name = "Smallbone2013 - Serine biosynthesis",
    d = 2, n = 6,
    S = matrix(c(
      1, 0,
      2, 0,
      2, 1,
      1, 2,
      0, 2,
      0, 1
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      2, 0,
      1, 0,
      1, 2,
      2, 1,
      0, 1,
      0, 2
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000098",
    name = "Goldbeter1990_CalciumSpike_CICR",
    d = 2, n = 6,
    S = matrix(c(
      0, 0,
      0, 0,
      0, 1,
      1, 0,
      1, 0,
      0, 1
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      0, 1,
      0, 1,
      1, 0,
      0, 1,
      0, 1,
      0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000117",
    name = "Dupont1991_CaOscillation",
    d = 2, n = 6,
    S = matrix(c(
      0, 0,
      0, 0,
      1, 0,
      0, 1,
      1, 0,
      0, 1
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      0, 1,
      0, 1,
      0, 1,
      1, 0,
      0, 1,
      0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001060",
    name = "Frank2021 - Macrophage polarization",
    d = 2, n = 7,
    S = matrix(c(
      0, 0,
      0, 0,
      1, 0,
      0, 0,
      0, 0,
      0, 0,
      0, 1
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      1, 0,
      1, 0,
      0, 0,
      0, 1,
      0, 1,
      0, 1,
      0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000316",
    name = "Shen-Orr2002_FeedForward_AND_gate",
    d = 2, n = 8,
    S = matrix(c(
      0, 0,
      1, 0,
      1, 0,
      0, 0,
      1, 0,
      1, 1,
      0, 1,
      0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      1, 0,
      0, 0,
      0, 0,
      1, 0,
      1, 1,
      1, 0,
      0, 0,
      0, 1
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001080",
    name = "DeBoeck2021 - Modular approach to modeling the cell cycle, 5 ODE model with 3 bistable switches",
    d = 2, n = 8,
    S = matrix(c(
      0, 0,
      0, 1,
      0, 1,
      0, 0,
      0, 0,
      1, 0,
      1, 0,
      0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 1,
      0, 0,
      0, 0,
      0, 1,
      1, 0,
      0, 0,
      0, 0,
      1, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000749",
    name = "Reppas2015 - tumor control via alternating immunostimulating and immunosuppressive phases",
    d = 2, n = 9,
    S = matrix(c(
      0, 0,
      0, 1,
      0, 1,
      0, 0,
      0, 1,
      1, 1,
      1, 0,
      0, 0,
      0, 1
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      0, 1,
      0, 0,
      0, 0,
      0, 1,
      1, 1,
      0, 1,
      0, 0,
      1, 0,
      0, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000956",
    name = "Bertozzi2020 - SIR model of scenarios of COVID-19 spread in CA and NY",
    d = 3, n = 2,
    S = matrix(c(
      0, 0, 1,
      1, 0, 0
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 1, 0
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000963",
    name = "Weitz2020 - SIR model of COVID-19 transmission with shielding",
    d = 3, n = 2,
    S = matrix(c(
      0, 1, 1,
      1, 0, 0
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      1, 1, 0,
      0, 1, 0
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001045",
    name = "Smith&Moore2004 - The SIR model for the spread of HongKong Flu",
    d = 3, n = 2,
    S = matrix(c(
      0, 0, 1,
      1, 0, 0
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 1, 0
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000917",
    name = "Phillips2007_AscendingArousalSystem_SleepWakeDynamics",
    d = 3, n = 3,
    S = matrix(c(
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 0, 1,
      0, 1, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000128",
    name = "Bertram2006_Endothelin",
    d = 3, n = 3,
    S = matrix(c(
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 0, 1,
      0, 1, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000025",
    name = "Smolen2002_CircClock",
    d = 3, n = 4,
    S = matrix(c(
      0, 0, 1,
      0, 0, 1,
      1, 0, 0,
      0, 1, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1, 0, 1,
      0, 1, 1,
      0, 0, 0,
      0, 0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000027",
    name = "Markevich2004 - MAPK double phosphorylation,  ordered Michaelis-Menton",
    d = 3, n = 4,
    S = matrix(c(
      1, 0, 0,
      1, 1, 0,
      1, 0, 1,
      0, 1, 1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      1, 0, 1,
      1, 1, 0,
      1, 0, 1
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000031",
    name = "Markevich2004_MAPK_orderedMM2kinases",
    d = 3, n = 4,
    S = matrix(c(
      1, 0, 0,
      1, 1, 0,
      1, 0, 1,
      0, 1, 1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      1, 0, 1,
      1, 1, 0,
      1, 0, 1
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000253",
    name = "Teusink1998_Glycolysis_TurboDesign",
    d = 3, n = 4,
    S = matrix(c(
      1, 0, 0,
      1, 0, 1,
      0, 1, 0,
      1, 0, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 0, 1,
      0, 1, 0,
      4, 0, 0,
      0, 0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000319",
    name = "Decroly1982_Enzymatic_Oscillator",
    d = 3, n = 4,
    S = matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 50, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "BM", id = "BIOMD0000000866",
    name = "Simon2019 - NIK-dependent p100 processing into p52, Michaelis-Menten, SBML 2v4",
    d = 3, n = 4,
    S = matrix(c(
      0, 0, 0,
      1, 1, 0,
      0, 0, 1,
      0, 1, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      1, 0, 1,
      0, 0, 0,
      0, 0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000546",
    name = "Miao2010 - Innate and adaptive immune responses to primary Influenza A Virus infection_1_1",
    d = 3, n = 5,
    S = matrix(c(
      1, 0, 1,
      0, 0, 0,
      0, 1, 0,
      0, 0, 1,
      0, 1, 0
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      0, 1, 1,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 1, 1
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000789",
    name = "Jenner2018 - treatment of oncolytic virus",
    d = 3, n = 5,
    S = matrix(c(
      1, 0, 0,
      0, 0, 1,
      0, 0, 0,
      0, 1, 1,
      1, 0, 0
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      1, 0, 1,
      0, 0, 0,
      0, 1, 0,
      1, 0, 1,
      0, 0, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000814",
    name = "Perez-Garcia19 - Computational design of improved standardized chemotherapy protocols for grade 2 oligodendrogliomas",
    d = 3, n = 5,
    S = matrix(c(
      1, 0, 0,
      0, 1, 1,
      0, 1, 1,
      1, 0, 1,
      0, 1, 0
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      1, 0, 1,
      1, 1, 0,
      0, 1, 0,
      0, 0, 1,
      0, 0, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000860",
    name = "Proctor2017- Role of microRNAs in osteoarthritis (Positive Feedforward Incoherent By MicroRNA)_1",
    d = 3, n = 5,
    S = matrix(c(
      1, 0, 0,
      0, 0, 1,
      1, 0, 0,
      0, 1, 0,
      0, 1, 1
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      1, 0, 1,
      0, 0, 0,
      1, 1, 0,
      0, 0, 0,
      0, 0, 1
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000927",
    name = "Grigolon2018-Responses to auxin signals",
    d = 3, n = 5,
    S = matrix(c(
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      1, 1, 0,
      0, 1, 1
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      0, 0, 1,
      0, 1, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001056",
    name = "Chulian2021 - feedback signalling in B lymphopoeisis",
    d = 3, n = 5,
    S = matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 1, 0,
      0, 0, 1
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 1, 0,
      0, 0, 1,
      0, 0, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000712",
    name = "Manchanda2014 - Effect on Immune System by 4 different Influenza A virus strains",
    d = 3, n = 6,
    S = matrix(c(
      0, 0, 0,
      1, 0, 1,
      0, 0, 1,
      1, 0, 0,
      0, 0, 0,
      0, 1, 0
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      0, 0, 1,
      1, 0, 0,
      1, 0, 1,
      0, 0, 0,
      0, 1, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000804",
    name = "Koenders2015 - multiple myeloma",
    d = 3, n = 6,
    S = matrix(c(
      1, 0, 1,
      0, 1, 1,
      0, 1, 0,
      0, 1, 0,
      1, 0, 0,
      0, 0, 1
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      1, 1, 1,
      1, 1, 1,
      0, 1, 1,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000900",
    name = "Bianca2013 - Persistence analysis in a Kolmogorov-type model for cancer-immune system competition",
    d = 3, n = 6,
    S = matrix(c(
      0, 0, 0,
      1, 1, 1,
      0, 1, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      1, 0, 1,
      1, 1, 0,
      0, 0, 0,
      0, 1, 1,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000931",
    name = "Voliotis2019-GnRH Pulse Generation",
    d = 3, n = 6,
    S = matrix(c(
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001048",
    name = "Siddhartha2002 - Kinetic modelling of cancer therapies",
    d = 3, n = 6,
    S = matrix(c(
      0, 0, 1,
      0, 1, 0,
      0, 1, 0,
      0, 0, 0,
      0, 0, 0,
      1, 1, 0
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      0, 0, 0,
      1, 1, 0,
      1, 0, 0,
      0, 0, 1,
      1, 0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000308",
    name = "Tyson2003_NegFB_Oscillator",
    d = 3, n = 6,
    S = matrix(c(
      0, 0, 0,
      1, 1, 0,
      0, 1, 0,
      0, 0, 1,
      0, 0, 1,
      1, 0, 0
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      1, 0, 0,
      0, 1, 1,
      0, 0, 0,
      1, 0, 1,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000461",
    name = "Liebal2012 - B.subtilis transcription inhibition model",
    d = 3, n = 6,
    S = matrix(c(
      0, 1, 0,
      0, 2, 0,
      1, 2, 2,
      2, 1, 2,
      0, 2, 1,
      0, 1, 2
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      0, 2, 0,
      0, 1, 0,
      2, 1, 2,
      1, 2, 2,
      0, 1, 2,
      0, 2, 1
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000419",
    name = "Ratushny2012_SPF_I",
    d = 3, n = 6,
    S = matrix(c(
      0, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 1, 0,
      0, 0, 0,
      0, 0, 1
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 0, 0,
      0, 1, 0,
      0, 0, 0,
      0, 0, 1,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000158",
    name = "Zatorsky2006_p53_Model2",
    d = 3, n = 7,
    S = matrix(c(
      0, 0, 0,
      1, 0, 0,
      1, 1, 0,
      1, 0, 0,
      0, 0, 1,
      0, 1, 0,
      0, 0, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 0, 0,
      0, 1, 0,
      1, 0, 1,
      0, 1, 0,
      0, 0, 0,
      0, 1, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000845",
    name = "Gulbudak2019.1 - Heterogeneous viral strategies promote coexistence in virus-microbe systems (Lytic)",
    d = 3, n = 7,
    S = matrix(c(
      0, 0, 0,
      0, 1, 1,
      0, 1, 0,
      1, 0, 0,
      1, 0, 0,
      1, 0, 0,
      0, 0, 1
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      1, 0, 1,
      0, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000846",
    name = "Gulbudak2019.2 - Heterogeneous viral strategies promote coexistence in virus-microbe systems (Chronic)",
    d = 3, n = 7,
    S = matrix(c(
      0, 0, 0,
      0, 1, 1,
      0, 1, 0,
      1, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 0, 1
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      1, 0, 0,
      1, 0, 1,
      0, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001007",
    name = "Zhang2007 - Mechanism of DNA damage response (Model1)",
    d = 3, n = 7,
    S = matrix(c(
      1, 0, 0,
      0, 0, 1,
      0, 0, 1,
      1, 0, 0,
      0, 1, 0,
      1, 0, 0,
      0, 1, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      1, 0, 1,
      0, 0, 0,
      1, 0, 1,
      0, 1, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001043",
    name = "Wodarz2001 - Viruses as antitumor weapons",
    d = 3, n = 7,
    S = matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 1, 0,
      1, 0, 1,
      0, 1, 0,
      1, 0, 0,
      0, 0, 1
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      1, 1, 0,
      0, 0, 0,
      1, 1, 0,
      0, 0, 1,
      1, 0, 0,
      1, 0, 1,
      0, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000157",
    name = "Zatorsky2006_p53_Model4",
    d = 3, n = 7,
    S = matrix(c(
      0, 0, 0,
      1, 0, 0,
      1, 1, 0,
      1, 0, 0,
      0, 0, 1,
      0, 1, 0,
      0, 0, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 0, 0,
      0, 1, 0,
      1, 0, 1,
      0, 1, 0,
      0, 0, 0,
      0, 1, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000258",
    name = "Ortega2006 - bistability from double phosphorylation in signal transduction",
    d = 3, n = 7,
    S = matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 1, 1,
      1, 0, 1,
      1, 1, 0,
      1, 0, 1,
      0, 0, 1
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      1, 0, 0,
      1, 0, 1,
      0, 1, 1,
      1, 0, 1,
      1, 1, 0,
      0, 1, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000520",
    name = "Smallbone2013 - Colon Crypt cycle - Version 0",
    d = 3, n = 7,
    S = matrix(c(
      2, 0, 0,
      2, 0, 0,
      2, 0, 0,
      0, 2, 0,
      0, 2, 0,
      0, 2, 0,
      0, 0, 2
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      2, 1, 0,
      3, 0, 0,
      0, 1, 0,
      0, 2, 1,
      0, 3, 0,
      0, 0, 1
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000236",
    name = "Westermark2003_Pancreatic_GlycOsc_extended",
    d = 3, n = 7,
    S = matrix(c(
      0, 0, 0,
      0, 0, 1,
      0, 0, 1,
      0, 1, 0,
      0, 1, 0,
      2, 0, 0,
      1, 0, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 0, 1,
      0, 0, 0,
      0, 1, 0,
      0, 0, 1,
      2, 0, 0,
      0, 1, 0,
      0, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001019",
    name = "Barros2021 - CARTmath, Mathematical Model of CAR-T Immunotherapy in HDLM-2 cell line",
    d = 3, n = 8,
    S = matrix(c(
      1, 0, 1,
      0, 1, 0,
      0, 1, 0,
      0, 1, 1,
      0, 1, 0,
      1, 0, 0,
      0, 0, 0,
      0, 1, 1
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 1, 1,
      0, 2, 0,
      0, 0, 0,
      0, 0, 1,
      1, 1, 0,
      0, 0, 0,
      0, 0, 1,
      0, 1, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001020",
    name = "Barros2021 - CARTmath, Mathematical Model of CAR-T Immunotherapy in Raji Cell Line",
    d = 3, n = 8,
    S = matrix(c(
      1, 0, 1,
      0, 1, 0,
      0, 1, 0,
      0, 1, 1,
      0, 1, 0,
      1, 0, 0,
      0, 0, 0,
      0, 1, 1
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 1, 1,
      0, 2, 0,
      0, 0, 0,
      0, 0, 1,
      1, 1, 0,
      0, 0, 0,
      0, 0, 1,
      0, 1, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000454",
    name = "Smallbone2013 - Metabolic Control Analysis - Example 1",
    d = 3, n = 8,
    S = matrix(c(
      1, 2, 1,
      2, 1, 2,
      0, 1, 2,
      0, 2, 1,
      2, 0, 0,
      1, 0, 0,
      2, 0, 0,
      1, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      2, 1, 2,
      1, 2, 1,
      0, 2, 1,
      0, 1, 2,
      1, 0, 0,
      2, 0, 0,
      1, 0, 0,
      2, 0, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000784",
    name = "Lopez2014 - A Validated Mathematical Model of Tumor Growth Including Tumor-Host Interaction and Cell-Mediated Immune Response",
    d = 3, n = 9,
    S = matrix(c(
      0, 0, 0,
      1, 1, 0,
      1, 0, 0,
      0, 0, 0,
      1, 1, 0,
      0, 0, 0,
      0, 0, 1,
      1, 0, 0,
      1, 0, 1
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 0, 0,
      0, 1, 0,
      1, 0, 0,
      0, 0, 1,
      0, 0, 0,
      1, 0, 1,
      1, 0, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000337",
    name = "Pfeiffer2001_ATP-ProducingPathways_CooperationCompetition",
    d = 3, n = 10,
    S = matrix(c(
      0, 0, 0,
      0, 0, 1,
      0, 0, 1,
      10, 0, 0,
      0, 0, 1,
      0, 1, 0,
      1, 0, 0,
      0, 0, 0,
      0, 1, 0,
      0, 0, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 0, 1,
      0, 0, 0,
      10, 0, 0,
      0, 0, 1,
      0, 1, 0,
      0, 0, 1,
      0, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 1, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000455",
    name = "Smallbone2013 - Metabolic Control Analysis - Example 2",
    d = 3, n = 10,
    S = matrix(c(
      1, 2, 1,
      2, 1, 2,
      0, 1, 2,
      0, 2, 1,
      2, 0, 0,
      1, 0, 0,
      2, 0, 0,
      1, 0, 0,
      0, 0, 2,
      0, 0, 1
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      2, 1, 2,
      1, 2, 1,
      0, 2, 1,
      0, 1, 2,
      1, 0, 0,
      2, 0, 0,
      1, 0, 0,
      2, 0, 0,
      0, 0, 1,
      0, 0, 2
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000317",
    name = "Shen-Orr2002_Single_Input_Module",
    d = 3, n = 12,
    S = matrix(c(
      0, 0, 0,
      1, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 1, 0,
      0, 1, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 1,
      0, 0, 1,
      0, 0, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 0,
      0, 0, 0,
      0, 1, 0,
      0, 0, 1,
      0, 0, 0,
      0, 0, 0,
      0, 0, 1
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001052",
    name = "Alharbi2020 - Tumor and immune system competition",
    d = 3, n = 14,
    S = matrix(c(
      0, 0, 0,
      0, 1, 0,
      1, 1, 0,
      0, 1, 1,
      0, 0, 0,
      0, 0, 1,
      0, 1, 0,
      1, 0, 1,
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1,
      1, 1, 0,
      1, 0, 1
    ), nrow = 14, byrow = TRUE),
    P = matrix(c(
      0, 1, 0,
      0, 0, 0,
      1, 0, 0,
      0, 0, 1,
      0, 0, 1,
      0, 0, 0,
      0, 1, 1,
      1, 0, 0,
      1, 0, 0,
      0, 0, 0,
      1, 1, 0,
      1, 0, 1,
      0, 1, 0,
      0, 0, 1
    ), nrow = 14, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000166",
    name = "Zhu2007_TF_modulated_by_Calcium",
    d = 3, n = 18,
    S = matrix(c(
      0, 0, 0,
      1, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 1,
      0, 0, 0,
      0, 0, 1,
      0, 0, 1,
      0, 1, 0,
      0, 1, 0,
      0, 0, 1,
      0, 1, 0,
      0, 0, 1,
      0, 0, 1,
      0, 0, 0
    ), nrow = 18, byrow = TRUE),
    P = matrix(c(
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      1, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 1,
      0, 0, 0,
      0, 0, 1,
      0, 0, 0,
      0, 1, 0,
      0, 0, 1,
      0, 0, 1,
      0, 1, 0,
      0, 0, 1,
      0, 1, 0,
      0, 0, 0,
      0, 0, 1
    ), nrow = 18, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000978",
    name = "Mukandavire2020 - SEIR model of early COVID-19 transmission in South Africa",
    d = 4, n = 3,
    S = matrix(c(
      0, 1, 0, 1,
      1, 0, 0, 0,
      0, 1, 0, 0
    ), nrow = 3, byrow = TRUE),
    P = matrix(c(
      1, 1, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000150",
    name = "Morris2002_CellCycle_CDK2Cyclin",
    d = 4, n = 4,
    S = matrix(c(
      0, 0, 1, 1,
      1, 0, 0, 0,
      1, 0, 0, 0,
      0, 1, 0, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0,
      0, 0, 1, 1,
      0, 1, 0, 0,
      1, 0, 0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000923",
    name = "Liò2012_Modelling osteomyelitis_Control Model",
    d = 4, n = 4,
    S = matrix(c(
      0, 0, 1, 1,
      1, 0, 0, 1,
      1, 0, 1, 0,
      0, 0, 0, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 1, 1, 1,
      1, 0, 1, 1,
      1, 0, 1, 1,
      1, 0, 0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000976",
    name = "Ghanbari2020 - forecasting the second wave of COVID-19 in Iran",
    d = 4, n = 5,
    S = matrix(c(
      0, 1, 0, 1,
      1, 0, 0, 1,
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 1, 0, 0
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      1, 1, 0, 0,
      1, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 1, 0,
      0, 0, 0, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001032",
    name = "Al-Tuwairqi2020 - Dynamics of cancer radiovirotherapy - Phase II treatment",
    d = 4, n = 7,
    S = matrix(c(
      0, 0, 0, 0,
      0, 1, 0, 1,
      0, 0, 1, 0,
      0, 1, 0, 1,
      0, 0, 0, 1,
      0, 0, 1, 0,
      1, 0, 0, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1,
      0, 1, 1, 0,
      0, 1, 0, 0,
      0, 0, 0, 1,
      1, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000029",
    name = "Markevich2004_MAPK_phosphoRandomMM",
    d = 4, n = 7,
    S = matrix(c(
      1, 1, 0, 0,
      1, 1, 1, 0,
      1, 0, 1, 0,
      1, 1, 1, 0,
      1, 0, 1, 1,
      0, 1, 1, 1,
      0, 1, 1, 1
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 1, 1, 0,
      1, 1, 0, 1,
      0, 1, 1, 0,
      1, 0, 1, 1,
      1, 1, 1, 0,
      1, 0, 1, 1,
      1, 1, 0, 1
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000865",
    name = "Nikolaev2019 - Immunobiochemical reconstruction of influenza lung infection-melanoma skin cancer interactions",
    d = 4, n = 8,
    S = matrix(c(
      1, 0, 1, 0,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      1, 0, 1, 0,
      0, 0, 1, 0,
      1, 0, 0, 0,
      0, 1, 1, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      1, 1, 1, 0,
      0, 0, 0, 0,
      1, 0, 0, 1,
      0, 0, 0, 0,
      1, 0, 2, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      1, 1, 1, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001047",
    name = "Collier1996 - Delta Notch intercellular signalling and lateral inhibition",
    d = 4, n = 8,
    S = matrix(c(
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000044",
    name = "Borghans1997 - Calcium Oscillation - Model 2",
    d = 4, n = 8,
    S = matrix(c(
      0, 1, 0, 0,
      0, 0, 0, 1,
      1, 0, 1, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      1, 0, 0, 1,
      1, 0, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1,
      0, 0, 1, 0,
      1, 0, 0, 1,
      0, 0, 0, 1,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000045",
    name = "Borghans1997 - Calcium Oscillation - Model 3",
    d = 4, n = 8,
    S = matrix(c(
      1, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 1, 0, 0,
      0, 0, 0, 1,
      0, 1, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 1,
      1, 0, 0, 0,
      0, 0, 0, 1,
      0, 1, 0, 0,
      0, 0, 0, 1
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000113",
    name = "Dupont1992_Ca_dpt_protein_phospho",
    d = 4, n = 8,
    S = matrix(c(
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 1, 0, 0,
      1, 1, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1,
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 1,
      0, 0, 0, 0,
      1, 1, 0, 0,
      0, 1, 0, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000462",
    name = "Proctor2012 - Role of Amyloid-beta dimers in aggregation formation",
    d = 4, n = 8,
    S = matrix(c(
      0, 0, 0, 0,
      0, 0, 2, 2,
      0, 0, 3, 0,
      2, 0, 0, 0,
      3, 0, 0, 0,
      0, 2, 2, 0,
      0, 2, 0, 0,
      0, 0, 0, 2
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0,
      0, 0, 1, 2,
      1, 0, 1, 0,
      1, 0, 2, 0,
      1, 4, 0, 0,
      0, 3, 1, 0,
      0, 1, 1, 0,
      0, 0, 0, 1
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "BM", id = "BIOMD0000000869",
    name = "Simon2019 - NIK-dependent p100 processing into p52 and IkBd degradation, Michaelis-Menten, SBML 2v4",
    d = 4, n = 8,
    S = matrix(c(
      0, 0, 0, 0,
      0, 1, 1, 0,
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 2, 0,
      1, 0, 0, 0,
      1, 0, 0, 0,
      1, 1, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0,
      0, 1, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 2, 0,
      0, 0, 0, 0,
      0, 1, 0, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000809",
    name = "Malinzi2018 - tumour-immune interaction model",
    d = 4, n = 9,
    S = matrix(c(
      0, 0, 0, 0,
      0, 0, 1, 0,
      0, 1, 1, 0,
      0, 0, 0, 0,
      0, 1, 1, 0,
      0, 1, 1, 0,
      1, 0, 0, 0,
      0, 1, 1, 0,
      0, 0, 0, 1
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0,
      0, 1, 1, 0,
      0, 0, 1, 0,
      0, 0, 1, 0,
      0, 1, 0, 0,
      1, 1, 1, 0,
      0, 0, 0, 0,
      0, 1, 1, 1,
      0, 0, 0, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000965",
    name = "LeBeau1999 - IP3-dependent intracellular calcium oscillations due to agonist stimulation from Cholecytokinin",
    d = 4, n = 9,
    S = matrix(c(
      0, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 1, 0,
      1, 0, 0, 0,
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 0
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0, 0,
      0, 0, 1, 1,
      0, 0, 0, 0,
      0, 0, 0, 1
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000967",
    name = "McLean1991 - Behaviour of HIV in the presence of zidovudine",
    d = 4, n = 9,
    S = matrix(c(
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 1, 0, 0,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      1, 0, 0, 0
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0,
      0, 0, 1, 1,
      0, 1, 0, 1,
      1, 0, 0, 1,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      1, 0, 0, 1
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000517",
    name = "Smallbone2013 - Colon Crypt cycle - Version 3",
    d = 4, n = 9,
    S = matrix(c(
      2, 0, 0, 2,
      1, 0, 0, 0,
      1, 0, 0, 0,
      0, 2, 0, 2,
      0, 1, 0, 0,
      0, 1, 0, 0,
      0, 0, 2, 2,
      1, 0, 0, 0,
      0, 0, 0, 2
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 2,
      1, 1, 0, 0,
      2, 0, 0, 0,
      0, 1, 0, 2,
      0, 1, 1, 0,
      0, 2, 0, 0,
      0, 0, 1, 2,
      1, 0, 0, 1,
      0, 0, 0, 1
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000518",
    name = "Smallbone2013 - Colon Crypt cycle - Version 2",
    d = 4, n = 9,
    S = matrix(c(
      2, 0, 0, 2,
      1, 0, 0, 0,
      1, 0, 0, 0,
      0, 2, 0, 2,
      0, 1, 0, 0,
      0, 1, 0, 0,
      0, 0, 2, 2,
      1, 0, 0, 0,
      0, 0, 0, 2
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 2,
      1, 1, 0, 0,
      2, 0, 0, 0,
      0, 1, 0, 2,
      0, 1, 1, 0,
      0, 2, 0, 0,
      0, 0, 1, 2,
      1, 0, 0, 1,
      0, 0, 0, 1
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000120",
    name = "Chan2004_TCell_receptor_activation",
    d = 4, n = 10,
    S = matrix(c(
      0, 0, 0, 0,
      1, 0, 1, 0,
      0, 1, 0, 0,
      0, 1, 0, 0,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      1, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 1, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0,
      0, 1, 1, 0,
      1, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 1, 0,
      1, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000775",
    name = "Iarosz2015 - brain tumor",
    d = 4, n = 10,
    S = matrix(c(
      0, 0, 0, 0,
      1, 1, 0, 0,
      0, 1, 0, 1,
      0, 0, 0, 0,
      1, 1, 0, 0,
      1, 0, 0, 1,
      0, 1, 0, 0,
      0, 0, 1, 1,
      0, 0, 0, 0,
      0, 0, 0, 1
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0, 1,
      0, 1, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 1,
      0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001009",
    name = "Zhang2007 - Mechanism of DNA damage response (Model2)",
    d = 4, n = 10,
    S = matrix(c(
      0, 0, 1, 0,
      1, 0, 0, 0,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 1, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      1, 0, 1, 0,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000805",
    name = "Al-Husari2013 - pH and lactate in tumor",
    d = 4, n = 11,
    S = matrix(c(
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 1, 0, 0,
      0, 1, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 1, 0, 0
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0,
      1, 0, 0, 0,
      1, 0, 0, 0,
      1, 0, 1, 0,
      0, 1, 0, 0,
      0, 1, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 1, 0, 1
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000058",
    name = "Bindschadler2001_coupled_Ca_oscillators",
    d = 4, n = 11,
    S = matrix(c(
      0, 0, 1, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 1, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 1, 0, 0
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      1, 0, 1, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 0,
      0, 1, 0, 1,
      0, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      1, 0, 0, 0
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000290",
    name = "Alexander2010_Tcell_Regulation_Sys2",
    d = 4, n = 12,
    S = matrix(c(
      0, 0, 1, 0,
      0, 0, 1, 0,
      0, 1, 0, 0,
      1, 0, 0, 0,
      1, 1, 0, 0,
      1, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      0, 1, 0, 0,
      0, 0, 1, 0,
      1, 0, 0, 0,
      1, 0, 0, 1
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0,
      1, 0, 1, 0,
      0, 1, 1, 0,
      1, 0, 0, 1,
      1, 1, 0, 1,
      1, 1, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000456",
    name = "Smallbone2013 - Metabolic Control Analysis - Example 3",
    d = 4, n = 12,
    S = matrix(c(
      1, 2, 1, 0,
      2, 1, 2, 0,
      0, 1, 2, 0,
      0, 2, 1, 0,
      2, 0, 0, 0,
      1, 0, 0, 0,
      2, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 1,
      0, 0, 0, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      2, 1, 2, 0,
      1, 2, 1, 0,
      0, 2, 1, 0,
      0, 1, 2, 0,
      1, 0, 0, 0,
      2, 0, 0, 0,
      1, 0, 0, 0,
      2, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 1
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000928",
    name = "Baker2017 - The role of cytokines, MMPs and fibronectin fragments osteoarthritis",
    d = 4, n = 13,
    S = matrix(c(
      0, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 1, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      0, 1, 0, 0,
      0, 1, 0, 0,
      1, 1, 0, 0,
      0, 1, 0, 1,
      1, 0, 0, 0
    ), nrow = 13, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1,
      0, 0, 1, 1,
      0, 0, 0, 0,
      0, 0, 1, 0,
      1, 0, 1, 0,
      0, 0, 0, 0,
      1, 1, 0, 0,
      0, 1, 0, 1,
      0, 0, 0, 0,
      1, 1, 0, 0,
      2, 1, 0, 0,
      1, 1, 0, 1,
      0, 0, 0, 0
    ), nrow = 13, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000327",
    name = "Whitcomb2004_Bicarbonate_Pancreas",
    d = 4, n = 20,
    S = matrix(c(
      0, 0, 0, 0,
      2, 0, 0, 1,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 1, 0,
      0, 1, 1, 0,
      1, 0, 0, 0,
      0, 0, 1, 0,
      1, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 1, 0, 0,
      0, 1, 0, 0,
      0, 0, 0, 0
    ), nrow = 20, byrow = TRUE),
    P = matrix(c(
      2, 0, 0, 1,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 0,
      1, 0, 0, 0,
      0, 1, 1, 0,
      1, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 1,
      1, 0, 0, 0,
      0, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 1, 0, 0
    ), nrow = 20, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000104",
    name = "Klipp2002_MetabolicOptimization_linearPathway(n=2)",
    d = 5, n = 2,
    S = matrix(c(
      1, 0, 1, 0, 0,
      0, 1, 0, 1, 0
    ), nrow = 2, byrow = TRUE),
    P = matrix(c(
      0, 1, 1, 0, 0,
      0, 0, 0, 1, 1
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000194",
    name = "Ibrahim2008_Cdc20_Sequestring_Template_Model",
    d = 5, n = 4,
    S = matrix(c(
      0, 0, 1, 0, 1,
      0, 0, 0, 1, 0,
      1, 0, 0, 1, 0,
      0, 1, 0, 0, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1, 0,
      0, 0, 1, 0, 1,
      0, 1, 1, 0, 0,
      1, 0, 0, 0, 1
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000125",
    name = "Komarova2005_TheoreticalFramework_BasicArchitecture",
    d = 5, n = 7,
    S = matrix(c(
      1, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 1, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 1
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      1, 1, 0, 0, 0,
      0, 1, 0, 1, 0,
      0, 1, 1, 0, 0,
      0, 1, 0, 0, 1,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000043",
    name = "Borghans1997 - Calcium Oscillation - Model 1",
    d = 5, n = 7,
    S = matrix(c(
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 1, 1, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1,
      0, 0, 1, 0, 1,
      0, 1, 0, 0, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 1,
      0, 0, 0, 1, 0,
      0, 0, 1, 0, 1,
      0, 0, 0, 0, 1,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 1,
      0, 0, 1, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000231",
    name = "Valero2006_Adenine_TernaryCycle",
    d = 5, n = 8,
    S = matrix(c(
      0, 0, 1, 0, 0,
      0, 1, 0, 0, 0,
      0, 1, 1, 0, 0,
      2, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 1, 0, 1,
      0, 0, 0, 1, 1,
      0, 0, 0, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      2, 0, 0, 0, 0,
      0, 1, 1, 0, 0,
      0, 0, 1, 0, 1,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 1, 1
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000719",
    name = "Tsai2014 - Cell cycle duration control by oscillatory Dynamics  in Early Xenopus laevis Embryos",
    d = 5, n = 9,
    S = matrix(c(
      0, 0, 0, 0, 0,
      1, 0, 0, 1, 0,
      1, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      0, 1, 0, 0, 1,
      1, 0, 0, 0, 0
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 1, 0, 0,
      0, 0, 1, 0, 1,
      0, 0, 0, 0, 0,
      1, 1, 0, 0, 1,
      0, 0, 0, 0, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000121",
    name = "Clancy2001_Kchannel",
    d = 5, n = 10,
    S = matrix(c(
      0, 0, 1, 0, 0,
      0, 1, 0, 0, 0,
      0, 1, 0, 0, 0,
      1, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 1,
      0, 0, 0, 1, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 1, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 0, 0, 1,
      1, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 1, 0,
      1, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000203",
    name = "Chickarmane2006 - Stem cell switch reversible",
    d = 5, n = 10,
    S = matrix(c(
      1, 0, 1, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 1,
      0, 0, 1, 0, 0,
      1, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      1, 0, 1, 0, 0,
      0, 0, 0, 1, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      1, 1, 1, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 1, 0, 1,
      0, 0, 0, 0, 0,
      1, 0, 1, 1, 0,
      0, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000204",
    name = "Chickarmane2006 - Stem cell switch irreversible",
    d = 5, n = 10,
    S = matrix(c(
      1, 0, 1, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 1,
      0, 0, 1, 0, 0,
      1, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      1, 0, 1, 0, 0,
      0, 0, 0, 1, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      1, 1, 1, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 1, 0, 1,
      0, 0, 0, 0, 0,
      1, 0, 1, 1, 0,
      0, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000016",
    name = "Goldbeter1995_CircClock",
    d = 5, n = 10,
    S = matrix(c(
      0, 0, 0, 0, 1,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1,
      1, 0, 0, 0, 0,
      0, 0, 0, 1, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 1,
      1, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000295",
    name = "Akman2008_Circadian_Clock_Model1",
    d = 5, n = 10,
    S = matrix(c(
      0, 0, 1, 1, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 1,
      1, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      0, 1, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 0, 1, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 1, 1,
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 1,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 1,
      0, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 1, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000791",
    name = "Wilson2012 - tumor vaccine efficacy",
    d = 5, n = 11,
    S = matrix(c(
      0, 0, 0, 0, 0,
      1, 1, 0, 1, 0,
      0, 0, 0, 1, 1,
      0, 0, 0, 1, 0,
      1, 0, 0, 0, 0,
      1, 0, 0, 1, 0,
      0, 1, 0, 0, 0,
      0, 1, 1, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 1
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1, 0,
      1, 1, 0, 0, 0,
      0, 0, 0, 0, 1,
      1, 0, 0, 1, 0,
      0, 0, 0, 0, 0,
      1, 1, 0, 1, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000197",
    name = "Bartholome2007_MDCKII",
    d = 5, n = 11,
    S = matrix(c(
      1, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 1
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      0, 1, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      1, 0, 0, 0, 0
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000148",
    name = "Komarova2003_BoneRemodeling",
    d = 5, n = 12,
    S = matrix(c(
      0, 0, 1, 0, 0,
      1, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      1, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      1, 1, 0, 0, 1,
      1, 1, 0, 0, 0,
      0, 0, 1, 1, 0,
      0, 0, 1, 1, 1
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      1, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      1, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      1, 1, 0, 0, 0,
      1, 1, 0, 0, 1,
      0, 0, 1, 1, 1,
      0, 0, 1, 1, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000828",
    name = "Jung2019 - Regulating glioblastoma signaling pathways and anti-invasion therapy - core control model",
    d = 5, n = 12,
    S = matrix(c(
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 1, 0, 0, 0,
      1, 0, 1, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 2, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 1, 0, 1,
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      1, 0, 0, 0, 1,
      0, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      1, 0, 0, 1, 0,
      0, 0, 0, 1, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000325",
    name = "Palini2011_Minimal_2_Feedback_Model",
    d = 5, n = 12,
    S = matrix(c(
      1, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 1, 0,
      0, 1, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 1, 1, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 1,
      1, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 1, 0, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 1, 0,
      0, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      0, 1, 1, 0, 0,
      1, 1, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 1, 0, 0,
      0, 0, 0, 0, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000817",
    name = "Gevertz2018 - cancer treatment with oncolytic viruses and dendritic cell injections minimal model",
    d = 5, n = 13,
    S = matrix(c(
      0, 0, 0, 0, 0,
      0, 1, 0, 1, 1,
      0, 0, 1, 0, 1,
      0, 1, 0, 0, 0,
      0, 1, 0, 1, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 1, 0, 0, 0,
      0, 1, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 1,
      0, 1, 0, 1, 0,
      0, 1, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 1, 1, 0, 0,
      0, 1, 0, 1, 0,
      1, 0, 0, 1, 0,
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000059",
    name = "Fridlyand2003_Calcium_flux",
    d = 5, n = 16,
    S = matrix(c(
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 1, 1, 0, 0,
      0, 1, 0, 0, 0
    ), nrow = 16, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 1, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 1,
      1, 1, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0
    ), nrow = 16, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000023",
    name = "Rohwer2001_Sucrose",
    d = 5, n = 22,
    S = matrix(c(
      0, 0, 0, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      1, 1, 0, 0, 0,
      1, 0, 1, 0, 0,
      1, 1, 0, 0, 0,
      0, 1, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 2, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 1,
      0, 0, 0, 1, 0,
      1, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 1, 0,
      1, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 0
    ), nrow = 22, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 0, 0, 0,
      1, 0, 1, 0, 0,
      1, 1, 0, 0, 0,
      0, 1, 1, 0, 0,
      1, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 2, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 1, 0,
      1, 0, 1, 0, 0,
      1, 1, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 1, 0
    ), nrow = 22, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000305",
    name = "Kolomeisky2003_MyosinV_Processivity",
    d = 6, n = 4,
    S = matrix(c(
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0
    ), nrow = 4, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 1,
      0, 1, 1, 0, 0, 0,
      1, 0, 0, 1, 0, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000624",
    name = "Sluka2016 - Acetaminophen metabolism",
    d = 6, n = 5,
    S = matrix(c(
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0
    ), nrow = 5, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000207",
    name = "Romond1999_CellCycle",
    d = 6, n = 6,
    S = matrix(c(
      0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      0, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000830",
    name = "GiantsosAdams2013 - Growth of glycocalyx under static conditions",
    d = 6, n = 7,
    S = matrix(c(
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000986",
    name = "Aubry1995 - Multi-compartment model of fluid-phase endocytosis kinetics in Dictyostelium discoideum",
    d = 6, n = 8,
    S = matrix(c(
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000864",
    name = "Proctor2017- Role of microRNAs in osteoarthritis (Negative Feedback By MicroRNA)",
    d = 6, n = 9,
    S = matrix(c(
      0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 1, 1, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001015",
    name = "Jarrah2014 - mathematical model of the immune response in muscle degeneration and subsequent regeneration in Duchenne muscular dystrophy in mdx mice",
    d = 6, n = 9,
    S = matrix(c(
      0, 1, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 2, 0,
      0, 1, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      0, 1, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 1, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000005",
    name = "Tyson1991 - Cell Cycle 6 var",
    d = 6, n = 9,
    S = matrix(c(
      0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001053",
    name = "Garde2020 - metabolic oscillations in Bacillus subtilis biofilms",
    d = 6, n = 10,
    S = matrix(c(
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 1,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 1, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000116",
    name = "McClean2007_CrossTalk",
    d = 6, n = 10,
    S = matrix(c(
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 1,
      0, 0, 1, 0, 0, 1
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001042",
    name = "Makhlouf2020 - No treatment model of the role of CD4 T cells in tumor-immune interactions",
    d = 6, n = 12,
    S = matrix(c(
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 1,
      0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 1, 1,
      0, 0, 1, 1, 1, 1,
      0, 1, 0, 0, 1, 1,
      0, 0, 0, 1, 0, 1,
      1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 1, 0,
      0, 0, 1, 0, 1, 1,
      0, 0, 0, 0, 0, 1,
      0, 1, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1,
      1, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001072",
    name = "Phillips2013 - physiologically based modeling explaining Mammalian rest/activity patterns",
    d = 6, n = 12,
    S = matrix(c(
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      1, 1, 0, 1, 0, 0,
      0, 0, 1, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 1, 1, 0, 0,
      0, 1, 1, 1, 0, 0,
      1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0,
      0, 0, 0, 0, 2, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 1, 0,
      0, 0, 0, 1, 1, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000209",
    name = "Chickarmane2008 - Stem cell lineage determination",
    d = 6, n = 12,
    S = matrix(c(
      1, 0, 1, 1, 0, 1,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 1, 1,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 1, 0,
      0, 1, 0, 0, 0, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      1, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 1,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 1, 1,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000213",
    name = "Nijhout2004_Folate_Cycle",
    d = 6, n = 12,
    S = matrix(c(
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000053",
    name = "Ferreira2003_CML_generation2",
    d = 6, n = 12,
    S = matrix(c(
      0, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 1, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1,
      0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000242",
    name = "Bai2003_G1phaseRegulation",
    d = 6, n = 12,
    S = matrix(c(
      0, 0, 0, 0, 0, 0,
      2, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 2, 0, 0, 1, 0,
      0, 0, 1, 1, 1, 0,
      0, 0, 0, 1, 0, 1,
      1, 0, 2, 0, 0, 0,
      0, 1, 2, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 2,
      0, 1, 0, 0, 1, 1,
      0, 0, 0, 0, 2, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 1,
      0, 1, 0, 0, 1, 0,
      0, 0, 1, 2, 1, 0,
      0, 0, 1, 0, 0, 0,
      1, 0, 1, 0, 0, 1,
      0, 1, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 2,
      0, 0, 0, 0, 1, 1,
      0, 1, 0, 0, 2, 1,
      0, 0, 0, 0, 1, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000101",
    name = "Vilar2006_TGFbeta",
    d = 6, n = 13,
    S = matrix(c(
      1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0
    ), nrow = 13, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000351",
    name = "Vernoux2011_AuxinSignaling_AuxinSingleStepInput",
    d = 6, n = 14,
    S = matrix(c(
      1, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 1,
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 0, 1, 0, 0
    ), nrow = 14, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 2, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0,
      1, 0, 0, 0, 0, 0,
      1, 1, 0, 1, 1, 0
    ), nrow = 14, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000352",
    name = "Vernoux2011_AuxinSignaling_AuxinFluctuating",
    d = 6, n = 14,
    S = matrix(c(
      1, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 1,
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 0, 1, 0, 0
    ), nrow = 14, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 2, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0,
      1, 0, 0, 0, 0, 0,
      1, 1, 0, 1, 1, 0
    ), nrow = 14, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000240",
    name = "Veening2008_DegU_Regulation",
    d = 6, n = 14,
    S = matrix(c(
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 2, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0
    ), nrow = 14, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 2, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 1
    ), nrow = 14, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000208",
    name = "Deineko2003_CellCycle",
    d = 6, n = 15,
    S = matrix(c(
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1
    ), nrow = 15, byrow = TRUE),
    P = matrix(c(
      1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0
    ), nrow = 15, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000816",
    name = "Gevertz2018 - Cancer Treatment with Oncolytic Viruses and Dendritic Cell injections original model",
    d = 6, n = 16,
    S = matrix(c(
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 1,
      0, 0, 0, 1, 0, 1,
      0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0
    ), nrow = 16, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 1, 0,
      0, 1, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 0,
      0, 0, 1, 0, 1, 0,
      0, 1, 1, 0, 0, 0,
      0, 1, 0, 0, 1, 0,
      1, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 16, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000245",
    name = "Lei2001_Yeast_Aerobic_Metabolism",
    d = 6, n = 16,
    S = matrix(c(
      0, 1, 0, 1, 0, 1,
      0, 0, 0, 1, 1, 1,
      0, 0, 0, 0, 1, 1,
      1, 1, 0, 0, 0, 1,
      0, 0, 1, 1, 0, 1,
      0, 1, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 1,
      0, 0, 1, 1, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1
    ), nrow = 16, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 1, 1,
      0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 1,
      1, 0, 1, 0, 0, 1,
      0, 0, 0, 1, 0, 1,
      1, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 2,
      0, 0, 0, 1, 0, 2,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 16, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000806",
    name = "Eftimie2019-Macrophages Plasticity",
    d = 6, n = 25,
    S = matrix(c(
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1,
      1, 0, 0, 0, 1, 0,
      0, 0, 1, 1, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      1, 2, 0, 0, 0, 0,
      0, 2, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 2,
      0, 0, 1, 1, 0, 2,
      1, 0, 0, 0, 0, 2,
      0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 1, 0,
      0, 0, 0, 1, 0, 1,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 1, 1, 0, 0,
      1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0
    ), nrow = 25, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 0, 0,
      0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 1,
      1, 1, 0, 0, 0, 0,
      0, 1, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1,
      0, 0, 1, 1, 0, 1,
      1, 0, 0, 0, 0, 1,
      0, 1, 1, 0, 0, 1,
      0, 0, 1, 0, 1, 0,
      0, 0, 1, 1, 0, 0,
      0, 0, 0, 1, 1, 0,
      0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0,
      0, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      1, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0
    ), nrow = 25, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000966",
    name = "Cui2008 - in vitro transcriptional response of zinc homeostasis system in Escherichia coli",
    d = 7, n = 6,
    S = matrix(c(
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 0
    ), nrow = 6, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0,
      1, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 1, 0
    ), nrow = 6, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000820",
    name = "West2019 - Cellular interactions constrain tumor growth",
    d = 7, n = 7,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000077",
    name = "Blum2000_LHsecretion_1",
    d = 7, n = 8,
    S = matrix(c(
      0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 2, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 2, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 1, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000983",
    name = "Zongo2020 - model of COVID-19 transmission dynamics under containment measures in France",
    d = 7, n = 10,
    S = matrix(c(
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 1, 1, 0, 0, 0, 1,
      0, 1, 1, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0,
      1, 1, 1, 0, 0, 0, 0,
      0, 1, 1, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000168",
    name = "Obeyesekere1999_CellCycle",
    d = 7, n = 10,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0,
      2, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 2, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 1,
      0, 1, 0, 0, 1, 2, 0,
      1, 0, 0, 0, 2, 0, 0,
      0, 0, 1, 0, 2, 0, 0,
      0, 1, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 2
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0,
      0, 1, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 1, 1,
      0, 1, 0, 0, 2, 1, 0,
      1, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 1, 0, 0,
      0, 1, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 1
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000248",
    name = "Lai2007_O2_Transport_Metabolism",
    d = 7, n = 10,
    S = matrix(c(
      0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      6, 0, 0, 1, 0, 0, 1,
      0, 6, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 1, 0, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 6, 0, 0, 0, 0, 1,
      6, 0, 0, 1, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 1, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001054",
    name = "Pearce2021 - Fibrin Polymerization",
    d = 7, n = 11,
    S = matrix(c(
      0, 0, 0, 0, 1, 0, 1,
      0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 1,
      1, 0, 0, 0, 0, 0, 0
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 1,
      0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 1, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 1
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000190",
    name = "Rodriguez-Caso2006_Polyamine_Metabolism",
    d = 7, n = 11,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 1, 0, 0,
      0, 1, 0, 1, 0, 0, 0,
      0, 1, 0, 1, 0, 0, 0,
      0, 1, 0, 1, 0, 1, 1,
      0, 0, 0, 1, 0, 1, 1,
      1, 0, 1, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0, 0, 0, 0,
      1, 0, 1, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 1, 0,
      0, 1, 1, 1, 0, 0, 1,
      0, 1, 0, 1, 0, 1, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000210",
    name = "Chickarmane2008 - Stem cell lineage - NANOG GATA-6 switch",
    d = 7, n = 12,
    S = matrix(c(
      1, 0, 1, 1, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 0,
      0, 1, 0, 0, 0, 0, 0
    ), nrow = 12, byrow = TRUE),
    P = matrix(c(
      1, 0, 1, 1, 1, 0, 1,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0, 1,
      0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0
    ), nrow = 12, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000908",
    name = "dePillis2013 - Mathematical modeling of regulatory T cell effects on renal cell carcinoma treatment",
    d = 7, n = 13,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0, 1,
      1, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 1,
      1, 1, 0, 1, 0, 0, 1,
      0, 1, 1, 0, 1, 0, 1,
      1, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0
    ), nrow = 13, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 1, 0, 0,
      1, 1, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1,
      1, 1, 1, 1, 0, 0, 1,
      0, 1, 0, 0, 1, 0, 1,
      1, 1, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000196",
    name = "Srividhya2006_CellCycle",
    d = 7, n = 13,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0, 0,
      0, 1, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 1, 0, 0,
      0, 1, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE),
    P = matrix(c(
      0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 1, 0, 0, 0,
      0, 1, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000232",
    name = "Nazaret2009_TCA_RC_ATP",
    d = 7, n = 13,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 1,
      0, 1, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 1, 0, 0,
      0, 0, 0, 1, 2, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000318",
    name = "Yao2008_Rb_E2F_Switch",
    d = 7, n = 17,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 1, 0, 0,
      1, 1, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0
    ), nrow = 17, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 0,
      0, 1, 1, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      1, 1, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0,
      1, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0
    ), nrow = 17, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000872",
    name = "Verma2016 - HIV and HPV co-infection, T-cell response",
    d = 7, n = 18,
    S = matrix(c(
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 0, 1,
      0, 1, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 2, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 1, 0, 0,
      0, 1, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0
    ), nrow = 18, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 3, 0, 0, 0, 0, 0,
      2, 0, 1, 0, 1, 0, 0,
      0, 1, 0, 2, 0, 0, 0,
      0, 1, 1, 0, 0, 0, 0
    ), nrow = 18, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001006",
    name = "Ciliberto2005 - Steady states and oscillations in the p53/Mdm2 network",
    d = 7, n = 21,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0
    ), nrow = 21, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 1, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 0
    ), nrow = 21, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000918",
    name = "Schwarz2018-Cdk Activity Threshold Determines Passage through the Restriction Point",
    d = 8, n = 7,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 1,
      1, 1, 0, 1, 0, 1, 1, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 1, 0, 0, 0,
      1, 1, 0, 0, 0, 1, 1, 0,
      1, 1, 0, 0, 0, 0, 0, 0
    ), nrow = 7, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1, 0, 0, 0, 1,
      1, 1, 1, 1, 0, 1, 1, 0,
      1, 0, 0, 1, 0, 0, 0, 1,
      0, 1, 1, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 1, 1, 0, 0,
      1, 1, 0, 0, 1, 1, 1, 0,
      1, 1, 0, 0, 0, 0, 1, 0
    ), nrow = 7, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000483",
    name = "Cao2008 - Network of a toggle switch",
    d = 8, n = 8,
    S = matrix(c(
      0, 0, 2, 0, 1, 0, 0, 0,
      0, 0, 0, 2, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 2, 0,
      0, 0, 0, 0, 0, 0, 0, 2,
      0, 0, 0, 2, 0, 0, 3, 0,
      0, 0, 2, 0, 0, 0, 0, 3,
      2, 0, 0, 0, 0, 0, 0, 0,
      0, 2, 0, 0, 0, 0, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 2, 0, 0, 0, 1, 0,
      0, 0, 0, 2, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 1,
      0, 1, 0, 1, 0, 0, 1, 0,
      1, 0, 1, 0, 0, 0, 0, 1,
      1, 0, 1, 0, 0, 0, 0, 2,
      0, 1, 0, 1, 0, 0, 2, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000193",
    name = "Ibrahim2008_MCC_assembly_model_KDM",
    d = 8, n = 9,
    S = matrix(c(
      0, 0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 1
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 1, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000591",
    name = "Boehm2014 - isoform-specific dimerization of pSTAT5A and pSTAT5B",
    d = 8, n = 9,
    S = matrix(c(
      3, 0, 0, 0, 0, 0, 0, 0,
      2, 2, 0, 0, 0, 0, 0, 0,
      0, 3, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 2, 0, 0,
      0, 0, 0, 0, 0, 0, 2, 0,
      0, 0, 0, 0, 0, 0, 0, 2,
      0, 0, 2, 0, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0, 0, 0,
      0, 0, 0, 0, 2, 0, 0, 0
    ), nrow = 9, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 1, 0, 0,
      1, 1, 0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0, 0, 1,
      2, 0, 1, 0, 0, 0, 0, 0,
      1, 1, 0, 1, 0, 0, 0, 0,
      0, 2, 0, 0, 1, 0, 0, 0
    ), nrow = 9, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000001103",
    name = "Palaniappan2021 - Cell free modelling of second generation Toehold switches",
    d = 8, n = 10,
    S = matrix(c(
      0, 0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000862",
    name = "Proctor2017- Role of microRNAs in osteoarthritis (Positive Feedback By Micro RNA)",
    d = 8, n = 11,
    S = matrix(c(
      0, 1, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000960",
    name = "Paiva2020 - SEIAHRD model of transmission dynamics of COVID-19",
    d = 8, n = 11,
    S = matrix(c(
      1, 0, 0, 0, 1, 1, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 1, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 1, 0, 1, 0, 0, 0, 0
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000955",
    name = "Giordano2020 - SIDARTHE model of COVID-19 spread in Italy",
    d = 8, n = 13,
    S = matrix(c(
      1, 1, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1
    ), nrow = 13, byrow = TRUE),
    P = matrix(c(
      1, 1, 0, 0, 1, 1, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000971",
    name = "Tang2020 - Estimation of transmission risk of COVID-19 and impact of public health interventions",
    d = 8, n = 13,
    S = matrix(c(
      1, 0, 0, 0, 1, 0, 1, 0,
      1, 0, 0, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE),
    P = matrix(c(
      1, 1, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 13, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000632",
    name = "Kollarovic2016 - Cell fate decision at G1-S transition",
    d = 8, n = 14,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 2, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      2, 0, 0, 0, 0, 0, 0, 0,
      2, 2, 0, 0, 0, 0, 0, 0,
      0, 0, 2, 0, 0, 0, 0, 0,
      0, 0, 2, 1, 0, 0, 2, 0,
      0, 0, 0, 2, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 2, 0, 0, 0,
      0, 0, 0, 0, 0, 2, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 2, 0,
      0, 0, 0, 0, 0, 0, 0, 2
    ), nrow = 14, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 2, 0, 0, 2, 0,
      0, 0, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 2,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1
    ), nrow = 14, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000167",
    name = "Mayya2005_STATmodule",
    d = 8, n = 14,
    S = matrix(c(
      0, 0, 0, 2, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 1, 0, 0, 0,
      0, 0, 2, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0
    ), nrow = 14, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 1, 0,
      0, 0, 1, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 2, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1
    ), nrow = 14, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000301",
    name = "Friedland2009_Ara_RTC3_counter",
    d = 8, n = 15,
    S = matrix(c(
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0
    ), nrow = 15, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 1, 1,
      0, 1, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 15, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000084",
    name = "Hornberg2005_ERKcascade",
    d = 8, n = 16,
    S = matrix(c(
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0
    ), nrow = 16, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 0, 0,
      0, 0, 0, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1
    ), nrow = 16, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000502",
    name = "Messiha2013 - Pentose phosphate pathway model",
    d = 8, n = 16,
    S = matrix(c(
      0, 0, 1, 2, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 2, 0, 0,
      0, 0, 0, 0, 2, 1, 0, 0,
      0, 0, 0, 0, 0, 2, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 2,
      0, 2, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 2, 0,
      2, 0, 0, 0, 0, 0, 1, 0,
      2, 0, 0, 0, 2, 0, 2, 2,
      1, 0, 0, 0, 2, 0, 2, 1,
      2, 0, 0, 0, 2, 0, 1, 2,
      2, 0, 0, 0, 1, 0, 2, 1,
      0, 1, 1, 0, 0, 0, 0, 0,
      0, 0, 2, 0, 0, 0, 0, 0,
      2, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 2, 0, 0, 0
    ), nrow = 16, byrow = TRUE),
    P = matrix(c(
      0, 0, 2, 1, 0, 2, 0, 0,
      0, 0, 0, 0, 2, 1, 0, 0,
      0, 0, 0, 0, 1, 2, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 2,
      0, 0, 0, 0, 0, 2, 0, 1,
      0, 1, 0, 2, 0, 0, 0, 0,
      2, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 2, 0,
      1, 0, 0, 0, 2, 0, 2, 1,
      2, 0, 0, 0, 2, 0, 2, 2,
      2, 0, 0, 0, 1, 0, 2, 1,
      2, 0, 0, 0, 2, 0, 1, 2,
      0, 2, 2, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0
    ), nrow = 16, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000065",
    name = "Yildirim2003_Lac_Operon",
    d = 8, n = 16,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0
    ), nrow = 16, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 1, 0
    ), nrow = 16, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000822",
    name = "Dorvash2019 - Dynamic modeling of signal transduction by mTOR complexes in cancer",
    d = 8, n = 18,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 1
    ), nrow = 18, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0
    ), nrow = 18, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000495",
    name = "Sen2013 - Phospholipid Synthesis in P.knowlesi",
    d = 8, n = 18,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0,
      3, 0, 0, 0, 0, 0, 0, 0,
      3, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 3, 0, 0, 0, 0,
      0, 0, 3, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 3,
      0, 0, 3, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 3, 0,
      0, 0, 0, 0, 0, 3, 0, 0,
      0, 0, 0, 0, 0, 0, 3, 0,
      0, 0, 0, 0, 0, 3, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 3,
      0, 0, 0, 0, 0, 0, 1, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 3, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0, 0, 0,
      0, 0, 0, 3, 0, 0, 0, 0
    ), nrow = 18, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0, 0, 0,
      2, 0, 0, 0, 0, 0, 0, 1,
      2, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 2, 0, 0, 0, 0,
      0, 0, 2, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 2,
      0, 0, 2, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 2, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 2, 0,
      0, 0, 0, 0, 0, 2, 0, 1,
      0, 0, 0, 0, 0, 0, 2, 1,
      0, 0, 0, 0, 0, 2, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 2,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 2, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 3, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0, 0, 0
    ), nrow = 18, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000221",
    name = "Singh2006_TCA_Ecoli_acetate",
    d = 8, n = 22,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0
    ), nrow = 22, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0
    ), nrow = 22, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000222",
    name = "Singh2006_TCA_Ecoli_glucose",
    d = 8, n = 22,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0
    ), nrow = 22, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0
    ), nrow = 22, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000212",
    name = "Curien2009_Aspartate_Metabolism",
    d = 8, n = 25,
    S = matrix(c(
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1
    ), nrow = 25, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 2, 0, 0, 0,
      0, 0, 0, 0, 2, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 25, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000926",
    name = "Rhodes2019 - Immune-Mediated theory of Metastasis",
    d = 8, n = 26,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 26, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 1,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0
    ), nrow = 26, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000803",
    name = "Park2019 - IL7 receptor signaling in T cells",
    d = 9, n = 8,
    S = matrix(c(
      0, 0, 0, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 8, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 1, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 0, 0, 0, 0
    ), nrow = 8, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000844",
    name = "Viertel2019 - A Computational model of the mammalian external tufted cell",
    d = 9, n = 10,
    S = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 10, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000423",
    name = "Nyman2012_InsulinSignalling",
    d = 9, n = 11,
    S = matrix(c(
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 1, 0, 1, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1
    ), nrow = 11, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 1, 0, 1, 0, 1,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0
    ), nrow = 11, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000584",
    name = "Mandlik2015 - Tristable genetic circuit of Leishmania",
    d = 9, n = 14,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 2, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 2, 0, 0, 0,
      2, 2, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 2, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 2, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 2, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 2, 0, 0, 0, 0
    ), nrow = 14, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      1, 1, 1, 0, 0, 0, 0, 0, 0,
      1, 1, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 14, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000981",
    name = "Wan2020 - risk estimation and prediction of the transmission of COVID-19 in maninland China excluding Hubei province",
    d = 9, n = 15,
    S = matrix(c(
      1, 0, 0, 0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 2, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0
    ), nrow = 15, byrow = TRUE),
    P = matrix(c(
      1, 0, 1, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 2, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 2, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 2, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0
    ), nrow = 15, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000320",
    name = "Grange2001 - PK interaction of L-dopa and benserazide",
    d = 9, n = 16,
    S = matrix(c(
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0
    ), nrow = 16, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 16, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000249",
    name = "Restif2006 - Whooping cough",
    d = 9, n = 20,
    S = matrix(c(
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      1, 1, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 1, 1, 1, 0, 0, 0, 1,
      1, 1, 0, 0, 1, 0, 1, 0, 0,
      0, 0, 1, 1, 1, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0
    ), nrow = 20, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      2, 1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 2, 1, 1, 0, 0, 0, 0,
      1, 2, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 2, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1
    ), nrow = 20, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000620",
    name = "Palmer2014 - Effect of IL-1β-Blocking therapies in T2DM - Disease Condition",
    d = 9, n = 20,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 1, 0, 0,
      0, 0, 1, 1, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 0, 1, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0
    ), nrow = 20, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 1, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 1, 1, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 20, byrow = TRUE)
  ),
  list(
    source = "HILL", id = "BIOMD0000000621",
    name = "Palmer2014 - Effect of IL-1β-Blocking therapies in T2DM - Healthy Condition",
    d = 9, n = 20,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 1, 0, 0,
      0, 0, 1, 1, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 0, 1, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0
    ), nrow = 20, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 1, 0, 1, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 1, 1, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 20, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000126",
    name = "Clancy2002_CardiacSodiumChannel_WT",
    d = 9, n = 22,
    S = matrix(c(
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 22, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1
    ), nrow = 22, byrow = TRUE)
  ),
  list(
    source = "HARV", id = "BIOMD0000000696",
    name = "Boada2016 - Incoherent type 1 feed-forward loop (I1-FFL)",
    d = 9, n = 23,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0, 0
    ), nrow = 23, byrow = TRUE),
    P = matrix(c(
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 1
    ), nrow = 23, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000107",
    name = "Novak1993 - Cell cycle M-phase control",
    d = 9, n = 23,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 23, byrow = TRUE),
    P = matrix(c(
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 23, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000218",
    name = "Singh2006_TCA_mtu_model2",
    d = 9, n = 28,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 28, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 28, byrow = TRUE)
  ),
  list(
    source = "MM", id = "BIOMD0000000219",
    name = "Singh2006_TCA_mtu_model1",
    d = 9, n = 30,
    S = matrix(c(
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 30, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0, 0
    ), nrow = 30, byrow = TRUE)
  ),
  list(
    source = "BM", id = "BIOMD0000000622",
    name = "NguyenLK2011 - Ubiquitination dynamics in Ring1B/Bmi1 system",
    d = 10, n = 19,
    S = matrix(c(
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
      1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 1, 0, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 19, byrow = TRUE),
    P = matrix(c(
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 1, 1, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 19, byrow = TRUE)
  )
)
