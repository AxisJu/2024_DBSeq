# Computation on supplied in-memory inputs.
compute <- function(supplied_inputs) {
  suppressPackageStartupMessages({
      library(ggplot2)
      library(dplyr)
      library(RColorBrewer)
  })
  colors_celltype_level3 <- c(Astro = "#a58946", EPCs = "#594a26", OPC = "#5953ff", Oligo = "#201e5a", 
      VCs = "#858881", CHPCs = "#cfd4c9", Micro = "#a87c5a", ExN_RSP_L23IT = "#ff1a71", ExN_RSP_L45IT = "#ed1986", 
      ExN_CLA = "#ba1369", ExN_CA3 = "#860d4c", ExN_CA1 = "#53082f", ExN_RSP_L6CT = "#61e2a4", ExN_HPF_CajalRetzius = "#d00000", 
      ExN_DG = "#16f2f2", InN_Immature = "#1b4332", InN_RSP_MGE = "#f954ee", InN_STRns_MGE = "#70e000", 
      InN_GPe_MGE = "#56ad00", InN_GPi_MGE = "#2f6000", InN_Nonspecific = "#f0a0ff", InN_STRns_LGE = "#d899ff", 
      InN_GPe_LGE = "#b199ff", InN_CP_D1 = "#8d7acc", InN_CP_D2 = "#695b99", InN_OT = "#4e4372", InN_LSc = "#3283fe", 
      InN_HYa = "#450099", InN_RT = "#ff6600", ExN_HYa = "#aa0dfe", ExN_TRS = "#7609b1", ExN_BST = "#420564", 
      ExN_MH = "#faa307", ExN_LH = "#c68105", ExN_THns = "#0d47a1", ExN_ATN = "#1460ff", ExN_RE = "#1340ff", 
      ExN_PF = "#0a4093", ExN_CM = "#08306d", InN_SCsg = "#9ef01a")
  global_ct_colors <- c(Glu.N = "#ba1369", GABA.N = "#56ad00", Astro. = "#a58946", Astro = "#a58946", OPCs = "#5953ff", 
      OPC = "#5953ff", Oligo. = "#201e5a", Oligo = "#201e5a", Endo. = "#858881", VCs = "#858881", Microglia = "#a87c5a", 
      Micro = "#a87c5a", `Mo/M<U+03C6>` = "#825f45", Epi. = "#594a26", EPCs = "#594a26", CHPCs = "#cfd4c9", 
      ExN = "#ba1369", InN = "#56ad00")
  colors_group <- c(DBS_I = "#ffd700", DBS_C = "#fff3b2", DBS = "#ffd700", Sham_I = "#ac00cc", Sham = "#ac00cc", 
      PTZ = "#eb7fff", PTZ_M1 = "#eb7fff", Saline_I = "#18a799", Saline = "#18a799", Control = "#18a799", 
      Epilepsy = "#ac00cc")
  colors_donors <- c(Saline_M1 = "#2a9d8f", Saline_M2 = "#18a799", Saline_M3 = "#13756b", Sham_M1 = "#c77dff", 
      Sham_M2 = "#ac00cc", Sham_M3 = "#7b0091", DBS_M1 = "#ffe45e", DBS_M2 = "#ffd700", DBS_M3 = "#ffb703", 
      DBS_M4 = "#fb8500", DBS_M5 = "#e09f3e", PTZ_M1 = "#eb7fff")
  region_colors <- c(RSPd = "#ff1a71", CA3 = "#860d4c", CA1 = "#53082f", DG = "#16f2f2", GPe = "#b199ff", 
      CP = "#8d7acc", LSc = "#3283fe", RT = "#ff6600", TRS = "#7609b1", MH = "#faa307", LH = "#c68105", 
      ATN = "#1460ff", RE = "#1340ff", PF = "#0a4093", CM = "#08306d")
  colors_brainregion <- c(HY = "#e74538", CNU = "#9bd5f4", Isocortex = "#18a799", HPF = "#84c551", TH = "#e0abce", 
      Midbrain = "#fcc396", `Non-neuronal` = "#000000", Nonspecific = "#cccccc")
  pal_jaccard <- c("#f4f9f4", "#7fcdbb", "#2c7fb8", "#253494")
  pal_cohend <- c(rev(brewer.pal(9, "Blues")[3:8]), "white", brewer.pal(9, "Reds")[3:8])
  invisible(as.list(environment()))
}
