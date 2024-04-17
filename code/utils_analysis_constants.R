# markers_all = list(
#   hepatocyte = c("Alb", "Ttr", "Apoa1", "Serpina1c", "Fabp1", "Echs1", "Glul", "Acly", "Asl", "Cyp2e1", "Cyp2f2", "Ass1", "Mup3", "Pck1", "G6pc", "Apoa1", "Ass1", "G6pc", "Mup3"),
#   hybrid_hepatocyte = c("Meg3", "Igfbp3", "Trp53inp1"),
#   midzonal_hepatocyte = c("Cyp2e1", "Cyp1a2", "Alb"),
#   periportal_hepatocyte = c("Apoa1", "Apoa2", "Apoa5", "Apob", "Alb", "Cyp2f2"),
#   duct_endothelial_cell = c("Sprr1a", "Krt7", "Fxyd3", "Krt15", "Upk1a", "Ly6d", "S100a6", "Cldn4", "Krt19"),
#   endothelial_activated = c("Ptprb", "Vwf", "Klf4", "Sdc1"),
#   endothelial_cell_of_hepatic_sinusoid = c("Eng", "Ehd3", "Cd300lg", "Ramp2", "Fam167b", "Pecam1", "Oit3", "Esam", "Sema6a", "Kdr", "Adam23", "Fcgr2b", "Gpr182"),
#   vascular_endothelial_cell = c("Lrrc8a", "Pecam1", "Kdr", "Tek", "Flt1", "Flt4", "Nrp1", "Nrp2", "Eng"),
#   epiblast = c("Utf1", "Epacm", "Pou5f1", "Dnmt3b"),
#   exE_endoderm = c("Ttr", "Apoa2", "Apoe", "Cystm1", "Emb"),
#   hepatoblast = c("Id3", "Mdk", "Gpc3", "Dlk1", "Afp", "Alb", "Sox9", "Sox11", "Hnf4a", "Krt18", "Krt8", "Hnf1b", "Hhex", "Met", "Anpep", "Cdh1"),
#   hepatoblast_to_hepatocyte = c("Ppara", "Rora", "Thrb", "Cux2", "Nr1i3", "Esr1", "Nr1h4", "Ahr", "Stat6", "Foxq1", "Nfia", "Zbtb16", "Сrebl2", "Bhlhe40", "Zbtb20", "Tox", "Nfib", "Zip791", "Tfcp2l1", "Zhx3", "Klf9"),
#   liver_progenitor_cell = c("Ascl2", "Krt19", "Dlk1", "Epcam"),
#   primitive_Streak = c("Eomes", "Pou5f1", "Nanog", "Cacna1a"),
#   stem_progenitor_cell = c("Cd34", "Cmtm7", "Spi1", "Klf4", "Prom1"),
#   B_cell = c('Cd79a', "Cd79b", "Cd74", "Cd19", "Fcmr", "Jchain", "Mzb1", "Igkc", "Cd22", "Ebf1", "Pax5", "Cxcl10", "Vpreb3", "Vpreb1"),
#   dendritic_cell = c("Irf5", "Irf8", "Ccl22", "Il12b", "Ccr7", "Id2", "Id3", "Siglech", "Ly6d", "Bst2", "Itgax", "Cd80", "Cd83"),
#   erythroblast = c("Hbb-bs", "Hba-a2", "Hba-a1", "Klf1", "Alad", "Blvrb", "Hmbs", "Rhd", "Gata1"),
#   gammaDelta_T_cell = c("Cxcr6", "Cd3g", "Cd163l1", "Itgae"),
#   granulocyte = c("S100a9", "S100a8", "Ccl4", "Prtn3", "Slpi", "Fcer1g"),
#   Kuppfer_cell = c("Adgre1", "Emr1", "Clec4f", "Cd68", "Irf7", "C1qa", "C1qb", "C1qc", "Csf1r"),
#   liver_capsule_macrophages = c("S100a4", "Itgax", "Crip1", "Ccr2", "Naaa", "Cx3cr1"),
#   macrophage = c("Ccl9", "Cd14", "Marco", "Adgre1", "Cd68", "Csf1r", "Ptprc", "Cd52", "H2-Aa"),
#   monocytes = c("Itgam", "Ly6g", "Ly6c1", "Ccr2"),
#   natural_killer_cell = c("Zap70", "Il2rb", "Nkg7", "Cxcr6", "Itga1"),
#   neutrophil_proinflammatory = c("Ccl3", "Ccl4", "Cxcl2", "Csf1", "Adgre1", "Cd5l", "Clec4f", "Timd4", "Folr2", "C1qa", "C1qb", "C1qc", "Vsig4", "Xcr1", "Cd209a", "Siglech", "Chil3", "F13a1", "S100a4", "Lgals3", "Gda"),
#   neutrophils = c("Csf3r", "Sepx1", "Retnlg", "S100a8", "S100a9", "Slpi", "Fos", "Junb", "Jun", "Zfp36", "Itgam", "Ly6g", "Elane", "Mpo", "Cebpe", "Cst7", "Ngp", "Lcn2", "Apoe", "Cd5l", "Timd4", "Folr2", "C1qa", "C1qb", "C1qc", "Vsig4", "Xcr1", "Cd209a", "Siglech", "Chil3", "F13a1", "S100a4", "Lgals3", "Gda"),
#   pDCs = c("Bcl11a", "Runx2", "Ccr9", "Siglech", "Spib", "Irf8"),
#   T_cell = c("Ccl5", "Cd7", "Cd3g", "Trbc1", "Trbc2", "Thy1", "Lat", "Cd3d", "Cdse", "Nkg7", "Ptprc", "Cd68", "Cd52"),
#   hepatic_stellate_cell = c("Ecm1", "Colec11", "Vipr1", "Hgf", "Rgs5", "Lrat", "Ngfr", "Reln", "Pth1r", "Col1a1", "Acta2", "Col3a1", "Dcn"),
#   hepatocellular_carcinoma = c("Trf", "Serpina1a", "Orm1", "Hnf4a", "Fbp1", "Mat1a", "Sult1a1", "Apoa4", "Apob", "Cyp1a2", "Cyp3a11", "Cyp3a25", "Ugt1a1", "Fmo5", "Scd1", "Afp", "Gpc3"),
#   stellate_activated = c("Ecm1", "Ccl2", "Col3a1", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"),
#   stellate_activated_MYCi = c("Ecm1", "Trp53inp1", "Sfpq", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"),
#   stellate_fibrotic = c("Ecm1", "Lrat", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"),
#   stellate_quiescent = c("Ecm1", "Tagln", "Lrat", "Pdgfra", "Pdgfrb", "Rgs5", "Col14a1", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"),
#   biliary_epithelial_cell = c("Epcam", "Krt19", "Spp1", "Hnf1b", "Prom1", "St14", "Foxj1", "Cftr", "Sctr"),
#   cholangiocyte = c("Car2", "Cd44", "Bcl11a", "Epcam", "Krt19", "Sox9", "Krt7", "Spp1"),
#   cycling_cell = c("Foxm1", "Ccna2", "Ccnb1", "Ccnb2", "Ccne2", "Cdk1", "Stmn1"),
#   epithelial_cell = c("Krt8", "Krt19", "Mmp7", "Krt18", "Muc1", "Krt23", "Epcam", "Cldn6", "Cldn7", "Cdh1"),
#   fibroblast = c("Gsn", "Clec3b", "Dpt", "Cd34", "Mfap4", "Entpd2", "Fbln2", "Col15a1", "Ccnb2", "Apoa2", "C3"),
#   mast_cell = c("Cpa3", "Hdc", "Srgn", "Tpsb2", "Poln", "Pga5", "Fosb", "Cma1"),
#   megakaryocyte = c("Pf4", "Ppbp", "Cd9", "Itga2b", "Plek", "Cxcr4", "Itgb3"),
#   mesenchymal_cell = c("Vim", "Col1a2", "Mest", "Mmp2", "Pdgfra", "Ncam1", "Lhx2", "Fn1", "Cdh2", "Sparc"),
#   vascular_smooth_muscle_cell = c("Tagln", "Tpm2", "Acta2", "Cnn1", "Pln", "Myh11")
# )


markers_all = list(
    hepatocyte= c("Alb", "Ttr", "Apoa1", "Serpina1c", "Fabp1", "Echs1", "Glul", "Acly", "Asl", "Cyp2e1", "Cyp2f2", "Ass1", "Mup3", "Pck1", "G6pc", "Apoa1", "Ass1", "G6pc", "Mup3"),
    hybrid_hepatocyte= c("Meg3", "Igfbp3", "Trp53inp1"),
    midzonal_hepatocyte= c("Oat", "Cyp2e1", "Cyp1a2", "Alb"),
    periportal_hepatocyte= c("Sds", "Sdsl", "Hal", "Apoa1", "Apoa2", "Apoa5", "Apob", "Alb", "Cyp2f2"),
    hepatic_stellate_cell= c("Ecm1", "Colec11", "Vipr1", "Hgf", "Rgs5", "Lrat", "Ngfr", "Reln", "Pth1r", "Col1a1", "Acta2", "Col3a1", "Dcn"),
    hepatocellular_carcinoma= c("Trf", "Serpina1a", "Orm1", "Hnf4a", "Fbp1", "Mat1a", "Sult1a1", "Apoa4", "Apob", "Cyp1a2", "Cyp3a11", "Cyp3a25", "Ugt1a1", "Fmo5", "Scd1", "Afp", "Gpc3"),
    stellate_activated= c("Ecm1", "Ccl2", "Col3a1", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"),
    stellate_activated_MYCi= c("Ecm1", "Trp53inp1", "Sfpq", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"),
    stellate_fibrotic= c("Ecm1", "Lrat", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"),
    stellate_quiescent= c("Ecm1", "Tagln", "Lrat", "Pdgfra", "Pdgfrb", "Rgs5", "Col14a1", "Col1a1", "Ereg", "Acta2", "Dcn", "Colec11", "Cxcl12", "Sod3", "Angptl6", "Rgs5", "Reln", "Tmem56", "Rbp1", "G0s2", "Rarres2", "Acta2", "Tagln"),
    Endothelial= c(
        "Ptprb", "Vwf", "Klf4", "Sdc1", "Eng", "Ehd3", "Cd300lg", "Ramp2", "Fam167b", "Pecam1", "Oit3", "Esam",
        "Sema6a", "Kdr", "Adam23", "Fcgr2b", "Gpr182", "Lrrc8a", "Tek", "Flt1", "Flt4", "Nrp1", "Nrp2"
    ), 
    B_Cell= c(
        "Cd79a", "Cd79b", "Cd74", "Cd19", "Fcmr", "Jchain", "Mzb1", "Igkc", "Cd22", "Ebf1", "Pax5", "Cxcl10",
        "Vpreb3", "Vpreb1"
    ), 
    Erythroblast= c(
        "Hbb-bs", "Hba-a2", "Hba-a1", "Klf1", "Alad", "Blvrb", "Hmbs", "Rhd", "Gata1"
    ), 
    Kupffer_Cell= c(
        "Adgre1", "Emr1", "Clec4f", "Cd68", "Irf7", "C1qa", "C1qb", "C1qc", "Csf1r"
    ), 
    Macrophage= c(
        "Ccl9", "Cd14", "Marco", "Adgre1", "Cd68", "Csf1r", "Ptprc", "Cd52", "H2-Aa"
    ), 
    Monocytes= c(
        "Itgam", "Ly6g", "Ly6c1", "Ccr2"
    ), 
    Natural_Killer= c(
        "Zap70", "Il2rb", "Nkg7", "Cxcr6", "Itga1"
    ), 
    T_Cell= c(
        "Ccl5", "Cd7", "Cd3g", "Trbc1", "Trbc2", "Thy1", "Lat", "Cd3d", "Cdse", "Nkg7", "Ptprc", "Cd68", "Cd52",
        "Cxcr6", "Cd163l1", "Itgae"
    ), 
    Plasmacytoid_DCs= c(
        "Bcl11a", "Runx2", "Ccr9", "Siglech", "Spib", "Irf8"
    ), 
    Conventional_DCs= c(
        "Irf5", "Irf8", "Ccl22", "Il12b", "Ccr7", "Id2", "Id3", "Siglech", "Ly6d", "Bst2", "Itgax", "Cd80", "Cd83"
    ), 
    Mesenchymal_Cells= c(
        "Ecm1", "Colec11", "Vim", "Col1a2", "Mest", "Mmp2", "Pdgfra", "Ncam1", "Lhx2", "Fn1", "Cdh2", "Sparc",
        "Tagln", "Tpm2", "Acta2", "Cnn1", "Pln", "Myh11"
    ), 
    Megakaryocytes= c(
        "Pf4", "Ppbp", "Cd9", "Itga2b", "Plek", "Cxcr4", "Itgb3"
    ), 
    Neutrophils= c('Sepx1', 'Retnlg', 'C1qa', 'Xcr1', 'Lcn2', 'Lgals3', 'Csf1', 'Ly6g', 'Csf3r', 'Slpi', 'Timd4',
        'Cd209a', 'Cebpe', 'Cxcl2', 'Ccl4', 'Elane', 'Siglech', 'Ngp', 'Jun', 'Junb', 'S100a4', 'F13a1', 'Mpo', 
        'S100a8', 'Cst7', 'Ccl3', 'Zfp36', 'Cd5l', 'Vsig4', 'Clec4f', 'Apoe', 'Gda', 'C1qc', 'Itgam', 'Chil3', 
        'S100a9', 'Folr2', 'Fos', 'Adgre1', 'C1qb'), 
        Cholangiocytes= c(
        "Epcam", "Krt19", "Spp1", "Hnf1b", "Prom1", "St14", "Foxj1", "Cftr", "Sctr", "Car2", "Cd44", "Sox9", "Krt7"
    ), 
    Fibroblast_like= c(
        "Gsn", "Clec3b", "Dpt", "Cd34", "Mfap4", "Entpd2", "Fbln2", "Col15a1", "Ccnb2", "Apoa2", "C3", "Cpa3",
        "Hdc", "Srgn", "Tpsb2", "Poln", "Pga5", "Fosb", "Cma1", "Tagln", "Tpm2", "Acta2", "Cnn1", "Pln", "Myh11"
    ), 
    Mast_Cells= c(
        "Cpa3", "Hdc", "Srgn", "Tpsb2", "Poln", "Pga5", "Fosb", "Cma1"
    ), 
    Epithelial_like= c(
        "Krt8", "Krt19", "Mmp7", "Krt18", "Muc1", "Krt23", "Epcam", "Cldn6", "Cldn7", "Cdh1", "Utf1", "Epacm", 
        "Pou5f1", "Dnmt3b"
    )
    # "Stellate Cells= c('Hgf', 'Reln', 'Mat1a', 'Cyp3a11', 'Fmo5', 'Pdgfra', 'Serpina1a', 'Hnf4a', 'Col1a1', 'G0s2', 'Rarres2', 
    #      'Dcn', 'Tagln', 'Ecm1', 'Sfpq', 'Ngfr', 'Pdgfrb', 'Rgs5', 'Ccl2', 'Sult1a1', 'Cxcl12', 'Pth1r', 'Cyp3a25', 
    #      'Lrat', 'Trf', 'Scd1', 'Col3a1', 'Ugt1a1', 'Tmem56', 'Angptl6', 'Gpc3', 'Fbp1', 'Colec11', 'Rbp1', 
    #      'Trp53inp1', 'Apob', 'Col14a1', 'Cyp1a2', 'Vipr1', 'Acta2', 'Orm1', 'Apoa4', 'Ereg', 'Sod3', 'Afp']
    # "Others= c(
    #     "Ccl5", "Cd7", "Cd3g", "Trbc1", "Trbc2", "Thy1", "Lat", "Cd3d", "Cdse", "Nkg7", "Ptprc", "Cd68", "Cd52",
    #     "Bcl11a", "Runx2", "Ccr9", "Siglech", "Spib", "Irf8", "Eomes", "Nanog", "Cacna1a", "Foxm1", "Ccna2", 
    #     "Ccnb1", "Ccnb2", "Ccne2", "Cdk1", "Stmn1", "S100a9", "S100a8", "Ccl4", "Prtn3", "Slpi", "Fcer1g", "Trf",
    #     "Serpina1a", "Orm1", "Hnf4a", "Fbp1", "Mat1a", "Sult1a1", "Apoa4", "Apob", "Cyp1a2", "Cyp3a11", "Cyp3a25",
    #     "Ugt1a1", "Fmo5", "Scd1", "Afp", "Gpc3", "Apoa5", "Sprr1a", "Fxyd3", "Krt15", "Upk1a", "S100a6", "Cldn4",
    #     "Ascl2", "Cmtm7", "Spi1", "Cd163l1", "Itgae", "Crip1", "Naaa", "Cx3cr1", "Csf3r", "Retnlg", "Fos", "Junb",
    #     "Jun", "Zfp36", "Elane", "Mpo", "Cebpe", "Cst7", "Ngp", "Lcn2", "Vipr1", "Hgf", "Rgs5", "Lrat", "Ngfr",
    #     "Reln", "Pth1r", "Col1a1", "Col3a1", "Dcn", "Ccl2", "Ereg", "Cxcl12", "Sod3", "Angptl6", "Tmem56", "Rbp1",
    #     "G0s2", "Rarres2", "Sfpq", "Pdgfrb", "Col14a1"
    # ]
)

type_cols = c(
  hepatocyte = "#FF5733",
  hybrid_hepatocyte = "#FF8247",
  midzonal_hepatocyte = "#FFA07A",
  periportal_hepatocyte = "#FFD700",
  endothelial_cell = "#32CD32",
  endothelial_activated = "#98FB98",
  endothelial_cell_of_hepatic_sinusoid = "#00CED1",
  vascular_endothelial_cell = "#1E90FF",
  epiblast = "#00BFFF",
  extra_embryonic_endoderm = "#0000FF",
  hepatoblast = "#8A2BE2",
  hepatoblast_to_hepatocyte = "#9370DB",
  liver_progenitor_cell = "#6A5ACD",
  primitive_streak = "#483D8B",
  stem_progenitor_cell = "#7B68EE",
  B_cell = "#9400D3",
  dendritic_cell = "#8B008B",
  erythroblast = "#800080",
  gamma_delta_T_cell = "#4B0082",
  granulocyte = "#9932CC",
  Kuppfer_cell = "#8A2BE2",
  liver_capsule_macrophages = "#9370DB",
  macrophage = "#6A5ACD",
  monocytes = "#483D8B",
  natural_killer_cell = "#7B68EE",
  neutrophil_proinflammatory = "#9400D3",
  neutrophils = "#8B008B",
  pDCs = "#800080",
  biliary_epithelial_cell = "#FF0000",
  cholangiocyte = "#FF4500",
  cycling_cell = "#FF6347",
  epithelial_cell = "#FF7F50",
  fibroblast = "#FFA500",
  mast_cell = "#FFD700",
  megakaryocyte = "#FFFF00",
  mesenchymal_cell = "#ADFF2F",
  vascular_smooth_muscle_cell = "#00FF00",
  stellate_fibrotic = "#00CED1",
  stellate_quiescent = "#1E90FF",
  others = "#8A2BE2"

)




