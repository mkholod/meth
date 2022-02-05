library("DESeq2")

if(species == "hs"){
  mapTab <-  mapTab1 <- readRDS("mapTab_HS.rds")
  mapTab2 <- readRDS("mapTab_HS_entrez.rds")
}

if(species == "mm"){
  mapTab <- mapTab1 <- readRDS("mapTab_mm.rds")
  mapTab2 <- readRDS("mapTab_entrezgene_id.rds")
  
}
    