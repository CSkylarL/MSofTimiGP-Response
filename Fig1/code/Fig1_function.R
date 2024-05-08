#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Li
# 10/06/2023
# Analyze spatial data
# Here are functions to analyze the IMC data modified from 
# https://zenodo.org/record/7990870 (AUTHOR: H R Ali)
# major modification:
  # Add `dir` parameter to specify the location of the data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From Header.R ################################################################

setDTthreads(detectCores())
projectSeed <- 343418784
set.seed(projectSeed)

# From Header.R ################################################################
# 1 ----------------------------------------------------------------------------
#' Capitalise first letter; similar to Stata's -proper-
#' @param string 
#' @return string in proper format
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  out <- paste(toupper(substring(c, 1,1)), substring(c, 2),
               sep="", collapse=" ")
  return(out)
}


#' shorthand for summary
su <- function(x, ...){summary(x, ...)}


#' shorthand for as.data.table
adt <- function(x, ...){as.data.table(x, ...)}


#' shorthand for sort unique elements
sortunique <- function(x, ...){sort(unique(x,...))}

#' Gets antibody panel
#' 
#' Returns formatted makers, named by metal
#' @return character vector of formatted names, named by metal tag
getPanel <- function(dir){
  panelDir <- paste0(dir, '/data/raw/panel')
  metalReadOrder <- fs::path(panelDir, 'metalReadOrder.csv')
  metalReadOrder <- readLines(metalReadOrder)
  panel <- fread(fs::path(panelDir, 'NeoTripFinalPanelToIMCTools1.csv'))
  replaceWith <- function(pattern, replacement){
    panel[grep(pattern, Target), Target := replacement]
  }
  replaceWith('Histone', 'H3')
  replaceWith('Cytokeratin 5/14', 'CK5/14')
  replaceWith('8/18', 'CK8/18')
  replaceWith('H2AX', 'pH2AX')  
  replaceWith('Rabbit', 'PD-L1 (SP142)')
  replaceWith('CD134', 'OX40')
  replaceWith('GATA-3', 'GATA3')
  replaceWith('TOX/TOX2', 'TOX')
  replaceWith('CD279', 'PD-1')
  replaceWith('FoxP3', 'FOXP3') 
  replaceWith('CD274', 'PD-L1 (73-10)')
  replaceWith('CD278', 'ICOS')
  replaceWith('PDGFRbeta', 'PDGFRB')
  replaceWith('Granzyme', 'GZMB')
  replaceWith('TCF1', 'TCF1')
  replaceWith('Podoplanin', 'PDPN')
  replaceWith('cleaved', 'c-PARP')
  replaceWith('panKeratin', 'panCK')
  replaceWith('EMPTY', 'Carboplatin')
  replaceWith('Ki-67', 'Ki67')
  replaceWith('Calponin 1', 'Calponin')
  panel <- panel[!duplicated(Target)]
  panelOut <- panel[, Target]
  names(panelOut) <- panel[, `Metal Tag`]
  panelOut <- panelOut[metalReadOrder]
  cat('Panel in imctools read order\n')
  
  noDNA <- grep('DNA', panelOut, invert = TRUE, value = TRUE)
  epithelial <- grep('SMA|CD68|CD45|CD20|CD163|CD8|OX40|CD11c|CD3|TOX|T-bet|PD-1|FOXP3|ICOS|CD4|PDGFRB|CD31|GZMB|PDPN|CD79a|MPO|Calponin|Caveolin',
                     noDNA, invert = TRUE, value = TRUE)
  TME <- grep('CK|panCK|^AR$', noDNA, invert = TRUE, value = TRUE)
  
  return(list(panel = panelOut, Epithelial = epithelial, TME = TME))
}


#'Load data.table of IDs
#'
#'@return data.table of all IDs 
getIDs <- function(dir){
  IDs <- read_fst(paste0(dir, '/data/derived/IDs.fst'), as.data.table = T)
}


#'Load data.table of single-cell data
#'
#'@param curated logical, defualt is FALSE, if TRUE loads cells with complete curation columns including cluster labels, colours and print order. Also Vimentin and Calponin counts corrected for Carboplatin. 
#'@param allCells logical, whether to include non-invasive epithelial cells in curated data 
#'@param ... passed to read_fst; useful to constrain data using -columns- arg
#'@return data.table of single-cell or summary data 
getCells <- function(...,dir, curated = FALSE, allCells = FALSE){
  cells <- read_fst(paste0(dir, '/data/derived/cells.fst'), as.data.table = TRUE, ...)
  if(allCells){
    cells <- cells[!is.na(CellClusters)]
    return(cells)
  }
  if(curated){ 
    cells <- cells[!is.na(CellClusters)][cellAnnotation %in% c('invasive','TME')] 
    return(cells)  
  }  
  return(cells)
}


#'Return a names vector of quantile thresholds for positivity based on image inspection using getCells(curated = TRUE, allCells = TRUE)
#'@return a named vector of quantile threshold for each marker
getPositiveThresholds <- function(){
  positiveQuantile <- c( 
    'Ki67' = 0.83,
    'CA9' = 0.88,
    'pH2AX' = 0.92,
    'CD163' = 0.87,
    'H3' = 0.1,
    'CD20' = 0.92,
    'CD56' = 0.97,
    'Helios' = 0.83,
    'CD8' = 0.81,
    'OX40' = 0.99,
    'CD11c' = 0.81,
    'CD3' = 0.71,
    'GATA3' = 0.81,
    'SMA' = 0.89,
    'TOX' = 0.94,
    'IDO' = 0.975,
    'AR' = 0.90,
    'FOXP3' = 0.93,
    'ICOS' = 0.93,
    'CD4' = 0.72,
    'TCF1' = 0.77,
    'PDGFRB' = 0.84,
    'CD31' = 0.88,
    'GZMB' = 0.93,
    'PDPN' = 0.87,
    'CD79a' = 0.85,
    'Carboplatin' = 0.87,
    'Vimentin' = 0.8,
    'Calponin' = 0.93,
    'CD15' = 0.93,
    'MPO' = 0.99,
    'CD68' = 0.87,
    'CD45' = 0.63,
    'CK5/14' = 0.87,
    'CK8/18' = 0.88,
    'T-bet' = 0.87,
    'PD-L1 (73-10)' = 0.9,
    'PD-L1 (SP142)' = 0.98,
    'c-PARP' = 0.985,
    'HLA-ABC' = 0.45,
    'Caveolin-1' = 0.9,
    'PD-1' = 0.89,
    'HLA-DR' = 0.7)
  return(positiveQuantile)
}

#' Make a data.table of all combinations of values from variables within a data.table
#' 
#' @param dat, data.table that contains the variables
#' @param vars, character vector of variable names to cross for all combinations
#' @return data.table of all combinations of \code{vars}
mkAllCombinations <- function(dat, vars){
  getUnique <- function(x) {
    out <- dat[, unique(.SD), .SDcols = x]  
    return(out[, get(x)])
  }
  toExpand <- lapply(vars, getUnique)
  allCombinations <- adt(expand.grid(toExpand))   
  setnames(allCombinations, names(allCombinations), vars)
  return(allCombinations)
}
# 2 ----------------------------------------------------------------------------


#'Load data.table of formatted clinical data.
#'   
#'@return data.table of clinical data (randomisation and response)
getClinical <- function(dir){
  clinical <- read_fst(paste0(dir, '/data/derived/clinical.fst'), as.data.table = TRUE)
  return(clinical)
}  


#'Get annotations for cell clusters
#'
#'@return 3 element list of cell annotations (order, colours and labels), separately for Epithelial and TME, and a summary data.table for 'All'
getCellClusters <- function(dir){
  datDir <- paste0(dir, '/data/derived/') 
  getMetaDat <- function(type) {
    file <- paste0(type, 'Clusters.csv')
    Clusters <- fread(file.path(datDir, file))
    setkey(Clusters, 'PrintOrder')
    Clusters <- Clusters[, .SD, .SDcols = c("ClusterID", "Label", "Colour", "PrintOrder")]
    labels <- setNames(Clusters[, unique(Label)], Clusters[, unique(PrintOrder)])
    cols <- setNames(Clusters[, unique(Colour)], Clusters[, unique(PrintOrder)])
    return(list(dat = Clusters, labels = labels, colours = cols))
  }
  types <- c('Ep', 'TME')
  out <- lapply(types, getMetaDat)
  names(out) <- types
  
  All <- rbind(out$Ep$dat[, isEpithelial := TRUE], out$TME$dat[, isEpithelial := FALSE])
  All[, keep := seq_len(.N) == 1L, by = .(PrintOrder, isEpithelial)]
  All <- All[(keep)]
  All <- All[, ClusterID := NULL][, keep := NULL]
  All[!(isEpithelial), isImmune := !grepl('Endothelial|Myo|Fibro|Stromal|CA9', Label)]
  All[(isEpithelial), Type := 'Epithelial']
  All[(isImmune), Type := 'Immune']
  All[is.na(Type), Type := 'Stromal']
  All[, isImmune := NULL]
  out[['All']] <- All
  return(out)
}


#'Get cell phenotype counts 
#'@return list of two items: 'Counts' is a data.table of cell phenotype counts aggregated per patient (always separately by BiopsyPhase and by epithelial v tme cells). And 'CountsPerImage', which is the same but aggregated over each image rather than patient. Where a cell type is not present the count is 0; one exception is where there is no epithelium at all, in which case no counts are returned for epithelial phenotypes. 
getCellCounts <- function(dir){
  # Complex function; challenge is to retrieve zero counts per image or patient for each cell, for each biopsy phase
  # While returning no counts for epithelial cells if they are completely absent from the image or all images for a patient (per biopsy phase)
  panel <- getPanel(dir = dir)
  cellAnnot <- getCellClusters(dir = dir)[['All']]
  mainCells <- getCells(dir = dir, curated = TRUE)
  mainCells <- mainCells[, .SD, .SDcols = setdiff(names(mainCells), c(panel$panel))]
  
  mkCountsBy <- function(ByVar){
    cells <- copy(mainCells)
    cells <- cells[, .SD, .SDcols = grep('Parent|Location|Child|PrintOrder|Colour|Area|BiopsyType|CellClusters', 
                                         names(cells), invert = TRUE, value = TRUE)]
    
    setkeyv(cells, c(ByVar, 'BiopsyPhase', 'isEpithelial'))
    cells[, hasEpithelium := sum(as.integer(isEpithelial)), by = c(ByVar, 'BiopsyPhase')]
    ByVars <- c(ByVar, 'BiopsyPhase', 'isEpithelial')
    cells[, TotalCells := .N, by = ByVars]
    cells[, CellCount := .N, by = c(ByVars, 'Label')]
    cells <- cells[, .SD[1], by = c(ByVars, 'Label')]
    
    allCombinations <- mkAllCombinations(cells, vars = c('Label', ByVar, 'BiopsyPhase'))
    cells <- merge(x = allCombinations, y = cells, by = names(allCombinations), all = TRUE)  
    cells[, isEpithelial := NULL]
    cells <- merge(x = cells, y = cellAnnot, by = 'Label', all = TRUE) 
    cells[, allNA := all(is.na(CellCount)), by = ByVars] # no sample at a phase
    cells <- cells[!(allNA)][, allNA := NULL]
    setkeyv(cells, c(ByVar, 'BiopsyPhase', 'isEpithelial', 'PrintOrder'))
    
    fillinVar <- function(var, byVars) {
      cells[order(get(var), decreasing = TRUE), eval(var) := get(var)[1], by = byVars]
    }
    
    fillinVar('hasEpithelium', byVars = c(ByVar, 'BiopsyPhase'))
    cells[hasEpithelium == 0L & isEpithelial == FALSE, hasEpithelium := 1L] # retain TME cells even if no cancer cells
    cells <- cells[hasEpithelium != 0L] # drop epithelial phenotypes in samples lacking any cancer cells 
    sapply(c(ByVar, 'BiopsyPhase', 'TotalCells'), fillinVar, by = ByVars)
    if(ByVar == 'ImageID') sapply(c('ImageNumber', 'PatientID'), fillinVar, by = ByVars)
    
    cells[is.na(CellCount), CellCount := 0L]
    minmax <- cells[, sum(CellCount / TotalCells), by = ByVars][, range(V1)]
    isOne <- function(x) x > 0.99 & x < 1.01
    stopifnot(all(sapply(minmax, isOne)))
    
    cells[, c('ObjectNumber', 'hasEpithelium')  := NULL]
    if(ByVar == 'PatientID') cells[, c('ImageNumber', 'ImageID') := NULL]
    setkeyv(cells, c(ByVar, 'BiopsyPhase', 'isEpithelial', 'PrintOrder'))
    suffix <- ifelse(ByVar == 'PatientID', 'PerPatient', 'PerImage')
    oldNames <- c('CellCount', 'TotalCells')
    newNames <- paste0(oldNames, suffix)
    setnames(cells, oldNames, newNames)
    return(cells)
  }
  
  out <- lapply(c('PatientID', 'ImageID'), mkCountsBy)
  names(out) <- c('Counts', 'CountsPerImage')  
  cat('Note : totals are computed separately for epithelial and TME cells\n')
  return(out)
}


#' Get cell phenotype densities by image area
#' 
#' @param wide, boolean: whether to return the data in wide format, defaults to FALSE
#' @param perImage, boolean: whether to return data per image (rather than per tumour), defaults to FALSE
#' @return data.table of cell phenotype densities
getDensities <- function(dir, wide = FALSE, perImage = FALSE){
  convexHullArea <- function(xCoords, yCoords){
    coordDat <- as.matrix(cbind(xCoords, yCoords))
    ChullOrder <- chull(x = xCoords, y = yCoords)
    polgyonOrder <- coordDat[ChullOrder,]
    polgyonOrder <- rbind(polgyonOrder, polgyonOrder[1,])
    polygonArea <- sp::Polygon(polgyonOrder, hole = F)@area
    return(polygonArea)
  }
  cells <- getCells(dir = dir, columns = c('ImageID', 'ObjectNumber', 
                                'Location_Center_X', 'Location_Center_Y'))
  IDs <- getCells(dir = dir, curated = TRUE, allCells = TRUE)[, unique(ImageID)]
  cells <- cells[ImageID %in% IDs]
  cells[, totalArea := 
          convexHullArea(Location_Center_X, Location_Center_Y), 
        by = ImageID]
  area <- cells[, .SD[1], by = ImageID]
  IDs <- getIDs(dir = dir)[, .(ImageID, PatientID, BiopsyPhase)]
  area <- merge(x = area, y = IDs, by = 'ImageID') 
  area[, totalArea := totalArea / 1e6] # in mm^2
  
  idVars <- if(perImage) c('ImageID', 'BiopsyPhase') else c('PatientID', 'BiopsyPhase')
  area[, totalArea := sum(totalArea), by = idVars]
  area <- area[, .SD[1], by = idVars]
  area <- area[, .SD, .SDcols = c(idVars, 'totalArea')]
  
  counts <- if(perImage) getCellCounts(dir = dir)$CountsPerImage else getCellCounts(dir = dir)$Counts
  counts <- merge(x = counts, y = area, by = idVars) 
  if(perImage) counts[, Density := CellCountPerImage / totalArea]
  else counts[, Density := CellCountPerPatient / totalArea]
  
  if(wide) counts <- dcast(counts, as.formula(paste(paste(idVars, collapse = ' + '), 'Label', sep = '~')), 
                           value.var = 'Density')
  cat('Densities based on area of polygon drawn using convex hull of all segmented cells (per mm^2)\n')
  return(counts)
}
# 3-----------------------------------------------------------------------------
#'Clip vector at 99th centile if it's >0L
#'@param numeric vector to clip
#'@param centile number between 0 and 1 to pass to \code{quantile}; default is 0.99 (99th centile)
#'@return numeric vector with values exceeding 99th centile clipped to 99th centile
clip99 <- function(vec, centile = 0.99){
  c99 <- quantile(vec, probs = centile)
  if(isTRUE(c99 > 0L)) vec[vec > c99] <- c99 
  return(vec)
}
