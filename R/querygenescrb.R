#' @title QueryGenes
#' @description It produces a data frame with the phenotypical terms from the genes selected and their entrez ids.
#' @param genes A vector with the genes symbol
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param organism A string parameter: "human" if you have human genes symbol and "mouse" if you have mouse genes symbol
#' @details This function has developed by Alejandro Cisterna Garcia as part of his PhD mentored by Juan Antonio Botia Blaya
#' @import readr
#' @import data.table
#' @import ggplot2
#' @import stats
#' @import plotly
#' @import ggpubr
#' @import dplyr
#' @import viridis
#' @import clusterProfiler
#' @import parallel
#' @import purrr
#' @import DT
#' @import Hmisc
#' @import pheatmap
#' @export
#'

QueryGenes = function(genes, database = "HPO", organism = "human"){
  
  if (organism == "human") {
    
    if (database == "UNIPROT") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/UNIPROTbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "CGI") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/CGIbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "PSYGENET") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/PSYGENETbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "ORPHANET") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/ORPHANETbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "CTD") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/CTD_humanbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "GENOMICS_ENGLAND") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/GENOMICS_ENGLANDbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "CLINGEN") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/clingenbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "DIS") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/disbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "HPO") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/hpobasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "CRB") {
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/cbrdatabaseok.csv"), header = T)
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(custom, genes, by = "entrez")
    }
    
    if (database == "MGD"  | database == "MOUSEDB") {
      # Read mouse data base
      musedata <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                               "/mgibasefinal.tsv"))
      
      names(musedata) <- c("entrez","symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info", "source"   )
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his MP
      onlyg <- inner_join(musedata, genes, by = "entrez")
      onlyg <- onlyg[,c(1,2,5,6)]
    }
    
    
  }
  
  
  if (organism == "mouse") {
    
    # genes = musedata[1:50,3]
    
    genes <- as.data.frame(genes)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)
    
    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genes = unique(onlygenes)
    
    if (database == "UNIPROT") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/UNIPROTbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "CGI") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/CGIbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "PSYGENET") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/PSYGENETbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "ORPHANET") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/ORPHANETbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "CTD") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/CTD_humanbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "GENOMICS_ENGLAND") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/GENOMICS_ENGLANDbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Symbol to entrez
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "CLINGEN") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/clingenbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "DIS") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/disbasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "HPO") {
      # Read human data base
      hpox <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/hpobasededatos.csv"))
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(hpox, genes, by = "entrez")
    }
    
    if (database == "CRB") {
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/cbrdatabaseok.csv"), header = T)
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his HP
      onlyg <- inner_join(custom, genes, by = "entrez")
    }
    
    if (database == "MGD"  | database == "MOUSEDB") {
      # Read mouse data base
      musedata <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                               "/mgibasefinal.tsv"))
      
      names(musedata) <- c("entrez","symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info", "source"   )
      
      # Convert vector to a dataframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # We will keep only with the genes in the query and his MP
      onlyg <- inner_join(musedata, genes, by = "entrez")
      onlyg <- onlyg[,c(1,2,5,6)]
    }
    
  }
  
  
  return(onlyg)
}

#' @title getdbnames
#' @description This function show the databases available in PhenoExam package.
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export
getdbnames = function(){

  sources <- c("HPO", "CRB", "MGD", "CLINGEN",
               "GENOMICS_ENGLAND", "CTD", "ORPHANET",
               "PSYGENET", "CGI", "UNIPROT")
  return(sources)
}


#' @title findoverlapgenes
#' @description This function get the genes that are linked to the phenotypes selected.
#' @param phenos A vector with the phenotypes id that you want.
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export

findoverlapgenes = function(phenos){

  # Un archivo base de datos con todas las DB
  sdb <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                      "/sdb.csv"), header = T)


  # Make an empty df with names
  rbjointri = data.frame(matrix(vector(), 0, 5,
                                dimnames=list(c(), c("entrez", "symbol",   "term_name", "term_id", "source"))),
                         stringsAsFactors=F)


  aparece <- length(phenos)
  # For each database name in a vector run its loadder
  for (ph in phenos) {

    fenos = dplyr::filter(sdb, sdb$term_id == ph)


    # join data frames generated
    rbjointri <- rbind(rbjointri, fenos)



  }

  genesaparecen <-   rbjointri %>% dplyr::group_by(symbol,entrez) %>% dplyr::tally()  %>%
    dplyr::filter(n >= aparece)

  genesoverlap <- unique(genesaparecen$symbol)

  return(genesoverlap)
}

#' @title ListGenesPheno
#' @description This function make a list with genes selected and his phenotypical terms
#' @param genes A vector with the genes symbol
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param organism A string parameter: "human" if you have human genes symbol and "mouse" if you have mouse genes symbol
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export

ListGenesPheno = function(genes, database = "HPO", organism = "human"){


  # We will keep only with the genes in the query and his HP
  onlyg <- QueryGenes(genes, database, organism)


    # We make a list with genes and his HP
    gbygen <-  split(onlyg$term_id , onlyg$symbol)



  return(gbygen)
}



#' @title FindPhenoFromGenes
#' @description Get the genes selected and only their phenotypical terms defined by the user.
#' @param genes A vector with genes symbol.
#' @param phenoid A vector with id terms to get.
#' @param database A string parameter with different options:  "HPO", "DIS", "CRB" and "MGD".
#' @param organism A string parameter: "human" if you have human genes symbol and "mouse" if you have mouse genes symbol
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export


FindPhenoFromGenes = function(genes, phenoid, database = "HPO", organism = "human"){

  # We will keep only with the genes in the query and his phenotypic terms
  onlyg <- QueryGenes(genes,database, organism )

  # Convert vector to a dataframe
  phenoid <- as.data.frame(phenoid)
  names(phenoid) <- "term_id"

  # We only keep genes and particular phenotypic term
  onlygenesandpehno <- inner_join(onlyg, phenoid, by = "term_id")

  return(onlygenesandpehno)
}



#' @title PhenoGeneNumber
#' @description This function counts the number of genes per phenotype term from the vector of genes.
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param genes A vector with human gene symbol.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export


PhenoGeneNumber = function(genes, database= "HPO", organism = "human"){
  
  if (organism == "human") {
    
    
    if (database == "UNIPROT") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/UNIPROTbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CGI") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/CGIbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "PSYGENET") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/PSYGENETbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "ORPHANET") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/ORPHANETbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CTD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/CTD_humanbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "GENOMICS_ENGLAND") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/GENOMICS_ENGLANDbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CLINGEN") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/clingenbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "DIS") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/disbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "HPO") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/hpobasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CRB") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/cbrdatabaseok.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
    }
    
    if (database == "MGD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "human_symbol"
      genes$HumanSymbol <- as.character(genes$human_symbol)
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$human_symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/mgibasefinal.tsv"), header = T, sep = "\t")
      names(custom) <- c("entrez","human_symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info" , "source"  )
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name ) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
    }
  }
  
  if (organism == "mouse") {
    
    genes <- as.data.frame(genes)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)
    
    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genes = unique(onlygenes)
    
    if (database == "UNIPROT") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/UNIPROTbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CGI") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/CGIbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "PSYGENET") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/PSYGENETbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "ORPHANET") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/ORPHANETbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CTD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/CTD_humanbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "GENOMICS_ENGLAND") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/GENOMICS_ENGLANDbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CLINGEN") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/clingenbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "DIS") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/disbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "HPO") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/hpobasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CRB") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/cbrdatabaseok.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
    }
    
    if (database == "MGD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "human_symbol"
      genes$HumanSymbol <- as.character(genes$human_symbol)
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$human_symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/mgibasefinal.tsv"), header = T, sep = "\t")
      names(custom) <- c("entrez","human_symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info" , "source"  )
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name ) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
    }
    
  }
  
  return(numportergenes)
}



#' @title PhenoGeneNumberOpt
#' @description This function counts the number of genes per phenotype term from the vector of genes.
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param genes A vector with human gene symbol.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export



PhenoGeneNumberOpt = function(genes, database= "HPO", organism = "human"){
  
  if (organism == "human") {
    
    
    
    if (database == "UNIPROT") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/UNIPROTbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CGI") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/CGIbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "PSYGENET") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/PSYGENETbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "ORPHANET") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/ORPHANETbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CTD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/CTD_humanbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "GENOMICS_ENGLAND") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/GENOMICS_ENGLANDbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CLINGEN") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/clingenbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "DIS") {
      
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/disbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "HPO") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/hpobasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
    }
    
    
    if (database == "CRB") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/cbrdatabaseok.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
    }
    
    if (database == "MGD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/mgibasefinal.tsv"), header = T, sep = "\t")
      names(custom) <- c("entrez","human_symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info" , "source"  )
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name ) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
    }
  }
  
  if (organism == "mouse") {
    
    # genes = musedata[1:50,3]
    
    genes <- as.data.frame(genes)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)
    
    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genes = unique(onlygenes)
    
    
    if (database == "UNIPROT") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/UNIPROTbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CGI") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/CGIbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "PSYGENET") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/PSYGENETbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "ORPHANET") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/ORPHANETbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CTD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/CTD_humanbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "GENOMICS_ENGLAND") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/GENOMICS_ENGLANDbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "CLINGEN") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/clingenbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "DIS") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/disbasededatos.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
      
    }
    
    if (database == "HPO") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/hpobasededatos.csv"), header = T)
      names(custom) <- c("entrez","symbol", "term_id", "term_name")
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
    }
    
    if (database == "CRB") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/cbrdatabaseok.csv"), header = T)
      names(custom) <- c("entrez","symbol", "term_name", "term_id")
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
    }
    
    if (database == "MGD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "human_symbol"
      genes$HumanSymbol <- as.character(genes$human_symbol)
      
      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$human_symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))
      
      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/mgibasefinal.tsv"), header = T, sep = "\t")
      names(custom) <- c("entrez","human_symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info"  , "source" )
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")
      
      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name ) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
      
    }
  }
  
  return(numportergenes)
}

#' @title RandomComparePheno
#' @description This function measures statistically significant phenotype similarities between gene sets and detects significant differential phenotypes for them.
#' @param geneset A vector with genes symbol.
#' @param genesetcompare A vector with gene symbol.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param p.value Cut of p value as they are relevant terms
#' @param nulltestnumber A number of random genes subset.
#' @param seed A seed.
#' @param url A logical parameter. T if you want interactive links
#' @param linux A logical parameter. T if you use linux SO, this is the default option
#' @param adjust A string parameter. You have two options: "bonferroni" is the default option and you can choose "fdr" for False Discovery Rate.
#' @param POR A string parameter. You have two options: Use the Jaccard index "jaccard" is the default option and you can choose the Forbes coefficient "forbes".
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export


RandomComparePheno = function(geneset, genesetcompare,  organism = "human", database = "HPO", p.value = 0.05, nulltestnumber = 100, seed = 1, conditional_case = F, url = F, linux=T,   adjust = "bonferroni", POR = "jaccard"){
  
  if (conditional_case == F) {
    
    # Get the Phenotype Enrichment Analysis from geneset
    referencia = PhenoEnrichGenes(genes=geneset, database, organism, url = url)
    enrichresult = as.data.frame(referencia$htmltabla)
    
    # Obtain the significative phenotypes terms from geneset
    if (adjust == "bonferroni") {
      relevantes <- dplyr::filter(enrichresult, enrichresult$bonferroni_pvalue <= p.value)
    }
    
    if (adjust == "fdr") {
      relevantes <- dplyr::filter(enrichresult, enrichresult$FDR_adjust <= p.value)
    }
    
    # Get the Phenotype Enrichment Analysis from genesetcompare
    targetenrich = PhenoEnrichGenes(genesetcompare, database, organism, url = url)
    targetenrich = as.data.frame(targetenrich$htmltabla)
    nphenotargetall <- nrow(targetenrich)
    
    # Obtain the significative phenotypes terms from genesetcompare
    if (adjust == "bonferroni") {
      targetrelevant <- dplyr::filter(targetenrich, targetenrich$bonferroni_pvalue <= p.value)
    }
    
    if (adjust == "fdr") {
      targetrelevant <- dplyr::filter(targetenrich, targetenrich$FDR_adjust <= p.value)
    }
    
    nphenotarget <- nrow(targetrelevant)
    
    # Obtain overlap
    phenounion <- full_join(enrichresult, targetenrich, by = "term_id")
    prunionn <- full_join(enrichresult,targetenrich, by = c("term_id","term_name","source"))
    nphenounion <- as.numeric(nrow(phenounion))
    
    # Inner phenotypes
    phenopresentall <- inner_join(enrichresult, targetenrich, by = "term_id") # Innner de fenotipos totales entre los dos datase
    ncompartidosall <- nrow(phenopresentall) # Número de fenotipos totales compartidos
    
    # RPOR interesccion entre union
    if (POR=="jaccard"){ 
      scoreshared <- ncompartidosall/nphenounion
    }
    
    if (POR=="forbes"){ 
      scoreshared <- forbes(enrichresult$term_id, targetenrich$term_id)
    }
    
    
    
    
    # Obtain union sifnificative
    if (adjust == "bonferroni") {
      set1join <- dplyr::filter(prunionn, prunionn$bonferroni_pvalue.x <= p.value)
    }
    
    if (adjust == "fdr") {  
      set1join <- dplyr::filter(prunionn, prunionn$FDR_adjust.x <= p.value)
    }
    
    set1join <-  set1join[,c(1:5,7:12,14:17)]
    names(set1join) <- c("term_id", "term_name","source", "bonferroni_pvalue_set1", "FDR_adjust_set1", "gene_overlap_set1", "overlap_ratio_set1", "pvalue_set1", "common_genes_set1","bonferroni_pvalue_set2","FDR_adjust_set2", "gene_overlap_set2", "overlap_ratio_set2", "pvalue_set2", "common_genes_set2")
    
    
    
    # Prepare for differentional analysis
    set1join$gene_overlap_set1[is.na(set1join$gene_overlap_set1)] <- 0
    set1join$gene_overlap_set2[is.na(set1join$gene_overlap_set2)] <- 0
    set1join[is.na(set1join)] <- 1
    
    # Obtain differentional phenotypes from geneset
    if (adjust == "bonferroni") {
      tabledif1u <- dplyr::filter(set1join, set1join$bonferroni_pvalue_set2 > p.value)
      # The same for genesetcomaper
      set2join <- dplyr::filter(prunionn, prunionn$bonferroni_pvalue.y <= p.value)
    }
    
    if (adjust == "fdr") {
      tabledif1u <- dplyr::filter(set1join, set1join$FDR_adjust_set2 > p.value)
      # The same for genesetcomaper
      set2join <- dplyr::filter(prunionn, prunionn$FDR_adjust.y <= p.value)
    }
    
    
    set2join <-  set2join[,c(1:5,7:12,14:17)]
    names(set2join) <- c("term_id", "term_name","source", "bonferroni_pvalue_set1", "FDR_adjust_set1", "gene_overlap_set1", "overlap_ratio_set1", "pvalue_set1", "common_genes_set1","bonferroni_pvalue_set2","FDR_adjust_set2", "gene_overlap_set2", "overlap_ratio_set2", "pvalue_set2", "common_genes_set2")
    
    set2join$gene_overlap_set1[is.na(set2join$gene_overlap_set1)] <- 0
    set2join$gene_overlap_set2[is.na(set2join$gene_overlap_set2)] <- 0
    set2join[is.na(set2join)] <- 1
    
    if (adjust == "bonferroni") {
      tabledif2u <- dplyr::filter(set2join, set2join$bonferroni_pvalue_set1 > p.value)
    }
    
    if (adjust == "fdr") {
      tabledif2u <- dplyr::filter(set2join, set2join$FDR_adjust_set1 > p.value)
    }
    # Order differentional phenotypes tables
    tabledif2u <- tabledif2u[order(tabledif2u$bonferroni_pvalue_set2),]
    tabledif1u <- tabledif1u[order(tabledif1u$bonferroni_pvalue_set1),]
    
    # Get fusion df term
    setfusion <- rbind(set1join,set2join)
    setfusion <- unique(setfusion)
    
    
    # get and comapre relevant phenotypes terms
    relevantesunion <- full_join(relevantes,targetrelevant, by = "term_id")
    runionn <- full_join(relevantes,targetrelevant, by = c("term_id","term_name", "source"))
    nrelevantesunion <- as.numeric(nrow(relevantesunion))
    relevantpresent <- inner_join(relevantes,targetrelevant, by = "term_id")
    ncompartidosr <- nrow(relevantpresent)
    
    # Obtain POR
    if (POR=="jaccard"){ 
      scoreenrich <- ncompartidosr/nrelevantesunion
    }
    
    # FALTA TOCAR ESTE
    if (POR=="forbes"){ 
      scoreenrich <- forbes(relevantes$term_id, targetrelevant$term_id)
    }
    
    # Inner with all phenotypes
    targetenrichttr <- targetenrich[,c(1:5,7:10)]
    names(targetenrichttr) <- c("term_id", "term_name", "source", "bonferroni_pvalue_set2", "FDR_adjust_set2","gene_overlap_set2", "overlap_ratio_set2", "pvalue_set2", "common_genes_set2")
    
    # Table with phenotypes
    triucecompare <- merge.data.frame(enrichresult,targetenrichttr, by = c("term_id","term_name", "source"))
    
    if (adjust == "bonferroni") {
      # Label if the phenotype is relevant in both datasets
      triucecompare <- triucecompare %>%
        mutate(shared_enrichment = ifelse(bonferroni_pvalue_set2 <= p.value & bonferroni_pvalue <= p.value, "Yes", "No"))
      
      # All shared phenotypes
      triucecompare <- arrange(triucecompare, desc(shared_enrichment), bonferroni_pvalue_set2)
    }
    
    
    if (adjust == "fdr") {
      # Label if the phenotype is relevant in both datasets
      triucecompare <- triucecompare %>%
        mutate(shared_enrichment = ifelse(FDR_adjust_set2 <= p.value & FDR_adjust <= p.value, "Yes", "No"))
      
      # All shared phenotypes
      triucecompare <- arrange(triucecompare, desc(shared_enrichment), FDR_adjust_set2)
    }
    
    
    
    targetrelevanttr <- targetrelevant[,c(1:5,7:10)]
    names(targetrelevanttr) <- c("term_id", "term_name", "source","bonferroni_pvalue_set2","FDR_adjust_set2","gene_overlap_set2", "overlap_ratio_set2", "pvalue_set2", "common_genes_set2")
    
    triucefinal <- merge.data.frame(relevantes,targetrelevanttr, by = c("term_id","term_name", "source"))
    
    triucefinal <- triucefinal[order(triucefinal$bonferroni_pvalue_set),]
    
    nspor <- nrow(triucefinal)
    
    # Regression analysis
    regresiontarget <- stats::lm(gene_overlap_set1  ~ gene_overlap_set2, data = setfusion, na.action=na.exclude)
    ortarget <- summary(regresiontarget)
    
    
    # Preparation for simulation analysis
    # We need all genes
    genAll <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/geneproteinlist.txt"))
    genAll <- as.data.frame(genAll)
    names(genAll) <- c("symbol","entrez")
    
    
    # We transform into a vector of entrez numbers
    genAll <- as.vector(genAll$entrez)
    
    # set seed
    set.seed(seed)
    
    # To split simulation analysis
    nulltestnumbernue =  nulltestnumber/2
    
    #We will keep a number of  a random set of genes with the same number of genes as the genes defined
    genesRandom <- replicate(nulltestnumbernue, sample(genAll, length(genesetcompare), replace=TRUE))
    genesRandomg2 <- replicate(nulltestnumbernue, sample(genAll, length(geneset), replace=TRUE))
    
    # Detect cores and work with a half of them
    nworkers <- (detectCores()/2)
    
    # Parallel simulation analysis
    if (linux==T){
      cl <- makeForkCluster(nworkers)
      clusterEvalQ(cl, library(PhenoExam))
      enrichrandomlist <- parApply(cl, genesRandom, 2, function(x) PhenoEnrichRandomOpt(x,database,organism="human"))
      enrichrandomlistg2 <- parApply(cl, genesRandomg2, 2, function(x) PhenoEnrichRandomOpt(x,database,organism="human"))
      stopCluster(cl)
    }else{  enrichrandomlist <- apply(genesRandom, 2, function(x) PhenoEnrichRandomOpt(x,database,organism="human"))
    enrichrandomlistg2 <- apply(genesRandomg2, 2, function(x) PhenoEnrichRandomOpt(x,database,organism="human"))}
    
    if (adjust == "bonferroni") {
      randomenrichrelevant <- lapply(enrichrandomlist, dplyr::filter, bonferroni_pvalue <= p.value)
      randomenrichrelevantg2 <- lapply(enrichrandomlistg2, dplyr::filter, bonferroni_pvalue <= p.value)
    }
    
    if (adjust == "fdr") {
      randomenrichrelevant <- lapply(enrichrandomlist, dplyr::filter, FDR_adjust <= p.value)
      randomenrichrelevantg2 <- lapply(enrichrandomlistg2, dplyr::filter, FDR_adjust <= p.value)
    }
    
    # Intersección RPOR
    phenopresent <- lapply(enrichrandomlist, inner_join, x = enrichresult, by = "term_name")
    ncompartidos <- lapply(phenopresent, nrow)
    
    phenopresentg2 <- lapply(enrichrandomlistg2, inner_join, x = targetenrich, by = "term_name")
    ncompartidosg2 <- lapply(phenopresentg2, nrow)
    
    
    # Unión RPOR
    phenopresentunion <- lapply(enrichrandomlist, full_join, x = enrichresult, by = "term_name")
    ncompartidosunion <- lapply(phenopresentunion, nrow)
    
    phenopresentuniong2 <- lapply(enrichrandomlistg2, full_join, x = targetenrich, by = "term_name")
    ncompartidosuniong2 <- lapply(phenopresentuniong2, nrow)
    
    SharedScore <- mapply(function(x,y) (x/y), x = ncompartidos, y = ncompartidosunion, SIMPLIFY = F)
    SharedScoreg2 <- mapply(function(x,y) (x/y), x = ncompartidosg2, y = ncompartidosuniong2, SIMPLIFY = F)
    
    
    
    # Aqui tengo que adaptarlo al forbes
    
    if (POR=="jaccard"){ 
      SharedScore <- c(SharedScore, SharedScoreg2)
    }
    
    if (POR=="forbes"){ 
      relev1 <- lapply(enrichrandomlist, function(x,y) {
        forbes(x$term_id,y)
      }, y = enrichresult$term_id)
      
      relev2 <- 
        lapply(enrichrandomlistg2, function(x,y) {
          forbes(x$term_id,y)
        }, y = targetenrich$term_id)
      
      
      SharedScore <- c(relev1, relev2)
    }
    
    
    # Intersección POR
    colnames <- c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust","gene_associated_in_db",  "gene_overlap", "overlap_ratio", "pvalue")
    randomenrichrelevant <- lapply(randomenrichrelevant, setNames, colnames)
    relevantpresent <- lapply(randomenrichrelevant, inner_join, x = relevantes, by = "term_name")
    ncompartidosr <- lapply(relevantpresent, nrow)
    
    randomenrichrelevantg2 <- lapply(randomenrichrelevantg2, setNames, colnames)
    relevantpresentg2 <- lapply(randomenrichrelevantg2, inner_join, x = targetrelevant, by = "term_name")
    ncompartidosrg2 <- lapply(relevantpresentg2, nrow)
    
    # Unión POR
    relevantpresentunion <- lapply(randomenrichrelevant, full_join, x = relevantes, by = "term_name")
    ncompartidosrunion <- lapply(relevantpresentunion, nrow)
    
    relevantpresentuniong2 <- lapply(randomenrichrelevantg2, full_join, x = targetrelevant, by = "term_name")
    ncompartidosruniong2 <- lapply(relevantpresentuniong2, nrow)
    
    EnrichScore <- mapply(function(x,y) (x/y), x = ncompartidosr, y = ncompartidosrunion, SIMPLIFY = F)
    EnrichScoreg2 <- mapply(function(x,y) (x/y), x = ncompartidosrg2, y = ncompartidosruniong2, SIMPLIFY = F)
    
    if (POR=="jaccard"){ 
      EnrichScore <- c(EnrichScore, EnrichScoreg2)
    }
    
    if (POR=="forbes"){ 
      porrelev1 <- lapply(randomenrichrelevant, function(x,y) {
        forbes(x$term_id,y)
      }, y = relevantes$term_id)
      
      porrelev2 <- 
        lapply(randomenrichrelevantg2, function(x,y) {
          forbes(x$term_id,y)
        }, y = targetrelevant$term_id)
      
      
      EnrichScore <- c(porrelev1, porrelev2)
    }
    
    
    
  }
  
  
  if (conditional_case == T) {
    # Obtener el análisis de enriquecimiento del dataset de referencia
    referencia = PhenoEnrichGenes(genes=geneset, database, organism, url = url)
    enrichresult = as.data.frame(referencia$htmltabla)
    
    
    # Obtain the significative phenotypes terms from geneset
    if (adjust == "bonferroni") {
      relevantes <- dplyr::filter(enrichresult, enrichresult$bonferroni_pvalue <= p.value)
    }
    
    if (adjust == "fdr") {
      relevantes <- dplyr::filter(enrichresult, enrichresult$FDR_adjust <= p.value)
    }
    
    
    
    
    # Todos los genes de referencia en los que no exista overlap con el dataset comparado deben salir
    # de la ecuación
    targetenrich = PhenoEnrichCompare(geneset, genesetcompare, database, organism, url = url)
    targetenrich = as.data.frame(targetenrich$htmltabla)
    nphenotargetall <- nrow(targetenrich)
    # names(targetenrich) <- c("term_id", "term_name", "bonferroni_pvalue","gene_associated_in_db", "total_size_ratio", "gene_overlap", "input_size_ratio", "input_total_fold", "overlap_ratio", "pvalue", "common_genes_set2")
    
    
    # Tomar los fenotipos estadisticamente significativos
    if (adjust == "bonferroni") {
      targetrelevant <- dplyr::filter(targetenrich, targetenrich$bonferroni_pvalue <= p.value)
    }
    
    if (adjust == "fdr") {
      targetrelevant <- dplyr::filter(targetenrich, targetenrich$FDR_adjust <= p.value)
    }
    
    nphenotarget <- nrow(targetrelevant)
    
    
    # Pregunta 1 overlap de términos aunque estos no sean relevantes
    phenounion <- full_join(enrichresult, targetenrich, by = "term_id")
    prunionn <- full_join(enrichresult,targetenrich, by = c("term_id","term_name", "source"))
    nphenounion <- as.numeric(nrow(phenounion))
    #nphenounion <- unique(phenounion$term_id)
    #nphenounion <- length(nphenounion) # Número de fenotipos totales UNICOS DE LA UNION
    phenopresentall <- inner_join(enrichresult, targetenrich, by = "term_id") # Innner de fenotipos totales entre los dos datase
    ncompartidosall <- nrow(phenopresentall) # Número de fenotipos totales compartidos
    
    if (POR=="jaccard"){ 
      scoreshared <- ncompartidosall/nphenounion
    }
    
    if (POR=="forbes"){ 
      scoreshared <- forbes(enrichresult$term_id, targetenrich$term_id)
    }
    
    
    
    if (adjust == "bonferroni") {
      set1join <- dplyr::filter(prunionn, prunionn$bonferroni_pvalue.x <= p.value)
    }
    
    if (adjust == "fdr") {  
      set1join <- dplyr::filter(prunionn, prunionn$FDR_adjust.x <= p.value)
    }
    
    set1join <-  set1join[,c(1:5,7:12,14:17)]
    names(set1join) <- c("term_id", "term_name","source", "bonferroni_pvalue_set1", "FDR_adjust_set1", "gene_overlap_set1", "overlap_ratio_set1", "pvalue_set1", "common_genes_set1","bonferroni_pvalue_set2","FDR_adjust_set2", "gene_overlap_set2", "overlap_ratio_set2", "pvalue_set2", "common_genes_set2")
    
    
    
    set1join$gene_overlap_set1[is.na(set1join$gene_overlap_set1)] <- 0
    set1join$gene_overlap_set2[is.na(set1join$gene_overlap_set2)] <- 0
    set1join[is.na(set1join)] <- 1
    
    if (adjust == "bonferroni") {
      tabledif1u <- dplyr::filter(set1join, set1join$bonferroni_pvalue_set2 > p.value)
      # The same for genesetcomaper
      set2join <- dplyr::filter(prunionn, prunionn$bonferroni_pvalue.y <= p.value)
    }
    
    if (adjust == "fdr") {
      tabledif1u <- dplyr::filter(set1join, set1join$FDR_adjust_set2 > p.value)
      # The same for genesetcomaper
      set2join <- dplyr::filter(prunionn, prunionn$FDR_adjust.y <= p.value)
    }
    
    
    set2join <-  set2join[,c(1:5,7:12,14:17)]
    names(set2join) <- c("term_id", "term_name","source", "bonferroni_pvalue_set1", "FDR_adjust_set1", "gene_overlap_set1", "overlap_ratio_set1", "pvalue_set1", "common_genes_set1","bonferroni_pvalue_set2","FDR_adjust_set2", "gene_overlap_set2", "overlap_ratio_set2", "pvalue_set2", "common_genes_set2")
    
    
    
    
    
    set2join$gene_overlap_set1[is.na(set2join$gene_overlap_set1)] <- 0
    set2join$gene_overlap_set2[is.na(set2join$gene_overlap_set2)] <- 0
    set2join[is.na(set2join)] <- 1
    
    
    
    
    if (adjust == "bonferroni") {
      tabledif2u <- dplyr::filter(set2join, set2join$bonferroni_pvalue_set1 > p.value)
    }
    
    if (adjust == "fdr") {
      tabledif2u <- dplyr::filter(set2join, set2join$FDR_adjust_set1 > p.value)
    }
    
    
    
    
    
    tabledif2u <- tabledif2u[order(tabledif2u$bonferroni_pvalue_set2),]
    tabledif1u <- tabledif1u[order(tabledif1u$bonferroni_pvalue_set1),]
    
    # Unir estos dos datasets
    setfusion <- rbind(set1join,set2join)
    setfusion <- unique(setfusion)
    
    
    # Comparación de fenotipos enriquecidos en ambos dataset, número de compartidos
    relevantesunion <- full_join(relevantes,targetrelevant, by = "term_id")
    nrelevantesunion <- as.numeric(nrow(relevantesunion))
    relevantpresent <- inner_join(relevantes,targetrelevant, by = "term_id")
    ncompartidosr <- nrow(relevantpresent)
    
    
    
    if (POR=="jaccard"){ 
      scoreenrich <- ncompartidosr/nrelevantesunion
    }
    
    if (POR=="forbes"){ 
      scoreenrich <- forbes(relevantes$term_id, targetrelevant$term_id)
    }
    
    
    
    # Inner with all phenotypes
    targetenrichttr <- targetenrich[,c(1:5,7:10)]
    names(targetenrichttr) <- c("term_id", "term_name", "source", "bonferroni_pvalue_set2", "FDR_adjust_set2","gene_overlap_set2", "overlap_ratio_set2", "pvalue_set2", "common_genes_set2")
    
    # Table with phenotypes
    triucecompare <- merge.data.frame(enrichresult,targetenrichttr, by = c("term_id","term_name", "source"))
    
    if (adjust == "bonferroni") {
      # Label if the phenotype is relevant in both datasets
      triucecompare <- triucecompare %>%
        mutate(shared_enrichment = ifelse(bonferroni_pvalue_set2 <= p.value & bonferroni_pvalue <= p.value, "Yes", "No"))
      
      # All shared phenotypes
      triucecompare <- arrange(triucecompare, desc(shared_enrichment), bonferroni_pvalue_set2)
    }
    
    
    if (adjust == "fdr") {
      # Label if the phenotype is relevant in both datasets
      triucecompare <- triucecompare %>%
        mutate(shared_enrichment = ifelse(FDR_adjust_set2 <= p.value & FDR_adjust <= p.value, "Yes", "No"))
      
      # All shared phenotypes
      triucecompare <- arrange(triucecompare, desc(shared_enrichment), FDR_adjust_set2)
    }
    
    
    
    targetrelevanttr <- targetrelevant[,c(1:5,7:10)]
    names(targetrelevanttr) <- c("term_id", "term_name", "source","bonferroni_pvalue_set2","FDR_adjust_set2","gene_overlap_set2", "overlap_ratio_set2", "pvalue_set2", "common_genes_set2")
    
    triucefinal <- merge.data.frame(relevantes,targetrelevanttr, by = c("term_id","term_name", "source"))
    
    triucefinal <- triucefinal[order(triucefinal$bonferroni_pvalue_set),]
    nspor <- nrow(triucefinal)
    
    # Esto nos da que cuanto más enriquecido en un sitio más enriquecido en el otro sitio
    regresiontarget <- stats::lm(gene_overlap_set1  ~ gene_overlap_set2, data = setfusion, na.action=na.exclude)
    ortarget <- summary(regresiontarget)
    
    
    
    
    
    if (organism == "human") {
      nonoverlapgenes <- foundgeneoverlap(geneset, genesetcompare)
    }
    if (organism == "mouse") {
      targetgenesm <- as.data.frame(genesetcompare)
      names(targetgenesm) <- "mouse_symbol"
      targetgenesm$mouse_symbol <- as.character(targetgenesm$mouse_symbol)
      
      # Read database
      homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                "/homology.csv"), header = T)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(targetgenesm, homologia, by = "mouse_symbol")
      onlygenes <- onlygenes$human_symbol
      onlygenes <- as.data.frame(onlygenes)
      names(onlygenes) <- "symbol"
      onlygenes$symbol <- as.character(onlygenes$symbol)
      targetgenesm = onlygenes
      
      referencegenesm <- as.data.frame(geneset)
      names(referencegenesm) <- "mouse_symbol"
      referencegenesm$mouse_symbol <- as.character(referencegenesm$mouse_symbol)
      
      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(referencegenesm, homologia, by = "mouse_symbol")
      onlygenes <- onlygenes$human_symbol
      onlygenes <- as.data.frame(onlygenes)
      names(onlygenes) <- "symbol"
      onlygenes$symbol <- as.character(onlygenes$symbol)
      referencegenesm = onlygenes
      nonoverlapgenes <- foundgeneoverlap(referencegenesm, targetgenesm)
    }
    
    genAll <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/geneproteinlist.txt"))
    genAll <- as.data.frame(genAll)
    names(genAll) <- c("symbol","entrez")
    
    
    # We will keep only the selected genes
    genAll <- dplyr::anti_join(genAll, nonoverlapgenes, by = "entrez")
    genAll <- as.vector(genAll$entrez)
    
    set.seed(seed)
    
    #We will keep a number of  a random set of genes with the same number of genes as the genes defined.
    nulltestnumbernue =  nulltestnumber/2
    genesRandom <- replicate(nulltestnumbernue, sample(genAll, length(genesetcompare), replace=TRUE))
    genesRandomg2 <- replicate(nulltestnumbernue, sample(genAll, length(geneset), replace=TRUE))
    
    # Detecto los corazones y utilizo la mitad
    nworkers <- (detectCores()/2)
    
    
    # ME HE QUEDADO CON PHENOENRICHCOMPARERANDOM
    # Paralelizo la opreación de calcular en enriquecimiento para cada set de genes aleatorios
    
    
    
    if (linux==T){
      cl <- makeForkCluster(nworkers)
      clusterEvalQ(cl, library(PhenoExam))
      enrichrandomlist <- parApply(cl, genesRandom, 2, function(x) PhenoEnrichCompareRandom(nonoverlapgenes, geneset, x, database,organism="human"))
      enrichrandomlistg2 <- parApply(cl, genesRandomg2, 2, function(x) PhenoEnrichCompareRandom(nonoverlapgenes, geneset, x,database,organism="human"))
      stopCluster(cl)
    }else{  enrichrandomlist <- apply(genesRandom, 2, function(x) PhenoEnrichCompareRandom(nonoverlapgenes, geneset, x,database,organism="human"))
    enrichrandomlistg2 <- apply(genesRandomg2, 2, function(x) PhenoEnrichCompareRandom(nonoverlapgenes, geneset,x,database,organism="human"))}
    
    
    
    if (adjust == "bonferroni") {
      randomenrichrelevant <- lapply(enrichrandomlist, dplyr::filter, bonferroni_pvalue <= p.value)
      randomenrichrelevantg2 <- lapply(enrichrandomlistg2, dplyr::filter, bonferroni_pvalue <= p.value)
    }
    
    if (adjust == "fdr") {
      randomenrichrelevant <- lapply(enrichrandomlist, dplyr::filter, FDR_adjust <= p.value)
      randomenrichrelevantg2 <- lapply(enrichrandomlistg2, dplyr::filter, FDR_adjust <= p.value)
    }
    
    
    
    # Intersección RPOR
    phenopresent <- lapply(enrichrandomlist, inner_join, x = enrichresult, by = "term_name")
    ncompartidos <- lapply(phenopresent, nrow)
    
    phenopresentg2 <- lapply(enrichrandomlistg2, inner_join, x = targetenrich, by = "term_name")
    ncompartidosg2 <- lapply(phenopresentg2, nrow)
    
    
    # Unión RPOR
    phenopresentunion <- lapply(enrichrandomlist, full_join, x = enrichresult, by = "term_name")
    ncompartidosunion <- lapply(phenopresentunion, nrow)
    
    phenopresentuniong2 <- lapply(enrichrandomlistg2, full_join, x = targetenrich, by = "term_name")
    ncompartidosuniong2 <- lapply(phenopresentuniong2, nrow)
    
    SharedScore <- mapply(function(x,y) (x/y), x = ncompartidos, y = ncompartidosunion, SIMPLIFY = F)
    SharedScoreg2 <- mapply(function(x,y) (x/y), x = ncompartidosg2, y = ncompartidosuniong2, SIMPLIFY = F)
    
    if (POR=="jaccard"){ 
      SharedScore <- c(SharedScore, SharedScoreg2)
    }
    
    if (POR=="forbes"){ 
      relev1 <- lapply(enrichrandomlist, function(x,y) {
        forbes(x$term_id,y)
      }, y = enrichresult$term_id)
      
      relev2 <- 
        lapply(enrichrandomlistg2, function(x,y) {
          forbes(x$term_id,y)
        }, y = targetenrich$term_id)
      
      
      SharedScore <- c(relev1, relev2)
    }
    
    
    # Intersección POR
    colnames <- c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust","gene_associated_in_db",  "gene_overlap", "overlap_ratio", "pvalue")
    randomenrichrelevant <- lapply(randomenrichrelevant, setNames, colnames)
    relevantpresent <- lapply(randomenrichrelevant, inner_join, x = relevantes, by = "term_name")
    ncompartidosr <- lapply(relevantpresent, nrow)
    
    randomenrichrelevantg2 <- lapply(randomenrichrelevantg2, setNames, colnames)
    relevantpresentg2 <- lapply(randomenrichrelevantg2, inner_join, x = targetrelevant, by = "term_name")
    ncompartidosrg2 <- lapply(relevantpresentg2, nrow)
    
    # Unión POR
    relevantpresentunion <- lapply(randomenrichrelevant, full_join, x = relevantes, by = "term_name")
    ncompartidosrunion <- lapply(relevantpresentunion, nrow)
    
    relevantpresentuniong2 <- lapply(randomenrichrelevantg2, full_join, x = targetrelevant, by = "term_name")
    ncompartidosruniong2 <- lapply(relevantpresentuniong2, nrow)
    
    EnrichScore <- mapply(function(x,y) (x/y), x = ncompartidosr, y = ncompartidosrunion, SIMPLIFY = F)
    EnrichScoreg2 <- mapply(function(x,y) (x/y), x = ncompartidosrg2, y = ncompartidosruniong2, SIMPLIFY = F)
    
    if (POR=="jaccard"){ 
      EnrichScore <- c(EnrichScore, EnrichScoreg2)
    }
    
    if (POR=="forbes"){ 
      porrelev1 <- lapply(randomenrichrelevant, function(x,y) {
        forbes(x$term_id,y)
      }, y = relevantes$term_id)
      
      porrelev2 <- 
        lapply(randomenrichrelevantg2, function(x,y) {
          forbes(x$term_id,y)
        }, y = targetrelevant$term_id)
      
      
      EnrichScore <- c(porrelev1, porrelev2)
    }
    
    
    
    
  }
  
  # We will keep with the phenoratios from random genes greater than the phenoratios from predictions genes
  # EnrichScore
  EnrichScoreMayores = which(EnrichScore > scoreenrich)
  EnrichScoreMayoresL = as.numeric(length(EnrichScoreMayores))
  
  # SharedScore
  SharedScoreMayores = which(SharedScore > scoreshared)
  SharedScoreMayoresL = as.numeric(length(SharedScoreMayores))
  
  # We are going to calculate the pvalue
  pvalorEnrichScore = (1 + EnrichScoreMayoresL) / (1 + as.numeric(length(EnrichScore)))
  pvalorSharedScoreMayores = (1 + SharedScoreMayoresL) / (1 + as.numeric(length(SharedScore)))
  
  
  # Plot  Relaxed Phenotypic Overlap Ratio (RPOR)
  ratiosRandomHP2 <- unlist(SharedScore, use.names = FALSE)
  pr <- as.data.frame(ratiosRandomHP2)
  rpre <- as.data.frame(scoreshared)
  pl1 <- ggplot2::ggplot(pr, ggplot2::aes(x= ratiosRandomHP2)) + ggplot2::geom_histogram(binwidth= mean(ratiosRandomHP2)/30, position="identity", alpha=0.7, fill="red") + ggplot2::geom_vline( ggplot2::aes(xintercept= scoreshared),
                                                                                                                                                                                                color="blue", linetype="dashed", size=0.7)+
    labs(title="Relaxed Phenotypic Overlap Ratio (RPOR) from random genes vs geneset2",x= "Relaxed Phenotypic Overlap Ratio (RPOR)", y = "Frecuency in random gene sets")
  
  a <- ggplot2::ggplot_build(pl1)$layout$panel_scales_y[[1]]$range$range
  
  plotshared <- pl1 +  ggplot2::geom_text( ggplot2::aes(scoreshared-(scoreshared/25), mean(a), label = "Relaxed Phenotypic Overlap Ratio (POR) from", vjust=0,angle=90, family="Times")) +  ggplot2::geom_text( ggplot2::aes(scoreshared-(scoreshared/55), mean(a) ,label = "geneset2", vjust=0,angle=90, family="Times"))
  
  
  # Plot Phenotype Overlap Ratio (POR)
  ratiosRandomHP2 <- unlist(EnrichScore, use.names = FALSE)
  
  pr <- as.data.frame(ratiosRandomHP2)
  rpre <- as.data.frame(scoreenrich)
  pl1 <- ggplot2::ggplot(pr, ggplot2::aes(x= ratiosRandomHP2)) + ggplot2::geom_histogram(binwidth= 0.005, position="identity", alpha=0.7, fill="red") + ggplot2::geom_vline( ggplot2::aes(xintercept= scoreenrich),
                                                                                                                                                                             color="blue", linetype="dashed", size=0.7)+
    labs(title="Phenotype Overlap Ratio (POR)
 from random genes vs geneset2",x= "Phenotype Overlap Ratio (POR)
", y = "Frecuency in random gene sets")
  
  b <- ggplot2::ggplot_build(pl1)$layout$panel_scales_y[[1]]$range$range
  
  plotenrichs <- pl1 +  ggplot2::geom_text( ggplot2::aes(scoreenrich-(scoreenrich/40), mean(b) ,label = "POR from geneset2",vjust=0,angle=90, family="Times"))
  
  
  
  
  # Interactive plot
  if (adjust == "bonferroni") {
    pdif <- setfusion %>%  dplyr::mutate(Enrich = case_when( setfusion$bonferroni_pvalue_set1 <= 0.05 & setfusion$bonferroni_pvalue_set2 <= 0.05  ~ "SharedEnrich",
                                                             setfusion$bonferroni_pvalue_set1 > 0.05  ~ "Relevant in geneset2",
                                                             setfusion$bonferroni_pvalue_set2 > 0.05  ~ "Relevant in geneset1" )) %>%
      
      # prepare text for tooltip
      dplyr::mutate(text = paste("term_id: ", term_id, "\nterm_name: ", term_name, "\npvalue_set1: ", bonferroni_pvalue_set1, "\npvalue_set2: ", bonferroni_pvalue_set2, "\nGeneassociation_geneset1: ", gene_overlap_set1, "\nGeneassociation_geneset2: ", gene_overlap_set2,  sep="")) %>%
      
      # Classic ggplot
      ggplot2::ggplot( ggplot2::aes(x=gene_overlap_set2, y=gene_overlap_set1,  color = Enrich, text=text)) +
      ggplot2::theme(text = element_text(size=8)) +
      ggplot2::geom_point(alpha = 0.5) +
      #viridis::scale_color_viridis(discrete=TRUE, guide=FALSE) +
      #ggplot2::theme(legend.position="none") +
      labs(title = "Phenotype relevance association analysis for gene sets", x = "Number of genes associated in gene set2", y = "Number of genes associated in gene set1") +
      theme(plot.title = element_text(size = 15),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
  }
  if (adjust == "fdr") {
    pdif <- setfusion %>%  dplyr::mutate(Enrich = case_when( setfusion$FDR_adjust_set1 <= 0.05 & setfusion$FDR_adjust_set2 <= 0.05  ~ "SharedEnrich",
                                                             setfusion$FDR_adjust_set1 > 0.05  ~ "Relevant in geneset2",
                                                             setfusion$FDR_adjust_set2 > 0.05  ~ "Relevant in geneset1" )) %>%
      
      # prepare text for tooltip
      dplyr::mutate(text = paste("term_id: ", term_id, "\nterm_name: ", term_name, "\npvalue_set1: ", FDR_adjust_set1, "\npvalue_set2: ", FDR_adjust_set2, "\nGeneassociation_geneset1: ", gene_overlap_set1, "\nGeneassociation_geneset2: ", gene_overlap_set2,  sep="")) %>%
      
      # Classic ggplot
      ggplot2::ggplot( ggplot2::aes(x=gene_overlap_set2, y=gene_overlap_set1,  color = Enrich, text=text)) +
      ggplot2::theme(text = element_text(size=8)) +
      ggplot2::geom_point(alpha = 0.5) +
      #viridis::scale_color_viridis(discrete=TRUE, guide=FALSE) +
      #ggplot2::theme(legend.position="none") +
      labs(title = "Phenotype relevance association analysis for gene sets", x = "Number of genes associated in gene set2", y = "Number of genes associated in gene set1") +
      theme(plot.title = element_text(size = 15),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
  }
  
  
  # turn ggplot interactive with plotly
  ppdif <- plotly::ggplotly(pdif, tooltip="text", parse = T) %>%
    add_lines(x = ~gene_overlap_set2, y = fitted(regresiontarget), name = "Regression line")
  
  
  
  
  
  if (pvalorEnrichScore <= 0.05) {
    mensaje <-  paste0("Gene set2 and gene set1 shared ",ncompartidosall, " phenotypic terms (out of ",nphenounion, " unique phenotypic terms in both), that yields a Relaxed Phenotypic Overlap Ratio (RPOR) of ",round(scoreshared,3), " (p < ",round(pvalorSharedScoreMayores,6),"). They shared ",nspor, " significant phenotypic terms (out of ",nrelevantesunion, " unique significant phenotypic terms in both), that yields a Phenotypic Overlap Ratio (POR) of ",round(scoreenrich,3), " (p < ",round(pvalorEnrichScore,6),"). Phenotype relevance association analysis for gene sets (i.e. whether the shared phenotypes are similar in relevance, i.e. in the number of genes associated with them, within each gene set) results in an adjusted R squared of ",round(ortarget$adj.r.squared,3)," ( p < ",ortarget$coefficients[2,4],") which suggests that an important portion of the common phenotypes are similar in relevance.  The p-values were obtained through randomization of ",nulltestnumber, " random gene sets.  Gene set 1 and gene set 2 are statistically significantly similar to each other in phenotypic terms.")
  }
  if (pvalorEnrichScore > 0.05) {
    mensaje <-  paste0("Gene set2 and gene set1 shared ",ncompartidosall, " phenotypic terms (out of ",nphenounion, " unique phenotypic terms in both), that yields a Relaxed Phenotypic Overlap Ratio (RPOR) of ",round(scoreshared,3), " (p < ",round(pvalorSharedScoreMayores,6),"). They shared ",nspor, " significant phenotypic terms (out of ",nrelevantesunion, " unique significant phenotypic terms in both), that yields a Phenotypic Overlap Ratio (POR) of ",round(scoreenrich,3), " (p < ",round(pvalorEnrichScore,6),"). Phenotype relevance association analysis for gene sets (i.e. whether the shared phenotypes are similar in relevance, i.e. in the number of genes associated with them, within each gene set) results in an adjusted R squared of ",round(ortarget$adj.r.squared,3)," ( p < ",ortarget$coefficients[2,4],") .  The p-values were obtained through randomization of ",nulltestnumber, " random gene sets.  Gene set 1 and gene set 2 are not statistically significantly similar to each other in phenotypic terms.")
  }
  
  newList <- list("phenomessage" = mensaje, "data.enrich" = triucefinal, "datatotal" = triucecompare, "RPOR" = plotshared, "POR" = plotenrichs,  "randomvalues" = pr,  "lm" = ortarget, "r.squared" = round(ortarget$adj.r.squared,5), "lmp.value" = ortarget$coefficients[2,4],  "plotdif" = ppdif, "tabledif1" = tabledif1u, "tabledif2" = tabledif2u)
  
  return(newList)
}

#' @title PhenoEnrichGenes
#' @description This function do a phenotype enrichment analysis and obtain a table and a graph
#' @param genes A vector with human gene symbol.
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @param plotn An integer with the number of term that you want in the plot. The default value is 30.
#' @param url A logical parameter. T if you want interactive links
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export



PhenoEnrichGenes = function(genes, database = "HPO", organism = "human", plotn=30, url = F){
  
  if (organism == "human") {
    
    # Make two empty df with names
    rbjointri = data.frame(matrix(vector(), 0, 10,
                                  dimnames=list(c(), c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust",
                                                       "genes_associated_in_db","gene_overlap", "overlap_ratio"
                                                       ,"pvalue","common_genes"))),
                           stringsAsFactors=F)
    rbjoinhtml <-rbjointri
    
    # For each database name in a vector run its loadder
    for (db in database) {
      
      # Get the loadder name
      getnamedb <- paste0("",db,"loadder")
      
      # Run the loadder
      peaload <- get(getnamedb)(genes,organism = "human", plotnumber = plotn, url = url)
      
      # join data frames generated
      rbjoinhtml <- rbind(rbjoinhtml, peaload$htmltabla)
      rbjointri <- rbind(rbjointri, peaload$alldata)
      
      # Order by pvalue
      rbjoinhtml = rbjoinhtml[order(rbjoinhtml$bonferroni_pvalue),] # htmltabla
      rbjointri = rbjointri[order(rbjointri$bonferroni_pvalue),] # alldata
      
      # To plot
      triucu = rbjointri[1:plotn,]
    }
  }
  
  if (organism == "mouse") {
    
    genes <- as.data.frame(genes)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)
    
    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genes = unique(onlygenes)
    genes = genes$symbol
    
    # Make two empty df with names
    rbjointri = data.frame(matrix(vector(), 0, 10,
                                  dimnames=list(c(), c("term_id", "term_name", "source", "bonferroni_pvalue","FDR_adjust",
                                                       "genes_associated_in_db","gene_overlap", "overlap_ratio"
                                                       ,"pvalue","common_genes"))),
                           stringsAsFactors=F)
    rbjoinhtml <-rbjointri
    
    # For each database name in a vector run its loadder
    for (db in database) {
      
      # Get the loadder name
      getnamedb <- paste0("",db,"loadder")
      
      # Run the loadder
      peaload <- get(getnamedb)(genes,organism = "human", plotnumber = plotn,  url = url)
      
      # join data frames generated
      rbjoinhtml <- rbind(rbjoinhtml, peaload$htmltabla)
      rbjointri <- rbind(rbjointri, peaload$alldata)
      
      # Order by pvalue
      rbjoinhtml = rbjoinhtml[order(rbjoinhtml$bonferroni_pvalue),] # htmltabla
      rbjointri = rbjointri[order(rbjointri$bonferroni_pvalue),] # alldata
      
      # To plot
      triucu = rbjointri[1:plotn,]
    }
  }
  
  # Esto es para que en la generación del gráfico no salga el no fenotipo
  #triucu <- filter(triucu, triucu$term_id != "CRB:XXX")
  #triucu <- filter(triucu, triucu$term_id != "wXXX") no pasa nada si
  
  # To plot
  p <- triucu %>%
    
    # prepare text for tooltip
    dplyr::mutate(text = paste("ID Number: ", term_id, "\nTerm: ", term_name,  "\nGene Ratio: ", overlap_ratio, "\nBonferroni pvalue: ", bonferroni_pvalue, "\nFDR_adjust: ", FDR_adjust, "\nCommon genes: ", common_genes)) %>%
    # Classic ggplot
    ggplot(  ggplot2::aes(x= gene_overlap , y= reorder(paste0("",term_id,": ",term_name,""), -bonferroni_pvalue),  fill = bonferroni_pvalue, text=text)) +
    theme (text = element_text(size= 8)) +
    geom_bar(stat="identity", alpha=0.7) +
    scale_fill_gradientn(colours = rainbow(2), trans= "log") +
    theme(legend.position="right") +
    labs(title = "Phenotype Enrichment Analysis", x = "Number of genes associated with the term", y = "", fill = "p.value") + #
    theme(plot.title = element_text(size = 15),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text = element_text(
            angle = 0))
  
  # turn ggplot interactive with plotly layout.title
  nuevvogg <- ggplotly(p, tooltip="text")  %>%
    config(mathjax = "cdn",
           toImageButtonOptions = list(
             format = "svg",
             filename = "myplot",
             width = 1000,
             height = 600
           )
    )
  
  newList <- list("alldata" = rbjointri, "htmltabla" = rbjoinhtml,   "graph" = nuevvogg)
  return(newList)
  
  
}




#' @title PhenoEnrichRandomOpt
#' @description This function makes a plot with enriched phenotype info.
#' @param genes A vector with human gene symbol.
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @param plotn An integer with the number of term that you want in the plot. The default value is 30.
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export



PhenoEnrichRandomOpt = function(genes,  database = "HPO", organism = "human"){
  
  
  if (organism == "human") {
    
    
    
    # Make one empty df with names
    rbjointri = data.frame(matrix(vector(), 0, 9,
                                  dimnames=list(c(), c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust",
                                                       "genes_associated_in_db","gene_overlap", "overlap_ratio"
                                                       ,"pvalue"))),
                           stringsAsFactors=F)
    
    
    # For each database name in a vector run its loadder
    for (db in database) {
      
      # Get the loadder name
      getnamedb <- paste0("",db,"loadderOpt")
      
      # Run the loadder
      peaload <- get(getnamedb)(genes,organism = "human", plotnumber = plotn)
      
      # join data frames generated
      rbjointri <- rbind(rbjointri, peaload)
      
      # Order by pvalue
      rbjointri = rbjointri[order(rbjointri$bonferroni_pvalue),] # alldata
      
    }
  }
  
  
  
  if (organism == "mouse") {
    genes <- as.data.frame(genes)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)
    
    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genes = unique(onlygenes)
    genes = genes$symbol
    
    
    # Make one empty df with names
    rbjointri = data.frame(matrix(vector(), 0, 9,
                                  dimnames=list(c(), c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust",
                                                       "genes_associated_in_db","gene_overlap", "overlap_ratio"
                                                       ,"pvalue"))),
                           stringsAsFactors=F)
    
    # For each database name in a vector run its loadder
    for (db in database) {
      
      # Get the loadder name
      getnamedb <- paste0("",db,"loadderOpt")
      
      # Run the loadder
      peaload <- get(getnamedb)(genes,organism = "human", plotnumber = plotn)
      
      # join data frames generated
      rbjointri <- rbind(rbjointri, peaload)
      
      # Order by pvalue
      rbjointri = rbjointri[order(rbjointri$bonferroni_pvalue),] # alldata
      
    }
    
    
  }
  
  newList <- rbjointri
  return(newList)
}




#' @title convertgenesymbolmessage
#' @description This function give a message with the homologues map
#' @param genes A vector with gene symbol.
#' @param organism A string parameter: "human" if you have human genes and "mouse if you have mouse genes
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export


convertgenesymbolmessage = function(genes, organism = "human"){

  if (organism == "human") {

    geneshomo <- as.data.frame(genes)
    names(geneshomo) <- "human_symbol"
    geneshomo$human_symbol <- as.character(geneshomo$human_symbol)
    geneshomo <- unique(geneshomo)
    ngeneshomo <- nrow(geneshomo)

    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)

    # We will keep only the selected genes
    onlygeneshomo <- dplyr::anti_join(geneshomo, homologia, by = "human_symbol")
    onlygeneshomo <- onlygeneshomo$human_symbol
    nonly <- length(onlygeneshomo)
    onlygeneshomo <- paste(onlygeneshomo,collapse=" ")


    porcenta <- round(((nonly/ngeneshomo)*100), digits = 2)


    if (porcenta == 0) {
      mensaje <-  paste0("All genes from the input have mice homologues genes")
    }
    if (porcenta != 0) {
      mensaje <-  paste0("PhenoExam could not find homologues of the human genes: ", onlygeneshomo," (",nonly,"/",ngeneshomo," genes) ",porcenta,"% of the input genes")

    }
    return(mensaje)
  }

  if (organism == "mouse") {
    genes <- as.data.frame(genes)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)
    genes <- unique(genes)
    ngenes <- nrow(genes)

    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)

    # We will keep only the selected genes
    onlygenes <- dplyr::anti_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$mouse_symbol
    nonly <- length(onlygenes)
    onlygenes <- paste(onlygenes,collapse=" ")


    porcenta <- round(((nonly/ngenes)*100), digits = 2)

    mensaje <-  paste0("PhenoExam could not find homologues of the mouse genes: ", onlygenes," (",nonly,"/",ngenes," genes) ",porcenta,"% of the input genes")

    return(mensaje)
  }
}


#' @title entrezmap
#' @description This function give a message with the entrez map
#' @param genes A vector with gene symbol.
#' @param human A logical parameter if you have human gene symbol
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export


entrezmap = function(genes, human = T){

  if (human) {
    mouse = F
    # Covert vector to a daframe
    genesentrez <- as.data.frame(genes)
    names(genesentrez) <- "Symbol"
    ngenes <- nrow(genesentrez)

    # Pasar estos genes a entrezid clusterprofiler
    genes.df <- bitr(genesentrez$Symbol, fromType = "SYMBOL",
                     toType = c("SYMBOL", "ENTREZID"),
                     OrgDb = "org.Hs.eg.db")
    genesentrez <- unique(genes.df)
    names(genesentrez) <- c("symbol","entrez")
    genesentrez$entrez <- as.numeric(as.character(genesentrez$entrez))

    # Read database
    custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                           "/geneproteinlist.txt"), header = T)
    names(custom) <- c("symbol","entrez")

    # We will keep only the selected genes
    onlygenesentrez <- dplyr::anti_join(genesentrez, custom, by = "entrez")

    # Ahora si tengo los dos y podría mostrar que genes no se mapean o entrar dentro
    symyentrez <- unique(onlygenesentrez[,1:2])

    genesnoentrez <- symyentrez$symbol
    nonly <- as.numeric(length(genesnoentrez))
    genesnoentrez <- paste(genesnoentrez,collapse=" ")


    porcenta <- round(((nonly/ngenes)*100), digits = 2)

    if (porcenta == 0) {
      mensaje <-  paste0("All genes from the input will be used")
    }
    if (porcenta != 0) {
      mensaje <-  paste0("PhenoExam could not use or find these genes: ", genesnoentrez," (",nonly,"/",ngenes," genes) ",porcenta,"% of the input genes")

    }

    return(mensaje)
  }
}


#' @title PhenoEnrichCompare
#' @description This function is used for RandomComparePheno.
#' @param genes A vector with human gene symbol.
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @param plotn An integer with the number of term that you want in the plot. The default value is 30.
#' @param url A logical parameter. T if you want interactive links
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export




# Función para obtener el enriquecimiento de términos HPO Y CRB
PhenoEnrichCompare = function(geneset, genesetcompare, database = "HPO", organism = "human", plotn=30, plotnumber = 1:plotn, url = F){
  
  
  
  
  if (organism == "human") {
    
    
    # Find gene overlap
    onlygenes <- foundgeneoverlap(geneset, genesetcompare)
    
    
    # Make two empty df with names
    rbjointri = data.frame(matrix(vector(), 0, 10,
                                  dimnames=list(c(), c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust",
                                                       "genes_associated_in_db","gene_overlap", "overlap_ratio"
                                                       ,"pvalue","common_genes"))),
                           stringsAsFactors=F)
    rbjoinhtml <-rbjointri
    
    # For each database name in a vector run its loadder
    for (db in database) {
      
      # Get the loadder name
      getnamedb <- paste0("",db,"loaddercompare")
      
      # Run the loadder
      peaload <- get(getnamedb)(onlygenes, genesetcompare,organism = "human", plotnumber = plotn, url = url)
      
      # join data frames generated
      rbjoinhtml <- rbind(rbjoinhtml, peaload$htmltabla)
      rbjointri <- rbind(rbjointri, peaload$alldata)
      
      # Order by pvalue
      rbjoinhtml = rbjoinhtml[order(rbjoinhtml$bonferroni_pvalue),] # htmltabla
      rbjointri = rbjointri[order(rbjointri$bonferroni_pvalue),] # alldata
      
      # To plot
      triucu = rbjointri[1:plotn,]
    }
  }
  
  
  
  if (organism == "mouse") {
    genes <- as.data.frame(genesetcompare)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)
    
    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genesetcompare = unique(onlygenes)
    
    genesr <- as.data.frame(geneset)
    names(genesr) <- "mouse_symbol"
    genesr$mouse_symbol <- as.character(genesr$mouse_symbol)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genesr, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genesr = unique(onlygenes)
    
    
    # Encontrar el no overlap entre unos y otros
    onlygenes <- foundgeneoverlap(genesr, genesetcompare)
    
    
    # Make two empty df with names
    rbjointri = data.frame(matrix(vector(), 0, 10,
                                  dimnames=list(c(), c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust",
                                                       "genes_associated_in_db","gene_overlap", "overlap_ratio"
                                                       ,"pvalue","common_genes"))),
                           stringsAsFactors=F)
    rbjoinhtml <-rbjointri
    
    # For each database name in a vector run its loadder
    for (db in database) {
      
      # Get the loadder name
      getnamedb <- paste0("",db,"loaddercompare")
      
      # Run the loadder
      peaload <- get(getnamedb)(onlygenes, genesetcompare,organism = "human", plotnumber = plotn, url = url)
      
      # join data frames generated
      rbjoinhtml <- rbind(rbjoinhtml, peaload$htmltabla)
      rbjointri <- rbind(rbjointri, peaload$alldata)
      
      # Order by pvalue
      rbjoinhtml = rbjoinhtml[order(rbjoinhtml$bonferroni_pvalue),] # htmltabla
      rbjointri = rbjointri[order(rbjointri$bonferroni_pvalue),] # alldata
      
      # To plot
      triucu = rbjointri[1:plotn,]
    }
    
    
    
    
    
  }
  
  newList <- list("alldata" = rbjointri, "htmltabla" = rbjoinhtml)
  return(newList)
}





#' @title PhenoEnrichGenespran
#' @description This function makes a plot with enriched phenotype info.
#' @param genes A vector with human gene symbol.
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @param plotn An integer with the number of term that you want in the plot. The default value is 30.
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export




# Función para obtener el enriquecimiento de términos HPO Y CRB
PhenoEnrichGenespran = function(genes, database = "HPO", organism = "human", plotn=30, plotnumber = 1:plotn){

  if (organism == "human") {

    # Si solo HPO
    if (database == "HPO") {


      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "HPO", organism)


      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))


      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))


      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))



      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]


    }

    # Si solo CRB
    if (database == "CRB") {
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))


      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(4,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]

    }

    # Si solo MGD
    if (database == "MGD"  | database == "MOUSEDB"){

      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")


      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "MGD",organism)


      # Make total_size_ratio per total
      numportertotalunir <- numportertotal

      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db,17895, length(genes), lower.tail = FALSE))



      triucusin = triucusin[order(triucusin$adjust_pvalue),]
      triucusin = triucusin[,c(1:2,6,3:5,7)]


      qglist <- QueryGenes(genes,database = "MGD", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]



    }

    # Si todas
    if (database == "ALL") {

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]


      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusinhpo = triucusin


      # Aqui empieza CRB

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))


      triucusin = triucusin[,c(1:2,6,3:5,7)]
      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(4,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusincrb = triucusin

      # Aqui empieza MGD
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")

      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "MGD",organism)



      numportertotalunir <- numportertotal


      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))





      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "MGD", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")


      triucusinmgd = triucusin

      triucusin = rbind(triucusinhpo,triucusincrb,triucusinmgd)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]


      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]


    }

    if (database == "CRBMGD") {



      # Aqui empieza CRB

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))


      triucusin = triucusin[,c(1:2,6,3:5,7)]
      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(4,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusincrb = triucusin

      # Aqui empieza MGD
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")

      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "MGD",organism)



      numportertotalunir <- numportertotal


      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))





      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "MGD", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")


      triucusinmgd = triucusin

      triucusin = rbind(triucusincrb,triucusinmgd)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]


      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]


    }

    if (database == "HPOMGD") {

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]


      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusinhpo = triucusin

      # Aqui empieza MGD
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")

      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "MGD",organism)



      numportertotalunir <- numportertotal


      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))





      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "MGD", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")


      triucusinmgd = triucusin

      triucusin = rbind(triucusinhpo,triucusinmgd)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]


      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]


    }

    if (database == "HPOCRB") {
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]
      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusinhpo = triucusin


      # Aqui empieza CRB

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_name")

      # Aqui guardo la fuente de datos de hpo
      triucusincrb = triucusin


      triucusin = rbind(triucusinhpo,triucusincrb)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]



    }

    # Si todas
    if (database == "HUMANDB") {
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]
      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusinhpo = triucusin


      # Aqui empieza CRB

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_name")

      # Aqui guardo la fuente de datos de hpo
      triucusincrb = triucusin


      triucusin = rbind(triucusinhpo,triucusincrb)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]



    }


    p <- triucu %>%

      # prepare text for tooltip
      dplyr::mutate(text = paste("ID Number: ", term_id, "\nTerm: ", term_name,  "\nGene Ratio: ", overlap_ratio, "\nAdjust pvalue: ", adjust_pvalue, "\nCommon genes: ", common_genes)) %>%
      # Classic ggplot
      ggplot(  ggplot2::aes(x= gene_overlap , y= reorder(paste0("",term_id,": ",term_name,""), -adjust_pvalue),  fill = adjust_pvalue, text=text)) +
      theme (text = element_text(size= 8)) +
      geom_bar(stat="identity", alpha=0.7) +
      scale_fill_gradientn(colours = rainbow(2), trans= "log") +
      theme(legend.position="right") +
      labs(title = "Phenotype Enrichment Analysis", x = "Number of genes associated with the term", y = "", fill = "p.value") + #
      theme(plot.title = element_text(size = 15),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            axis.text = element_text(
              angle = 0))

    # turn ggplot interactive with plotly layout.title
    nuevvogg <- ggplotly(p, tooltip="text") %>%
      config(mathjax = "cdn",
             toImageButtonOptions = list(
               format = "svg",
               filename = "myplot",
               width = 1000,
               height = 600
             )
      )



  }

  if (organism == "mouse") {
    genes <- as.data.frame(genes)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)

    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)

    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genes = onlygenes

    # Si solo HPO
    if (database == "HPO") {


      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "HPO", organism)


      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))


      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))


      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))



      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]


    }

    # Si solo CRB
    if (database == "CRB") {
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))


      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(4,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]

    }

    # Si solo MGD
    if (database == "MGD"  | database == "MOUSEDB"){

      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")


      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "MGD",organism)


      # Make total_size_ratio per total
      numportertotalunir <- numportertotal

      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db,17895, length(genes), lower.tail = FALSE))



      triucusin = triucusin[order(triucusin$adjust_pvalue),]
      triucusin = triucusin[,c(1:2,6,3:5,7)]


      qglist <- QueryGenes(genes,database = "MGD", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]



    }

    # Si todas
    if (database == "ALL") {

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]


      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusinhpo = triucusin


      # Aqui empieza CRB

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))


      triucusin = triucusin[,c(1:2,6,3:5,7)]
      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(4,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusincrb = triucusin

      # Aqui empieza MGD
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")

      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "MGD",organism)



      numportertotalunir <- numportertotal


      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))





      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "MGD", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")


      triucusinmgd = triucusin

      triucusin = rbind(triucusinhpo,triucusincrb,triucusinmgd)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]


      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]


    }

    if (database == "CRBMGD") {



      # Aqui empieza CRB

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))


      triucusin = triucusin[,c(1:2,6,3:5,7)]
      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(4,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusincrb = triucusin

      # Aqui empieza MGD
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")

      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "MGD",organism)



      numportertotalunir <- numportertotal


      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))





      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "MGD", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")


      triucusinmgd = triucusin

      triucusin = rbind(triucusincrb,triucusinmgd)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]


      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]


    }

    if (database == "HPOMGD") {

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]


      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusinhpo = triucusin

      # Aqui empieza MGD
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")

      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes, database = "MGD",organism)



      numportertotalunir <- numportertotal


      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))





      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 17895, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "MGD", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")


      triucusinmgd = triucusin

      triucusin = rbind(triucusinhpo,triucusinmgd)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]


      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]


    }

    # Si todas
    # Si todas
    if (database == "HUMANDB") {
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]
      qglist <- QueryGenes(genes,database = "HPO", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")

      # Aqui guardo la fuente de datos de hpo
      triucusinhpo = triucusin


      # Aqui empieza CRB

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumber(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, length(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]

      qglist <- QueryGenes(genes,database = "CRB", organism )
      urls <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",qglist$symbol,"")
      qglist$symbol <- paste0("[", qglist$symbol, "](", urls, ")")

      listadegenesqg <- qglist %>%
        group_by(term_id) %>%
        mutate(common_genes = paste0(symbol, collapse = ", "))


      listadegenesqg <- unique(listadegenesqg[,c(3,5)])

      triucusin <- inner_join(triucusin,listadegenesqg, by = "term_name")

      # Aqui guardo la fuente de datos de hpo
      triucusincrb = triucusin


      triucusin = rbind(triucusinhpo,triucusincrb)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      # Esto es para decidir el número de términos para plotear
      triucu = triucusin[plotnumber,]



    }


    p <- triucu %>%

      # prepare text for tooltip
      dplyr::mutate(text = paste("ID Number: ", term_id, "\nTerm: ", term_name,  "\nGene Ratio: ", overlap_ratio, "\nAdjust pvalue: ", adjust_pvalue, "\nCommon genes: ", common_genes)) %>%
      # Classic ggplot
      ggplot(  ggplot2::aes(x= gene_overlap , y= reorder(paste0("",term_id,": ",term_name,""), -adjust_pvalue),  fill = adjust_pvalue, text=text)) +
      theme (text = element_text(size= 8)) +
      geom_bar(stat="identity", alpha=0.7) +
      scale_fill_gradientn(colours = rainbow(2), trans= "log") +
      theme(legend.position="right") +
      labs(title = "Phenotype Enrichment Analysis", x = "Number of genes associated with the term", y = "", fill = "p.value") + #
      theme(plot.title = element_text(size = 15),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            axis.text = element_text(
              angle = 0))

    # turn ggplot interactive with plotly layout.title
    nuevvogg <- ggplotly(p, tooltip="text") %>%
      config(mathjax = "cdn",
             toImageButtonOptions = list(
               format = "svg",
               filename = "myplot",
               width = 1000,
               height = 600
             )
      )



  }

  newList <- list("alldata" = triucusin,   "graph" = nuevvogg)
  return(newList)
}



#' @title PhenoEnrichCompareRandom
#' @description This function makes a plot with enriched phenotype info.
#' @param genes A vector with human gene symbol.
#' @param database A string parameter with different options. You can define a vector with the databases.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @param plotn An integer with the number of term that you want in the plot. The default value is 30.
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export



PhenoEnrichCompareRandom = function(onlygenes, geneset, genesetcompare, database = "HPO", organism = "human", plotn=30, plotnumber = 1:plotn){
  
  
  if (organism == "human") {
    
    
    
    # Make two empty df with names
    rbjointri = data.frame(matrix(vector(), 0, 9,
                                  dimnames=list(c(), c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust",
                                                       "genes_associated_in_db","gene_overlap", "overlap_ratio"
                                                       ,"pvalue"))),
                           stringsAsFactors=F)
    
    # For each database name in a vector run its loadder
    for (db in database) {
      
      # Get the loadder name
      getnamedb <- paste0("",db,"loaddercompareOpt")
      
      # Run the loadder
      peaload <- get(getnamedb)(onlygenes, geneset, genesetcompare,organism = "human", plotnumber = plotn)
      
      # join data frames generated
      rbjointri <- rbind(rbjointri, peaload)
      
      # Order by pvalue
      rbjointri = rbjointri[order(rbjointri$bonferroni_pvalue),] # alldata
      
      # To plot
      triucu = rbjointri[1:plotn,]
    }
  }
  
  
  
  if (organism == "mouse") {
    genes <- as.data.frame(genesetcompare)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)
    
    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genesetcompare = unique(onlygenes)
    
    genesr <- as.data.frame(geneset)
    names(genesr) <- "mouse_symbol"
    genesr$mouse_symbol <- as.character(genesr$mouse_symbol)
    
    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genesr, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genesr = unique(onlygenes)
    
    
    # Encontrar el no overlap entre unos y otros
    onlygenes <- foundgeneoverlap(genesr, genesetcompare)
    
    
    # Make two empty df with names
    rbjointri = data.frame(matrix(vector(), 0, 9,
                                  dimnames=list(c(), c("term_id", "term_name", "source", "bonferroni_pvalue", "FDR_adjust",
                                                       "genes_associated_in_db","gene_overlap", "overlap_ratio"
                                                       ,"pvalue"))),
                           stringsAsFactors=F)
    
    # For each database name in a vector run its loadder
    for (db in database) {
      
      # Get the loadder name
      getnamedb <- paste0("",db,"loaddercompareOpt")
      
      # Run the loadder
      peaload <- get(getnamedb)(onlygenes, genesetcompare,organism = "human", plotnumber = plotn)
      
      # join data frames generated
      rbjointri <- rbind(rbjointri, peaload)
      
      # Order by pvalue
      rbjointri = rbjointri[order(rbjointri$bonferroni_pvalue),] # alldata
      
      # To plot
      triucu = rbjointri[1:plotn,]
    }
  }
  
  newList <- rbjointri
  return(newList)
}



#' @title foundgeneoverlap
#' @description This function give no overlapping genes
#' @param geneset A vector with human gene symbol.
#' @param genesetcompare A vector with human gene symbol.
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export




# Función para obtener el enriquecimiento de términos HPO Y CRB
foundgeneoverlap = function(geneset, genesetcompare){


    # Ver cuantos genes no hacen overlap para sacarlos del cálculo del target
    # Covert vector to a daframe
    referencegenesc <- as.data.frame(geneset)
    names(referencegenesc) <- "Symbol"

    # Pasar estos genes a entrezid clusterprofiler
    genes.df <- bitr(referencegenesc$Symbol, fromType = "SYMBOL",
                     toType = c("SYMBOL", "ENTREZID"),
                     OrgDb = "org.Hs.eg.db")
    referencegenesc <- unique(genes.df)
    names(referencegenesc) <- c("symbol","entrez")
    referencegenesc$entrez <- as.numeric(as.character(referencegenesc$entrez))


    targetgenesc <- as.data.frame(genesetcompare)
    names(targetgenesc) <- "Symbol"

    # Pasar estos genes a entrezid clusterprofiler
    genes.df <- bitr(targetgenesc$Symbol, fromType = "SYMBOL",
                     toType = c("SYMBOL", "ENTREZID"),
                     OrgDb = "org.Hs.eg.db")
    targetgenesc <- unique(genes.df)
    names(targetgenesc) <- c("symbol","entrez")
    targetgenesc$entrez <- as.numeric(as.character(targetgenesc$entrez))



    # We will keep only the selected genes
    onlygenes <- dplyr::anti_join(referencegenesc, targetgenesc, by = "entrez")




  return(onlygenes)
}





#



#' @title PhenoRelation
#' @description This function give no overlapping genes
#' @param genes A vector with human gene symbol.
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export

PhenoRelation = function(enrichresult, genes,  p.value = 0.05){


  # Tomar los fenotipos estadisticamente significativos
  relevantes <- dplyr::filter(enrichresult, enrichresult$adjust_pvalue <= p.value)


    p <- relevantes %>%  dplyr::mutate(db = case_when(stringr::str_detect(relevantes$term_id, "HP")  ~ "HPO",
                                                      stringr::str_detect(relevantes$term_id, "MP") ~ "MPO",
                                                      stringr::str_detect(relevantes$term_id, "CXXX") ~ "DIS",
                                                      stringr::str_detect(relevantes$term_id, "C\\d+") ~ "DIS",
                                                      stringr::str_detect(relevantes$term_id, "CRB")  ~ "CRB"))
    hpo <- dplyr::filter( p,  p$db == "HPO")
    hponr <- nrow(hpo)
    nhpo <-round(hponr * 0.1)
    hpo <- hpo[c(1:nhpo),]


    dis <- dplyr::filter( p,  p$db == "DIS")
    disnr <- nrow(dis)
    ndis <-round(disnr * 0.2)
    dis <- dis[c(1:ndis),]

    hpodis <- rbind(hpo,dis)
    hpodisnomb <- hpodis$term_id


  # Contar los fenotipos significativos
  numerpheno <- nrow(relevantes) # relevantes

  # Obtener el término de todos ellos
  nombresmatriz <- relevantes$term_id
  nombresfenotipos <- as.data.frame(nombresmatriz)
  names(nombresfenotipos) <- "term_id"

  # Crear una matriz con el número de fenotipos relevantes
  matpheno = matrix(data=NA, nrow= numerpheno, ncol=numerpheno)
  rownames(matpheno) <- nombresmatriz
  colnames(matpheno) <- nombresmatriz

  # Aqui ya empieza a depender de las bases de datos que estén marcadas, creo que esto puede ser siempre todas
  allquerygenes <- QueryGenes(genes, database = "ALL")


  # Mis genes
  genes <- as.data.frame(genes)
  names(genes) <- "Symbol"

  # Pasar estos genes a entrezid clusterprofiler
  genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = "org.Hs.eg.db")
  genes <- unique(genes.df$ENTREZID)
  genes <- as.data.frame(genes)
  names(genes) <- "entrez"
  genes$entrez <- as.numeric(as.character(genes$entrez))

  listadedf = list()
  for (i in 1:numerpheno) {

    # Leo el fenotipo
    phenotipo <- nombresmatriz[i]

    # Me quedo con todos los genes y fenotipos relacionados con esos de mis genes
    genesdeesefenotipo <- dplyr::filter(allquerygenes, allquerygenes$term_id == phenotipo)
    genesdeesefenotipo <- genesdeesefenotipo$entrez

    # De alguna forma parece que pierdo un gen del total?
    # Con estos genes ahora tengo que calcular el PhenoGeneNumber
    numerosdefenotiposhpo <- PhenoGeneNumberOpt(genesdeesefenotipo, database = "HPO")
    numerosdefenotiposmdg <- PhenoGeneNumberOpt(genesdeesefenotipo, database = "MGD")
    numerosdefenotiposdis <- PhenoGeneNumberOpt(genesdeesefenotipo, database = "DIS")
    numerosdefenotiposcrb <- PhenoGeneNumberOpt(genesdeesefenotipo, database = "CRB")

    # Junto los cálculos
    numerosdefenotipos <- rbind(numerosdefenotiposhpo,numerosdefenotiposdis,numerosdefenotiposmdg,numerosdefenotiposcrb)

    names(nombresfenotipos) <- "term_id"
    numerorelevantes <- inner_join(numerosdefenotipos,nombresfenotipos, by = "term_id")


    nopresentes <- anti_join(nombresfenotipos,numerorelevantes, by = "term_id")

    if (nrow(nopresentes)>0) {
      nopresentes$gen_associated_number <- 0
      nopresentes$term_id <- as.character(nopresentes$term_id)
      nopresentes$gen_associated_number <- as.integer(nopresentes$gen_associated_number)
      colnames(nopresentes) <- c("term_id","gene_overlap")
      numerorelevantes = numerorelevantes[,c(1,3)]

      resultadofenotipo <- rbind.data.frame(numerorelevantes,nopresentes)
      names(resultadofenotipo) <- c("term_id", paste0("",phenotipo,"") )

    }else{
      resultadofenotipo =numerorelevantes[,c(1,3)]
      names(resultadofenotipo) <- c("term_id", paste0("",phenotipo,"") )
    }



    listadedf[[i]] <-resultadofenotipo
  }

# Esto se lo come rapido
a <- Reduce(function(...) merge(..., by = "term_id", all=TRUE), listadedf)


  # Hacer que el data frame sea una matriz (Guardarla para enseñar porque tarda)
  mat1 = data.matrix(a)
  rownames(mat1) <-a$term_id
  mat1<-mat1[,-1]

  # Select mp overexpressed
  # select crb overexpressed

  matnew <-mat1[,colnames(mat1) %in% hpodisnomb]



  # El heatmap ya hace un cluster, podemos sacar ese cluster
  # COn el cluster tenemos que ver cuando caen cerca términos de diferentes bases de datos y mostrar eso
  # Tenemos colnames y rownames
  correlacion <- cor(mat1)

  matnew <-correlacion[rownames(correlacion) %in% hpodisnomb,colnames(correlacion) %in% hpodisnomb]

  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }

  res2<-rcorr(mat1)
  fmatda <- flattenCorrMatrix(res2$r, res2$P)

  heatmapfin <- pheatmap(matnew, cutree_rows = 4)

  newList <- list("cormatrix" = correlacion, "colum" = fmatda, "heatmap" = heatmapfin)
  return(newList)

}

#' @title PhenoGeneNumberSimuok
#' @description This function counts the number of genes per phenotype term from the vector of genes.
#' @param database A string parameter with different options: "ALL", "HUMANDB", "MOUSEDB", "HPO", "CRB" and "MGD"
#' @param genes A vector with human gene symbol.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export


PhenoGeneNumberSimuok = function(genes, database= "HPO", organism = "human"){

  if (organism == "human") {

    if (database == "HPO") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"


      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/hpobasededatos.csv"), header = T)
      names(custom) <- c("entrez","symbol", "term_id", "term_name")

      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "symbol")

      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")


    }

    if (database == "DIS") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"


      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/disbasededatos.csv"), header = T)
      custom <- custom[,c(1,2,5,6)]
      names(custom) <- c("entrez","symbol", "term_id", "term_name")

      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "symbol")

      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")


    }

    if (database == "CRB") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "symbol"





      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/cbrdatabaseok.csv"), header = T)
      names(custom) <- c("entrez","symbol", "term_name", "term_id")

      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "symbol")

      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
    }

    if (database == "MGD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "human_symbol"
      genes$HumanSymbol <- as.character(genes$human_symbol)



      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/mgibasefinal.tsv"), header = T, sep = "\t")
      names(custom) <- c("entrez","human_symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info"   )

      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "human_symbol")

      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name ) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")

    }
  }

  if (organism == "mouse") {

    # genes = musedata[1:50,3]

    genes <- as.data.frame(genes)
    names(genes) <- "mouse_symbol"
    genes$mouse_symbol <- as.character(genes$mouse_symbol)

    # Read database
    homologia <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                              "/homology.csv"), header = T)

    # We will keep only the selected genes
    onlygenes <- dplyr::inner_join(genes, homologia, by = "mouse_symbol")
    onlygenes <- onlygenes$human_symbol
    onlygenes <- as.data.frame(onlygenes)
    names(onlygenes) <- "symbol"
    onlygenes$symbol <- as.character(onlygenes$symbol)
    genes = onlygenes


    if (database == "HPO") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"

      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))



      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/hpobasededatos.csv"), header = T)
      names(custom) <- c("entrez","symbol", "term_id", "term_name")

      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")

      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
    }

    if (database == "CRB") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "Symbol"

      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$Symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))



      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/cbrdatabaseok.csv"), header = T)
      names(custom) <- c("entrez","symbol", "term_name", "term_id")

      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")

      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")
    }

    if (database == "MGD") {
      # Covert vector to a daframe
      genes <- as.data.frame(genes)
      names(genes) <- "human_symbol"
      genes$HumanSymbol <- as.character(genes$human_symbol)

      # Pasar estos genes a entrezid clusterprofiler
      genes.df <- bitr(genes$human_symbol, fromType = "SYMBOL",
                       toType = c("SYMBOL", "ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
      genes <- unique(genes.df$ENTREZID)
      genes <- as.data.frame(genes)
      names(genes) <- "entrez"
      genes$entrez <- as.numeric(as.character(genes$entrez))

      # Read database
      custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                             "/mgibasefinal.tsv"), header = T, sep = "\t")
      names(custom) <- c("entrez","human_symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info"   )

      # We will keep only the selected genes
      onlygenes <- dplyr::inner_join(custom, genes, by = "entrez")

      # We will count the number of genes per phenotype
      numportertotal <-   onlygenes %>% dplyr::group_by(term_id, term_name ) %>% dplyr::tally()
      numportergenes <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
      names(numportergenes) <- c("term_id", "term_name", "gene_overlap")

    }
  }

  return(numportergenes)
}

#' @title PhenoEnrichGenespransimu
#' @description This function analyses the similarity of phenotypic terms between two groups of genes
#' @param genes A vector with gene symbol, fixed genes or panel genes.
#' @param organism A string parameter: "human" if you have human genes and "mouse" if you have mouse genes
#' @param database A string parameter with different options: "ALL", "HUMANDB", "MOUSEDB", "HPO", "CRB" and "MGD"
#' @param p.value Cut of p value as they are relevant terms
#' @param nulltestnumber A number of random genes subset.
#' @param seed A seed.
#' @details This function has developed by Alejandro Cisterna García as part of his PhD mentored by Juan Antonio Botía Blaya
#' @export


PhenoEnrichGenespransimu = function(genes, database = "HPO", organism = "human"){


  if (organism == "human") {

    # Si solo HPO
    if (database == "HPO") {

      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumberSimuok(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, nrow(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, nrow(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]


    }

    # Si solo CRB
    if (database == "CRB") {
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumberSimuok(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, nrow(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, nrow(genes), lower.tail = FALSE))


      triucusin = triucusin[,c(1:2,6,3:5,7)]


    }

    # Si solo MGD
    if (database == "MGD"  | database == "MOUSEDB"){

      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")

      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumberSimuok(genes, database = "MGD",organism)



      numportertotalunir <- numportertotal


      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))





      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, nrow(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 17895, nrow(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]


    }

    if (database == "DIS") {

      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/disdatostotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumberSimuok(genes, database = "DIS", organism)


      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))


      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19722, nrow(genes), lower.tail = FALSE), method = "bonferroni")) %>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19722, nrow(genes), lower.tail = FALSE))


      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      triucusin = triucusin[,c(1:2,6,3:5,7)]

    }

    # Si todas
    if (database == "ALL") {

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumberSimuok(genes,database = "HPO",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))




      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248, nrow(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19248, nrow(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]


      # Aqui guardo la fuente de datos de hpo
      triucusinhpo = triucusin


      # Aqui empieza CRB

      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
      names(numportertotal) = c("term_id",   "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumberSimuok(genes,database = "CRB",organism)

      names(panel) = c( "term_id",   "term_name", "gene_overlap")

      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal



      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))






      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274, nrow(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19274, nrow(genes), lower.tail = FALSE))


      triucusin = triucusin[,c(1:2,6,3:5,7)]

      # Aqui guardo la fuente de datos de hpo
      triucusincrb = triucusin

      # Aqui empieza MGD
      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")

      names(numportertotal) <- c("term_id", "term_name", "genes_associated_in_db")

      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumberSimuok(genes, database = "MGD",organism)



      numportertotalunir <- numportertotal


      panelunir <- panel


      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))





      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)



      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895, nrow(genes), lower.tail = FALSE), method = "bonferroni"))%>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 17895, nrow(genes), lower.tail = FALSE))

      triucusin = triucusin[,c(1:2,6,3:5,7)]
      triucusinmgd = triucusin


      # Read data base
      numportertotal <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/disdatostotales.csv"), header = T)


      # We will count the number of genes per phenotype
      panel <- PhenoGeneNumberSimuok(genes, database = "DIS", organism)


      # Make total_size_ratio per total
      #numportertotalunir <- numportertotal[,c(1,2,4)]
      numportertotalunir <- numportertotal


      #panelunir <- panel[,c(1,2,4)]
      panelunir <- panel

      names(numportertotalunir) = c( "term_id",   "term_name", "genes_associated_in_db")

      triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))


      # Esto es para que almenos tenga anotados 10 genes en ese fenotipo
      triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)


      triucusin <- triucu

      # Add gene total_size_ratio col
      triucusin <- triucusin %>%
        dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))


      # Calculate p.value
      triucusin <- triucusin %>%
        dplyr::mutate(adjust_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19722, nrow(genes), lower.tail = FALSE), method = "bonferroni")) %>%
        dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, 19722, nrow(genes), lower.tail = FALSE))

      triucusin = triucusin[order(triucusin$adjust_pvalue),]

      triucusin = triucusin[,c(1:2,6,3:5,7)]



      triucusindis = triucusin



      triucusin = rbind(triucusinhpo,triucusincrb,triucusinmgd,triucusindis)


      triucusin = triucusin[order(triucusin$adjust_pvalue),]




    }


  }

  return(triucusin)

}



#' @title HPOloadder
#' @export



HPOloadder = function(genes, organism = "human", plotnumber = 30, url = F){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)
  # Number of genes for HPO
  ngeneshpo = 19248
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "HPO", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, ngeneshpo - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, ngeneshpo - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  
  
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "HPO", organism )
  
  if (url ) {
    
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
    
  }
  # Get the url for gene
  
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique names
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  #listadegenesqg[,c(3,6)]
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url ) {
    htmlttriucu$term_id <- paste0("<a href='https://hpo.jax.org/app/browse/term/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  
  # Return list
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


#' @title MGDloadder
#' @export
MGDloadder = function(genes, organism = "human", plotnumber = 30, url = F){
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")
  ngenesmgd = 17895
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "MGD", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, ngenesmgd - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, ngenesmgd - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "MGD", organism )
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,5)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='http://www.informatics.jax.org/vocab/mp_ontology/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}

#' @title CRBloadder
#' @export
CRBloadder = function(genes, organism = "human", plotnumber = 30, url = F){
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "CRB", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19274  - genes_associated_in_db , length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "CRB", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
    
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(4,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://crisprbrain.org/simple-screen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}

#' @title CLINGENloadder
#' @export
CLINGENloadder = function(genes, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/clingendatostotales.csv"), header = T)
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "CLINGEN", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19198  - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19198 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "CLINGEN", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


#' @title GENOMICS_ENGLANDloadder
#' @export
GENOMICS_ENGLANDloadder = function(genes, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/GENOMICS_ENGLANDdatostotales.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "GENOMICS_ENGLAND", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19230 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19230 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "GENOMICS_ENGLAND", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}

#' @title CTDloadder
#' @export
CTDloadder = function(genes, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/CTD_humandatostotales.csv"), header = T)
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "CTD", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19636 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19636 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "CTD", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}

#' @title ORPHANETloadder
#' @export
ORPHANETloadder = function(genes, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/ORPHANETdatostotales.csv"), header = T)
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "ORPHANET", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19262 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19262 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "ORPHANET", organism )
  
  # Get the url for gene
  if (url) {
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}

#' @title PSYGENETloadder
#' @export
PSYGENETloadder = function(genes, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/PSYGENETdatostotales.csv"), header = T)
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "PSYGENET", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19262 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19262 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "PSYGENET", organism )
  
  # Get the url for gene
  if (url) {
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}

#' @title CGIloadder
#' @export
CGIloadder = function(genes, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/CGIdatostotales.csv"), header = T)
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "CGI", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19198 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19198 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "CGI", organism )
  
  # Get the url for gene
  if (url) {
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}

#' @title UNIPROTloadder
#' @export
UNIPROTloadder = function(genes, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/UNIPROTdatostotales.csv"), header = T)
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genes, database = "UNIPROT", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19204 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19204 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genes,database = "UNIPROT", organism )
  
  # Get the url for gene
  if (url) {
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


#' @title HPOloadderOpt
#' @export

HPOloadderOpt = function(genes, organism = "human", plotnumber = 30){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/hpodatatotales.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "HPO", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19248 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19248 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
}

#' @title MGDloadderOpt
#' @export
MGDloadderOpt = function(genes, organism = "human", plotnumber = 30){
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/totalmponumber.tsv"), header = T, sep = "\t")
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "MGD", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 17895 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 17895 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
}

#' @title CRBloadderOpt
#' @export

CRBloadderOpt = function(genes, organism = "human", plotnumber = 30){
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/cbrtotalnumber.csv"), header = T)
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "CRB", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19274 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19274 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
}

#' @title CLINGENloadderOpt
#' @export


CLINGENloadderOpt = function(genes, organism = "human", plotnumber = 30){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/clingendatostotales.csv"), header = T)
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "CLINGEN", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19198 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19198 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
}

#' @title GENOMICS_ENGLANDloadderOpt
#' @export


GENOMICS_ENGLANDloadderOpt = function(genes, organism = "human", plotnumber = 30){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/GENOMICS_ENGLANDdatostotales.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "GENOMICS_ENGLAND", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19230 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19230 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
}


#' @title CTDloadderOpt
#' @export

CTDloadderOpt = function(genes, organism = "human", plotnumber = 30){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/CTD_humandatostotales.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "CTD", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19636 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19636 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
}


#' @title ORPHANETloadderOpt
#' @export


ORPHANETloadderOpt = function(genes, organism = "human", plotnumber = 30){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/ORPHANETdatostotales.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "ORPHANET", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19262 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19262 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
}

#' @title PSYGENETloadderOpt
#' @export

PSYGENETloadderOpt = function(genes, organism = "human", plotnumber = 30){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/PSYGENETdatostotales.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "PSYGENET", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19262 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19262 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
}


#' @title UNIPROTloadderOpt
#' @export
UNIPROTloadderOpt = function(genes, organism = "human", plotnumber = 30){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/UNIPROTdatostotales.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "UNIPROT", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19204 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19204 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
}



#' @title CGIloadderOpt
#' @export

CGIloadderOpt = function(genes, organism = "human", plotnumber = 30){
  
  # Read data base
  numportertotalunir <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                                     "/CGIdatostotales.csv"), header = T)
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genes, database = "CGI", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, 19198 - genes_associated_in_db, length(genes), lower.tail = FALSE), method = "bonferroni")) %>%
    dplyr::mutate(pvalue = phyper(gene_overlap -1, genes_associated_in_db, 19198 - genes_associated_in_db, length(genes), lower.tail = FALSE))
  
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
}

 #' @title HPOloaddercompare
 #' @export
HPOloaddercompare = function(onlygenes, genesetcompare, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/hpobasededatos.csv"), header = T)
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "HPO"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "HPO", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare, database = "HPO", organism )
  
  if (url) {
    
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
    
  }
  
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique names
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://hpo.jax.org/app/browse/term/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  
  # Return list
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


 #' @title MGDloaddercompare
 #' @export


MGDloaddercompare = function(onlygenes,genesetcompare, organism = "human", plotnumber = 30,  url = F){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/mgibasefinal.tsv"), header = T, sep = "\t")
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "MGD"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "MGD",organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene total_size_ratio col
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$entrez)), length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$entrez)), length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "MGD", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,5)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='http://www.informatics.jax.org/vocab/mp_ontology/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}



 #' @title CRBloaddercompare
 #' @export


CRBloaddercompare = function(onlygenes, genesetcompare, organism = "human", plotnumber = 30,  url = F){
  # Read data base
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/cbrdatabaseok.csv"), header = T)
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "CRB"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare,database = "CRB",organism)
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "CRB", organism )
  
  # Get the url for gene
  if (url) {
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(4,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://crisprbrain.org/simple-screen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


 #' @title CLINGENloaddercompare
 #' @export

CLINGENloaddercompare = function(onlygenes,genesetcompare, organism = "human", plotnumber = 30,  url = F){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/clingenbasededatos.csv"), header = T)
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "CLINGEN"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "CLINGEN", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "CLINGEN", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}



 #' @title GENOMICS_ENGLANDloaddercompare
 #' @export

GENOMICS_ENGLANDloaddercompare = function(onlygenes,genesetcompare, organism = "human", plotnumber = 30,  url = F){
  
  # Read data base
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/GENOMICS_ENGLANDbasededatos.csv"), header = T)
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "GENOMIC_ENGLAND"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "GENOMICS_ENGLAND", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "GENOMICS_ENGLAND", organism )
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


 #' @title CTDloaddercompare
 #' @export

CTDloaddercompare = function(onlygenes,genesetcompare, organism = "human", plotnumber = 30,  url = F){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/CTD_humanbasededatos.csv"), header = T)
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "CTD"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "CTD", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "CTD", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


 #' @title ORPHANETloaddercompare
 #' @export

ORPHANETloaddercompare = function(onlygenes,genesetcompare, organism = "human", plotnumber = 30,  url = F){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/ORPHANETbasededatos.csv"), header = T)
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  
  numportertotalunir$source <- "ORPHANET"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "ORPHANET", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "ORPHANET", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


 #' @title PSYGENETloaddercompare
 #' @export

PSYGENETloaddercompare = function(onlygenes,genesetcompare, organism = "human", plotnumber = 30,  url = F){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/PSYGENETbasededatos.csv"), header = T)
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "PSYGENET"
  
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "PSYGENET", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "PSYGENET", organism )
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}


 #' @title UNIPROTloaddercompare
 #' @export



UNIPROTloaddercompare = function(onlygenes,genesetcompare, organism = "human", plotnumber = 30,  url = F){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/UNIPROTbasededatos.csv"), header = T)
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  
  numportertotalunir$source <- "UNIPROT"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "UNIPROT", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "UNIPROT", organism )
  
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}




 #' @title CGIloaddercompare
 #' @export

CGIloaddercompare = function(onlygenes,genesetcompare, organism = "human", plotnumber = 30,  url = F){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/CGIbasededatos.csv"), header = T)
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  
  numportertotalunir$source <- "CGI"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumber(genesetcompare, database = "CGI", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue= p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
  # Get each gene per phenotype
  qglist <- QueryGenes(genesetcompare,database = "CGI", organism )
  if (url) {
    # Get the url for gene
    qglist$symbol <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", qglist$symbol,"'>", qglist$symbol,"</a>")
  }
  # Collapse genes per phenotype
  listadegenesqg <- qglist %>%
    group_by(term_id) %>%
    mutate(common_genes = paste0(symbol, collapse = ", "))
  
  # Unique
  listadegenesqg <- unique(listadegenesqg[,c(3,6)])
  
  # Join genes and dataframe
  triucusin <- inner_join(triucusin,listadegenesqg, by = "term_id")
  
  # Esto es para decidir el número de términos para plotear
  triucu = triucusin[1:plotnumber,]
  
  
  # Make interactive table with phenotype url
  htmlttriucu <- triucusin
  if (url) {
    htmlttriucu$term_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/medgen/", htmlttriucu$term_id,"'>", htmlttriucu$term_id,"</a>")
  }
  newList <- list("alldata" = triucusin, "htmltabla" = htmlttriucu)
  return(newList)
}

 #' @title HPOloaddercompareOpt
 #' @export

HPOloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  
  # Read data base
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/hpobasededatos.csv"), header = T)
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "HPO"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "HPO", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  return(triucusin)
}


 #' @title MGDloaddercompareOpt
 #' @export
MGDloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/mgibasefinal.tsv"), header = T, sep = "\t")
  
  names(custom) <- c("entrez","human_symbol", "mouse_symbol", "mgi", "term_id", "term_name", "info", "source"   )
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "MGD"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "MGD",organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene total_size_ratio col
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$entrez)), length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$entrez)), length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  return(triucusin)
}




 #' @title CRBloaddercompareOpt
 #' @export

CRBloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  # Read data base
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/cbrdatabaseok.csv"), header = T)
  
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "CRB"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare,database = "CRB",organism)
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save phenotypes with at least 10 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 10)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
  return(triucusin)
}



 #' @title CLINGENloaddercompareOpt
 #' @export
CLINGENloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/clingenbasededatos.csv"), header = T)
  
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "CLINGEN"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "CLINGEN", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
  return(triucusin)
}


 #' @title GENOMICS_ENGLANDloaddercompareOpt
 #' @export

GENOMICS_ENGLANDloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  
  # Read data base
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/GENOMICS_ENGLANDbasededatos.csv"), header = T)
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "GENOMICS_ENGLAND"
  
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "GENOMICS_ENGLAND", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  return(triucusin)
}


 #' @title CTDloaddercompareOpt
 #' @export

CTDloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/CTD_humanbasededatos.csv"), header = T)
  
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "CTD"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "CTD", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  return(triucusin)
}


 #' @title ORPHANETloaddercompareOpt
 #' @export

ORPHANETloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/ORPHANETbasededatos.csv"), header = T)
  
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "ORPHANET"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "ORPHANET", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
  return(triucusin)
}


 #' @title PSYGENETloaddercompareOpt
 #' @export

PSYGENETloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/PSYGENETbasededatos.csv"), header = T)
  
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "PSYGENET"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "PSYGENET", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
  return(triucusin)
}


 #' @title UNIPROTloaddercompareOpt
 #' @export
UNIPROTloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/UNIPROTbasededatos.csv"), header = T)
  
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "UNIPROT"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "UNIPROT", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
  return(triucusin)
}



 #' @title CGIloaddercompareOpt
 #' @export
CGIloaddercompareOpt = function(onlygenes, geneset, genesetcompare, organism = "human", plotnumber = 30){
  custom <- fread(paste0(path.package("PhenoExam", quiet = FALSE),
                         "/CGIbasededatos.csv"), header = T)
  
  
  
  custom <- dplyr::anti_join(custom, onlygenes, by = "entrez")
  
  # We will count the number of genes per phenotype
  numportertotal <-   custom %>% dplyr::group_by(term_id, term_name) %>% dplyr::tally()
  numportertotalunir <- dplyr::arrange(numportertotal, dplyr::desc(numportertotal$n))
  names(numportertotalunir) <- c("term_id", "term_name", "genes_associated_in_db")
  numportertotalunir$source <- "CGI"
  
  # We will count the number of genes per phenotype
  panelunir <- PhenoGeneNumberOpt(genesetcompare, database = "CGI", organism)
  
  # We join the dataframes
  triucu <- merge.data.frame(numportertotalunir,panelunir, by = c("term_id","term_name"))
  
  # Save diseases terms with at least 5 genes
  triucu <- dplyr::filter(triucu, triucu$genes_associated_in_db >= 5)
  
  # To have two DF with different names
  triucusin <- triucu
  
  # Add gene overlap_ratio
  triucusin <- triucusin %>%
    dplyr::mutate(overlap_ratio = paste0("(",gene_overlap,"/",genes_associated_in_db,")"))
  
  
  # Calculate p.value (genes in each database)
  # Calculate p.value
  triucusin <- triucusin %>%
    dplyr::mutate(bonferroni_pvalue = p.adjust(phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE), method = "bonferroni"))%>%
    dplyr::mutate(pvalue = phyper(gene_overlap - 1, genes_associated_in_db, length(unique(custom$symbol)) - genes_associated_in_db, length(genesetcompare), lower.tail = FALSE))
  
  triucusin <- triucusin %>%
    dplyr::mutate(FDR_adjust= p.adjust(p = triucusin$pvalue, method = "fdr"))
  # Order by pvalue
  triucusin = triucusin[order(triucusin$bonferroni_pvalue),]
  triucusin = triucusin[,c(1:2,4,7,9,3,5,6,8)]
  
  
  return(triucusin)
}


#' @title forbes
#' @export
forbes<-function(x,y,corrected)	{
  if (missing(corrected))	{
    corrected <- T
  }
  if (is.numeric(x) && is.numeric(y) && min(x) == 0 && min(y) == 0 && length(x) == length(y))	{
    a <- length(which((x * y) > 0))
    b <- length(which(x > 0)) - a
    c <- length(which(y > 0)) - a
  } else	{
    a <- length(na.omit(match(x,y)))
    b <- length(x) - a
    c <- length(y) - a
  }
  n <- a + b + c
  if (corrected == T)	{
    return(a * (n + sqrt(n))/((a + b) * (a + c) + a * sqrt(n) + (b * c)/2))
  } else	{
    return(a * n/((a + b) * (a + c)))
  }
}


