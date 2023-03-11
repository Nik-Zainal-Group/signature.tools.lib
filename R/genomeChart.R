
#' @export
genomeChart <- function(outfilename,
                        sample_name,
                        SNV_vcf_file = NULL,
                        SNV_tab_file = NULL,
                        Indels_vcf_file = NULL,
                        Indels_tab_file = NULL,
                        CNV_tab_file = NULL,
                        SV_bedpe_file = NULL,
                        plot_title = NULL,
                        runKataegis = TRUE,
                        genome.v = "hg19"){
  
  # set up some variables
  snvs_table <- NULL
  sbs_obj <- NULL
  indels_obj <- NULL
  CNV_table <- NULL
  sv_obj <- NULL
  
  # Loading SNVs if available
  if(!is.null(SNV_vcf_file)){
    if (file.exists(SNV_vcf_file)){
      snvs_table <- fromVcfToTable(vcfFilename = SNV_vcf_file,
                                   genome.v = genome.v)
    }else{
      message("[warning genomeChart] SNV_vcf_file file not found: ",SNV_vcf_file,". Ignoring and moving on.")
      SNV_vcf_file <- NULL
    }
  }
  if(!is.null(SNV_tab_file)){
    if (file.exists(SNV_tab_file)){
      if(is.null(snvs_table)){
        snvs_table <- read.table(file = SNV_tab_file,sep = "\t",header = TRUE,
                                 check.names = FALSE,stringsAsFactors = FALSE)
      }else{
        message("[warning genomeChart] SNV_tab_file ignored because SNVs already loaded from SNV_vcf_file.")
        SNV_tab_file <- NULL
      }
    }else{
      message("[warning genomeChart] SNV_tab_file file not found: ",SNV_tab_file,". Ignoring and moving on.")
      SNV_tab_file <- NULL
    }
  }
  
  # Loading Indels if available
  if(!is.null(Indels_vcf_file)){
    if (file.exists(Indels_vcf_file)){
      indels_obj <- vcfToIndelsClassification(indelsVCF.file = Indels_vcf_file,
                                              sampleID = sample_name,
                                              genome.v = genome.v)
    }else{
      message("[warning genomeChart] Indels_vcf_file file not found: ",Indels_vcf_file,". Ignoring and moving on.")
      Indels_vcf_file <- NULL
    }
  }
  if(!is.null(Indels_tab_file)){
    if (file.exists(Indels_tab_file)){
      if(is.null(indels_obj)){
        indels_obj <- tabToIndelsClassification(indel.data = read.table(file = Indels_tab_file,sep = "\t",header = TRUE,
                                                                        check.names = FALSE,stringsAsFactors = FALSE),
                                                sampleID = sample_name,
                                                genome.v = genome.v)
      }else{
        message("[warning genomeChart] Indels_tab_file ignored because Indels already loaded from Indels_vcf_file")
        Indels_tab_file <- NULL
      }
    }else{
      message("[warning genomeChart] Indels_tab_file file not found: ",Indels_tab_file,". Ignoring and moving on.")
      Indels_tab_file <- NULL
    }
  }
  
  # Loading CNV file
  if(!is.null(CNV_tab_file)){
    if (file.exists(CNV_tab_file)){
      CNV_table <- read.table(file = CNV_tab_file,sep = "\t",header = TRUE,
                              check.names = FALSE,stringsAsFactors = FALSE)
    }else{
      message("[warning genomeChart] CNV_tab_file file not found: ",CNV_tab_file,". Ignoring and moving on.")
      CNV_tab_file <- NULL
    }
  }
  
  # Loading SV bedpe file
  if(!is.null(SV_bedpe_file)){
    if (file.exists(SV_bedpe_file)){
      sv_obj <- bedpeToRearrCatalogue(sv_bedpe = readTable(file = SV_bedpe_file))
    }else{
      message("[warning genomeChart] SV_bedpe_file file not found: ",SV_bedpe_file,". Ignoring and moving on.")
      SV_bedpe_file <- NULL
    }
  }
  
  # check if there is anything at all to plot
  if(is.null(SNV_vcf_file) & is.null(SNV_tab_file) & is.null(Indels_vcf_file) & is.null(Indels_tab_file) & is.null(CNV_tab_file) & is.null(SV_bedpe_file)){
    message("[error genomeChart] no data found, nothing to plot. Quit.")
    return(NULL)
  }
  
  # get SBS catalogue if possible
  if(!is.null(snvs_table)){
    sbs_obj <- tabToSNVcatalogue(subs = snvs_table,
                                 genome.v = genome.v)
  }
  
  # check SBS catalogues in SV clustering regions
  clustering_regions_sbs_catalogues <- NULL
  clusteringSBScatalogue_all <- NULL
  if(!is.null(sv_obj$clustering_regions) & !is.null(snvs_table)){
    # now find snvs in SV clusters and build catalogues
    resSBSinSVclusters <- getSBScataloguesInSVclusters(clustering_regions = sv_obj$clustering_regions,
                                                       snvs_table = snvs_table,
                                                       sample_name = sample_name,
                                                       genome.v = genome.v)
    clustering_regions_sbs_catalogues <- resSBSinSVclusters$clustering_regions_sbs_catalogues
    clusteringSBScatalogue_all <- resSBSinSVclusters$clusteringSBScatalogue_all
    snvs_table <- resSBSinSVclusters$snvs_table
  }
  
  # run the kataegis algorithm if requested
  kataegis <- NULL
  if(!is.null(snvs_table) & runKataegis){
    kataegis <- findKataegis(snvs_table = snvs_table,
                             sample_name = sample_name)
  }
}

