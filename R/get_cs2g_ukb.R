get_cs2g_ukb <- function(dat){
    GROUP <- TRAIT <- NULL;
    
    #### Prepare metadata ####
    dat[, ID := gsub(" ", "-", DISEASE.TRAIT)]
    #### Remove non-ASCI characters ####
    dat[, ID := iconv(ID, "latin1", "ASCII", sub = "")]
    #### Create metadata ####
    obs <- dat[, .(
        DISEASE.TRAIT = unique(DISEASE.TRAIT),
        N_SNP = data.table::uniqueN(SNP),
        N_GENE = data.table::uniqueN(gene),
        cS2G_mean = mean(cS2G = TRUE),
        PIP_mean = mean(PIP, na.rm = TRUE)
    ),
    keyby = c("ID")
    ]
    obs[, GROUP := stringr::str_split(DISEASE.TRAIT, "_",
                                      n = 2, simplify = TRUE
    )[, 1]]
    obs[, TRAIT := stringr::str_split(DISEASE.TRAIT, "_",
                                      n = 2, simplify = TRUE
    )[, 2]]
    obs <- data.frame(obs,row.names = obs$ID)
    return(list(data=dat,
                obs=obs))
}