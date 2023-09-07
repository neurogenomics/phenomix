get_cs2g_gwascatalog <- function(dat){
    ID <- NULL
    #### Assign unique trait IDs to differentiate them after truncation ####
    dat[, trait_id := .GRP, by = .(DISEASE.TRAIT, PUBMEDID)]
    data.table::setnames(dat,"POS","BP")
    #### Remove non-ASCI characters ####
    dat[, ID :=
            paste(
                stringr::str_trunc(
                    iconv(
                        paste(
                            gsub(" ", "-", DISEASE.TRAIT),
                            PUBMEDID,
                            sep = "_"
                        ),
                        "latin1", "ASCII",
                        sub = ""
                    ),
                    width = 40
                ), # MOFA2 doesn't let you have IDs > 50 characters
                trait_id,
                sep = "_"
            )]
    #### Create metadata ####
    obs <- dat[, .(
        DISEASE.TRAIT = unique(DISEASE.TRAIT),
        PUBMEDID = unique(PUBMEDID),
        N_SNP = data.table::uniqueN(SNP),
        N_GENE = data.table::uniqueN(gene),
        cS2G_mean = mean(cS2G)
    ),
    keyby = "ID"
    ] 
    obs <- data.frame(obs,row.names = obs$ID)
    return(list(data=dat,
                obs=obs))
}
