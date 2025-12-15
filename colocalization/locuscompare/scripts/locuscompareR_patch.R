#!/us/bin/env Rscript

library(httr)
library(jsonlite)

#get_position <- function(x, genome = c("hg19", "hg38")) {
#    stopifnot("rsid" %in% colnames(x))
#
#    genome <- match.arg(genome)
#
#    # hg19 = GRCh37 を指定する
#    assembly <- ifelse(genome == "hg19", "GRCh37", "GRCh38")
#
#    fetch_one <- function(rs) {
#        url <- sprintf(
#            "https://rest.ensembl.org/variation/human/%s?content-type=application/json",
#            rs
#        )
#
#        res <- tryCatch(
#            httr::GET(url, timeout(5)),
#            error = function(e) return(NULL)
#        )
#
#        if (is.null(res) || httr::status_code(res) != 200) {
#            return(data.frame(rsid = rs, chr = NA, pos = NA))
#        }
#
#        js <- jsonlite::fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))
#
#        # mappings の中に chr, pos がある
#        if (!"mappings" %in% names(js)) {
#            return(data.frame(rsid = rs, chr = NA, pos = NA))
#        }
#
#        # 指定した assembly の位置を取得
#        m <- js$mappings
#        m <- m[m$assembly_name == assembly, ]
#
#        if (nrow(m) == 0) {
#            return(data.frame(rsid = rs, chr = NA, pos = NA))
#        }
#
#        return(data.frame(rsid = rs, chr = m$seq_region_name[1], pos = m$start[1]))
#    }
#
#    # 全 rsid に対して API コール
#    out <- do.call(
#        rbind,
#        lapply(x$rsid, fetch_one)
#    )
#
#    # merge
#    merged <- merge(x, out, by = "rsid", all.x = TRUE)
#
#    return(merged)
#}


get_position <- function(x, genome = c("hg19", "hg38")) {
    stopifnot("rsid" %in% colnames(x))

    genome <- match.arg(genome)
    assembly <- ifelse(genome == "hg19", "GRCh37", "GRCh38")

    # --- バッチ API ---
    fetch_batch <- function(rsids) {
        url <- "https://rest.ensembl.org/variation/homo_sapiens"

        res <- httr::POST(
            url,
            body = list(ids = rsids),
            encode = "json",
            httr::add_headers("Content-Type" = "application/json",
                              "Accept" = "application/json")
        )

        if (httr::status_code(res) != 200) {
            warning("Batch request failed")
            return(data.frame(rsid = rsids, chr = NA, pos = NA))
        }

        js <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))

        # rsID ごとに整形
        out <- lapply(names(js), function(id) {
            entry <- js[[id]]

            if (is.null(entry$mappings)) {
                return(data.frame(rsid = id, chr = NA, pos = NA))
            }

            m <- entry$mappings
            m <- m[m$assembly_name == assembly, ]
            if (nrow(m) == 0) {
                return(data.frame(rsid = id, chr = NA, pos = NA))
            }

            data.frame(
                rsid = id,
                chr = m$seq_region_name[1],
                pos = m$start[1]
            )
        })

        do.call(rbind, out)
    }

    # --- 200件ずつバッチ処理 ---
    rsids <- x$rsid
    batches <- split(rsids, ceiling(seq_along(rsids) / 200))
    n_batch <- length(batches)

    ## シンプルver.
    #batch_results <- do.call(
    #    rbind,
    #    lapply(batches, fetch_batch)
    #)

    # 進捗表示ver.
    all_results <- list()

    for (i in seq_along(batches)) {
        cat(sprintf("Processing batch %d / %d (size = %d)\n",
                    i, n_batch, length(batches[[i]])))
        flush.console()

        all_results[[i]] <- fetch_batch(batches[[i]])
    }

    batch_results <- do.call(rbind, all_results)

    # --- merge with original ---
    merged <- merge(x, batch_results, by = "rsid", all.x = TRUE)

    return(merged)
}


retrieve_LD <- function(chr, snp, population = "EUR") {

    pop_map <- list(
        EUR = "1000GENOMES:phase_3:EUR",
        EAS = "1000GENOMES:phase_3:EAS",
        AFR = "1000GENOMES:phase_3:AFR",
        SAS = "1000GENOMES:phase_3:SAS",
        AMR = "1000GENOMES:phase_3:AMR",
        GBR = "1000GENOMES:phase_3:GBR"
    )

    if (!(population %in% names(pop_map))) {
        stop("Population must be one of: EUR, EAS, AFR, SAS, AMR, GBR")
    }

    pop_code <- pop_map[[population]]

    url <- sprintf(
        "https://rest.ensembl.org/ld/human/%s/%s?content-type=application/json&population=%s",
        snp, pop_code, pop_code
    )

    res <- jsonlite::fromJSON(url)

    if (length(res) == 0) {
        warning(sprintf("No LD returned for %s", snp))
        return(data.frame(SNP_A = character(), SNP_B = character(), R2 = numeric()))
    }

    df <- data.frame(
        SNP_A = snp,
        SNP_B = res$variation2,
        R2    = as.numeric(res$r2),
        stringsAsFactors = FALSE
    )

    df <- rbind(
        df,
        data.frame(
            SNP_A = df$SNP_B,
            SNP_B = df$SNP_A,
            R2 = df$R2
        )
    )
    ## 自分自身の行を除外
    #df <- df[df$SNP_B != snp, ]
    #
    ## R² > 0.2 の LD のみ（元関数と同じ）
    #df <- df[df$R2 > 0.2, ]

    return(df)
}
