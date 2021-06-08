capture_v2 = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/dec4_pancancer_v2/wrk/metrics/byproduct")
by_product <- function(inputdir) {
  dfall = data.frame()
  for (f in list.files(inputdir, "byproduct.txt")) {
    tmp = unlist(strsplit(f, "\\."))
    sampleid = tmp[3]
    laneid = tmp[1]
    laneid = str_remove(laneid, "_split")
    df = fread(file.path(inputdir, f)) 
    #colnames(df) = c(names(df)[2:29], "V1")
    df$sample_id = sampleid
    df$lane = laneid
    dfall = rbind(dfall, df)
  }
  dfall
}
capture_df = by_product(capture_v2)
sum(capture_df$n_cds_highconf) / sum(capture_df$n_total)

capture_df %>% group_by(sample_id, lane) %>% summarise(n_cds_highconf = sum(n_cds_highconf), 
                                      n_low_conf = sum(n_cds_lowconf, n_single_lowconf), 
                                      n_singleinsert = sum(n_single_highconf), 
                                      n_double_ligation = sum(n_double_ligation),
                                      n_adp_dimer = sum(n_adp_dimer),
                                      n_intermol = sum(n_intermol),
                                      n_unmapped = sum(n_unmapped),
                                      n_insuf_trim = sum(n_insuf_trim), 
                                      n_total = sum(n_total)) %>% mutate (highconf = n_cds_highconf/n_total, lowconf=n_low_conf / n_total, doubleligation = n_double_ligation /n_total,
                                                                          singleinsert = n_singleinsert / n_total, adp_dimer = n_adp_dimer / n_total, 
                                                                          intermol = n_intermol / n_total, unmapped = n_unmapped /n_total, longinsert = n_insuf_trim /n_total) %>%
                                      select(sample_id, lane, n_cds_highconf, highconf, doubleligation, adp_dimer, intermol, unmapped) %>% write_csv(., path="Dec4.v2_byproduct.csv")