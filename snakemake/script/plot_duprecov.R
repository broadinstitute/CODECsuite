library(ggplot2)
fp = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/ffpe_wes_cds/pip/wrk/metrics/ffpe_tumor.duplex_yield_metrics.txt"
df = fread(reformat_on_mac(fp))
ggplot(df, aes(x=read_pairs, y=ds_duplexes)) + geom_line() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(28)

s = 0
for (i in 1:7) {
  s = s + (1/i)^2 
}
s

pi^2/6 - 1.511797


