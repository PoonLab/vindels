# README (within host figures)

## ins-plot-jan23.png & del-plot-jan23.png
Plot comparing nucleotide compositions between insertions with lengths a multiple of 3 vs those not a multiple of 3.
Error bars represent the 95% CIs.


## nt-indels-all
Error bars are 95% CIs.

## pnlgs-indel.png

Net changes in PNLG sites due to insertions and deletions. Error bars represent 95% CIs. Top legend for Insertions vs Deletions has been finicky. 

## indel-rates.png
_Revising this figure._ Error bars are currently the 95% CIs of the 20 replicates, but I plan to fix these with the new data. 

## indel-rtt-midpoints.png
This plot describes the root-to-tip midpoints of all tree edges that contained an insertion or deletion. This included all patients. Error bars are the 50% quantiles. Not sure if we should include this given the enormous uncertainty. 

## rates-ins-tipvnode.png  & rates-del-tipvnode.png
No error bars shown. I have a bootstrapping algorithm set up but I'm uncertain about the best way to generate these CIs.

## dicnucl-ins-dec.png  & dicnucl-del-dec.png
Randomization test was performed and yielded no significant results, even for p < 0.1. For this test, I sampled substrings of identical length from the v-loop sequences 100 times for every indel and compared observed values to this sampled distribution of nucleotide proportions.

## insertion / deletion timings
Histogram describing the distribution of insertion/deletion timings in days. I applied a correction factor to these values which adjusts for the loss in data as bin sizes get larger. Not sure how to display uncertainty / error measurements on this.
