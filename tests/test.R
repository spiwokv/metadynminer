library('metadynminer')
library('covr')

tfes<-fes(acealanme, imax=5000)
expect_equal(sum(tfes$fes), -2024857)

