load("~/Project/epistasis/epi/wholeData_Li_JUN.RData")
write.table(Data, 'wholeData.txt', sep = " ")
write.table(COV, 'cov_whole.txt', sep = " ")
