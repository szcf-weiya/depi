#!/usr/bin/env Rscript
source('ped2textFn.r')
source('KinshipFN.r')
xped = read.table('../data/X.PED')
xtxt = ped2textFn(xped)
xk = kinshipFn(xtxt)
