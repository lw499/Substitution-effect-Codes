truer = c(2.347211, -1.634608)

tmpdat = read.csv("seq_reg.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean
gformseq2 = tmpdat;
gformseqSE = apply(gformseq2, 2, var)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = (gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, bpgformseq*100, MSE*100)

tmpdat = read.csv("ipw.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean
gformseq2 = tmpdat;
gformseqSE = apply(gformseq2, 2, var)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = (gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, bpgformseq*100, MSE*100)


tmpdat = read.csv("dr.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean
gformseq2 = tmpdat;
gformseqSE = apply(gformseq2, 2, var)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = (gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, bpgformseq*100, MSE*100)
