#!/usr/bin/env Rscript
# Final S8 panels and Pearson statistics from an audited 293-individual source table.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: Rscript plot_s8_genotype_burden.R <S8_source_data.tsv> <output_dir>')
q <- args[[2]]
dir.create(q, recursive = TRUE, showWarnings = FALSE)
d <- read.delim(args[[1]])
required <- c('species','sample','FROH','target_allele_count','het_loci','hom_loci')
if (!all(required %in% names(d))) stop('Missing required S8 input columns')
counts <- table(d$species)
if (!identical(as.integer(counts[c('bamboo','teosinte','Aegilops')]), c(63L,90L,140L))) stop('Unexpected species sample counts')
if (any(d$target_allele_count != 2*d$hom_loci + d$het_loci)) stop('Inconsistent allele dosage')
sp <- c('bamboo','teosinte','Aegilops');metrics <- c('target_allele_count','het_loci','hom_loci')
res <- list();diags <- list()
for(s in sp){g<-subset(d,species==s);stopifnot(!anyDuplicated(g$sample),all(complete.cases(g)))
 for(k in metrics){ct<-cor.test(g$FROH,g[[k]],method='pearson');st<-cor.test(g$FROH,g[[k]],method='spearman',exact=FALSE);fit<-lm(g[[k]]~g$FROH);loo<-sapply(seq_len(nrow(g)),function(i)cor(g$FROH[-i],g[[k]][-i]));
 res[[length(res)+1]]<-data.frame(species=s,metric=k,n=nrow(g),r=unname(ct$estimate),p=ct$p.value,ci_low=ct$conf.int[1],ci_high=ct$conf.int[2],spearman_rho=unname(st$estimate),spearman_p=st$p.value,median=median(g[[k]]),min=min(g[[k]]),max=max(g[[k]]),loo_r_min=min(loo),loo_r_max=max(loo))
 diags[[length(diags)+1]]<-data.frame(species=s,metric=k,sample=g$sample,FROH=g$FROH,count=g[[k]],cooks_distance=cooks.distance(fit),leverage=hatvalues(fit),residual=residuals(fit))
 }}
r<-do.call(rbind,res);r$p_BH_9<-p.adjust(r$p,'BH');write.table(r,file.path(q,'correlation_results.tsv'),sep='\t',row.names=FALSE,quote=FALSE);write.table(do.call(rbind,diags),file.path(q,'influence_diagnostics.tsv'),sep='\t',row.names=FALSE,quote=FALSE);print(r,digits=10)
# Quantitative grid: genotype counts and their association with autozygosity.
# Unpolarized GERP >4, ALL regions, final S6 sample sets; no point removal or axis truncation.
cols<-c('#92B8B3','#509088','#214F53');names(cols)<-metrics
namesp<-c('Moso bamboo (n = 63)','Teosinte (n = 90)',"Tausch's goatgrass (n = 140)")
draw<-function(){
 layout(matrix(1:6,nrow=2,byrow=TRUE),heights=c(.95,1.15))
 par(mar=c(3.0,3.8,2.2,1.0),oma=c(1.1,.35,.2,.2),mgp=c(1.9,.55,0),tcl=-.25,family='Arial',cex=.68,cex.axis=.68)
 scales<-c(1e4,1e5,1e3)
 for(i in 1:3){
  g<-subset(d,species==sp[i]);hom<-g$hom_loci/scales[i];het<-g$het_loci/scales[i]
  xexpr<-list(expression('Number of loci (' %*% 10^4 * ')'),expression('Number of loci (' %*% 10^5 * ')'),expression('Number of loci (' %*% 10^3 * ')'))[[i]]
  if(i==1){
   # One common x-axis; the two genotype states occupy opposite halves of each bin.
   outer<-range(pretty(range(c(hom,het)),n=6));breaks<-seq(outer[1],outer[2],length.out=25)
   hh<-hist(hom,breaks=breaks,plot=FALSE);he<-hist(het,breaks=breaks,plot=FALSE)
   stopifnot(sum(hh$counts)==length(hom),sum(he$counts)==length(het))
   ymax<-ceiling((max(hh$counts,he$counts)+3)/5)*5
   xrange<-diff(range(breaks));plot(NA,xlim=c(min(breaks)-xrange*.02,max(breaks)+xrange*.02),ylim=c(-ymax*.035,ymax),xlab=xexpr,ylab='Individuals',main=namesp[i],bty='l',cex.main=1.05,xaxs='i',yaxs='i',yaxt='n')
   axis(2,at=seq(0,30,10),las=1,cex.axis=.68)
   width<-diff(breaks);gap<-width*.025
   rect(breaks[-length(breaks)]+gap,0,breaks[-1]-width/2-gap,hh$counts,col=cols['hom_loci'],border='gray35',lwd=.25)
   rect(breaks[-length(breaks)]+width/2+gap,0,breaks[-1]-gap,he$counts,col=cols['het_loci'],border='gray35',lwd=.25)
   abline(v=median(hom),col=cols['hom_loci'],lty=2,lwd=1.3);abline(v=median(het),col=cols['het_loci'],lty=2,lwd=1.3)
  }else{
   # Match the original broken-axis treatment, keeping each genotype's full range.
   first<-if(i==2)hom else het;second<-if(i==2)het else hom
   first_key<-if(i==2)'hom_loci' else 'het_loci';second_key<-if(i==2)'het_loci' else 'hom_loci'
   b1<-seq(min(first),max(first),length.out=25);b2<-seq(min(second),max(second),length.out=25)
   h1<-hist(first,breaks=b1,plot=FALSE);h2<-hist(second,breaks=b2,plot=FALSE)
   stopifnot(sum(h1$counts)==length(first),sum(h2$counts)==length(second))
   ytop<-ceiling((max(h1$counts,h2$counts)+3)/5)*5
   plot(NA,xlim=c(-.04,2.04),ylim=c(-ytop*.035,ytop),xaxt='n',yaxt='n',xlab=xexpr,ylab='',main=namesp[i],bty='l',cex.main=1.05,xaxs='i',yaxs='i')
   axis(2,at=seq(0,ytop,5),las=1,cex.axis=.68)
   tr1<-function(x)(x-min(b1))/diff(range(b1))*.96
   tr2<-function(x)1.04+(x-min(b2))/diff(range(b2))*.96
   e1<-tr1(b1);e2<-tr2(b2);gap1<-diff(e1)[1]*.025;gap2<-diff(e2)[1]*.025
   rect(head(e1,-1)+gap1,0,tail(e1,-1)-gap1,h1$counts,col=cols[first_key],border='gray35',lwd=.25)
   rect(head(e2,-1)+gap2,0,tail(e2,-1)-gap2,h2$counts,col=cols[second_key],border='gray35',lwd=.25)
   t1<-pretty(range(b1),n=3);t1<-t1[t1>=min(b1)&t1<=max(b1)]
   t2<-pretty(range(b2),n=3);t2<-t2[t2>=min(b2)&t2<=max(b2)]
   axis(1,at=c(tr1(t1),tr2(t2)),labels=c(t1,t2),cex.axis=.68)
   abline(v=tr1(median(first)),col=cols[first_key],lty=2,lwd=1.3)
   abline(v=tr2(median(second)),col=cols[second_key],lty=2,lwd=1.3)
   mtext('//',side=1,at=1,line=.15,cex=.72)
  }
  legend('topright',legend=c('Homozygous','Heterozygous'),fill=cols[c('hom_loci','het_loci')],bty='n',cex=.68)
  if(i==1)mtext('A',side=3,adj=0,line=.65,font=2,cex=1.05)
 }
 par(mar=c(3.2,3.8,2.0,1.0))
 for(i in 1:3){g<-subset(d,species==sp[i]);ymax<-max(g[,metrics])/1e5;xr<-range(g$FROH)
  plot(NA,xlim=xr,ylim=c(0,ymax*1.45),xlab=expression(italic(F)[ROH]),ylab=if(i==1)expression('Count ('%*%10^5*')')else '',bty='l')
  for(k in metrics){fit<-lm(g[[k]]/1e5~FROH,data=g);xx<-seq(xr[1],xr[2],length.out=150);pr<-predict(fit,newdata=data.frame(FROH=xx),interval='confidence');polygon(c(xx,rev(xx)),c(pr[,2],rev(pr[,3])),col=adjustcolor(cols[k],.12),border=NA);lines(xx,pr[,1],col=cols[k],lwd=1.2);points(g$FROH,g[[k]]/1e5,col=adjustcolor(cols[k],.72),pch=16,cex=.55)}
  rr<-r[r$species==sp[i],];lab<-paste0(c('Total','Heterozygous','Homozygous'),': r = ',sprintf('%.3f',rr$r),', P = ',formatC(rr$p,format='g',digits=3));legend('topright',legend=lab,text.col=cols,bty='n',cex=.62,x.intersp=.2,y.intersp=1.06);if(i==1)mtext('B',side=3,adj=0,line=.65,font=2,cex=1.05)
 }
 mtext('GERP >4 | Dashed lines: medians | Linear fits with 95% confidence intervals',side=1,outer=TRUE,line=.1,cex=.56)
}
png(file.path(q,'S8_corrected.png'),width=2400,height=1320,res=300,type='cairo',family='Arial');draw();dev.off()
cairo_pdf(file.path(q,'S8_corrected.pdf'),width=8,height=4.4,family='Arial');draw();dev.off()
png(file.path(q,'diagnostics.png'),width=2100,height=2100,res=220);par(mfrow=c(3,3),mar=c(4,4,2,1));for(s in sp){g<-subset(d,species==s);for(k in metrics){fit<-lm(g[[k]]~g$FROH);plot(fitted(fit),residuals(fit),main=paste(s,k),xlab='Fitted count',ylab='Residual',pch=16,cex=.6);abline(h=0,lty=2)}};dev.off()







