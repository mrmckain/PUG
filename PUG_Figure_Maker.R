args<-commandArgs(TRUE)
library(ape)
library(RColorBrewer)

WGD_tree<-read.tree(file=args[1]);
WGD_counts<-read.table(file=args[2], header=F);

total_taxa<-length(WGD_tree$tip.label);
total_nodes<-length(WGD_tree$node.label);
cols<-rainbow(10000,start=1/6,end=.5);
total80<-sum(WGD_counts[,2]);
max80<-max(WGD_counts[,2]);

tenper80<-0.1*max80;

cols[1:(10000*tenper80/max80)]="black";
current_set80<-matrix(nrow=(total_taxa+total_nodes),ncol=1)
current_set80[1:total_nodes]=0;
 for (i in 1:length(WGD_counts[,2])){
   		current_set80[which(WGD_tree$node.label==WGD_counts[i,1])+total_taxa]<-WGD_counts[i,2];
 }
current_set80[is.na(current_set80)]<-0
xlims<-c(0,max80);
breaks<-0:10000/10000*(xlims[2]-xlims[1])+xlims[1];
whichColor<-function(p,cols,breaks){
      i<-1
     while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
      cols[i]
}
colors<-sapply(current_set80,whichColor,cols=cols,breaks=breaks);
WGD_colors<-colors[WGD_tree$edge[,2]]

add.color.bar<-function (position=c(0,1), leg, cols, title = NULL, lims = c(0, 1), digits = 1, prompt = TRUE, lwd = 4, outline = TRUE, ...) 
{
	y <- position[2]
    x <- position[1]
   fsize <- 1
    
    X <- x + cbind(0:(length(cols) - 1)/length(cols), 1:length(cols)/length(cols)) * 
        (leg)
    Y <- cbind(rep(y, length(cols)), rep(y, length(cols)))
    if (outline) 
        lines(c(X[1, 1], X[nrow(X), 2]), c(Y[1, 1], Y[nrow(Y), 
            2]), lwd = lwd + 2, lend = 2)
    for (i in 1:length(cols)) lines(X[i, ], Y[i, ], col = cols[i], 
        lwd = lwd, lend = 2)
    text(x = x, y = y, round(lims[1], digits), pos = 1, cex = .75)
    text(x = x + leg/2, y = y, round(lims[2], digits), pos = 1, 
        cex = .75)

    text(x = x + leg, y = y, round(lims[3], digits), pos = 1, 
        cex = .75)
    if (is.null(title)) 
        title <- "P(state=1)"
    text(x = (2 * x + leg)/2, y = y, title, pos = 3, cex = .9)
    }

pdf(file=args[3], colormodel="cmyk")
    plot.phylo(WGD_tree,edge.width=4,edge.color=WGD_colors,lend=2,new=FALSE,label.offset=.5);
add.color.bar(c(1,total_taxa*.9),6,cols,title="Unique Gene Duplications",lims=c(0,floor(max80/2),max80),1,prompt=FALSE, cex=0.75);


dev.off()
