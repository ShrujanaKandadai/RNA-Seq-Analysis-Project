library(plotly)

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
#look at pca.res in environment
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color=treatment, shape = genotype) +
  geom_point(size=4) +
  #geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

