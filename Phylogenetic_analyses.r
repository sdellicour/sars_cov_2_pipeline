library(diagram)
library(HDInterval)
library(lubridate)
library(seraphim)
library(treeio)

# 1. Preparing the input files for the discrete phylogeographic analyses 
# 2. Analysing the outputs of the preliminary discrete phylogeographic analysis 
# 3. Preparing the continuous phylogeographic analysis (RRW, Cauchy model)
# 4. Running BEAST and building the maximum clade consensus (MCC) tree
# 5. Extracting spatio-temporal information embedded in MCC and posterior trees
# 6. Generating a dispersal history graph (mapped MCC trees, 80% HPD polygons)
# 7. Jackknife on phylogenetic branches (randomly keeping 75% of branches)
# 8. Comparing intra- and inter-province/municipality lineage migration events
# 9. Estimating and plotting dispersal statistics associated with lineages

analysis = "TreeTime_100620"; removeSuspiciousHomoplasies = FALSE
	
data1 = read.csv("Sequences_metadata/SARS-CoV-2_KULeuven_100620.csv", sep=";")
data2 = read.csv("Sequences_metadata/SARS-CoV-2_ULiegeSeq_020620.csv", sep=";")
writingFiles = FALSE; showingPlots = FALSE

# 1. Preparing the input files for the discrete phylogeographic analyses 

tree = read.tree(paste0(analysis,".tre"))
if (grepl("TreeTime",analysis) == TRUE)
	{
		seqIDs = tree$tip.label; countries = rep(NA, length(seqIDs)); collectionDates = rep(NA, length(seqIDs))
		for (i in 1:length(seqIDs))
			{
				if (grepl("hCoV-19",seqIDs[i]))
					{
						countries[i] = unlist(strsplit(seqIDs[i],"\\/"))[2]
					}	else	{
						countries[i] = unlist(strsplit(seqIDs[i],"\\/"))[1]
					}
				if (length(unlist(strsplit(seqIDs[i],"\\|"))) == 3)
					{
						collectionDates[i] = unlist(strsplit(seqIDs[i],"\\|"))[length(unlist(strsplit(seqIDs[i],"\\|")))]
					}
				if (length(unlist(strsplit(seqIDs[i],"\\|"))) == 4)
					{
						collectionDates[i] = unlist(strsplit(seqIDs[i],"\\|"))[length(unlist(strsplit(seqIDs[i],"\\|")))-1]
					}
			}
		tab = cbind(seqIDs,countries,collectionDates); colnames(tab) = c("Strain","Country","Collection Data")
		write.csv(tab, paste0(analysis,".csv"), row.names=F, quote=F)
		data = read.csv(paste0(analysis,".csv"))
	}
numberOfDifferentCountries = length(unique(data[,"Country"]))
if (showingPlots)
	{
		dev.new(width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, col="gray30", edge.color="gray30")
		for (i in 1:dim(tree$edge)[1])
			{
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("Belgium",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
					}
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("ULG-",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
	}
if (removeSuspiciousHomoplasies)
	{
		suspicious_homoplasies = read.csv("Suspicious_homopl.csv")[,"sequence"]
		suspicious_homoplasies = as.character(unique(suspicious_homoplasies[!is.na(suspicious_homoplasies)]))
		labs = unique(data[which(data[,"Country"]=="Belgium"),"Submitting.Lab"])
		sedIDs = gsub("Belgium\\/","",gsub("\\/2020","",data[which(data[,"Country"]=="Belgium"),"Strain"]))
		sedIDs_KUL = sedIDs[(!grepl("ULG-",sedIDs))&(!grepl("UGent-",sedIDs))]
		temp = data1[which(data1[,"GisAID"]=="OK"),]
		temp[,"sequence.name"] = gsub("SARS2-CoV\\/Belgium\\/Human\\/","",gsub("\\/2020","",temp[,"sequence.name"]))
		temp[,"sequence.name"] = gsub("SARS2-Cov\\/Belgium\\/Human\\/","",gsub("\\/2020","",temp[,"sequence.name"]))
		indices = which(!temp[,"sequence.name"]%in%sedIDs_KUL)
		discardedKULsequences = temp[indices,"sequence.name"]
		if (showingPlots)
			{
				dev.new(width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
				plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, col="gray30", edge.color="gray30")
				for (i in 1:dim(tree$edge)[1])
					{
						if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("Belgium",tree$tip.label[tree$edge[i,2]])))
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}
						if ((!tree$edge[i,2]%in%tree$edge[,1]) & (tree$tip.label[tree$edge[i,2]]%in%suspicious_homoplasies))
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="red")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}
					}
				add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
			}
		tree = drop.tip(tree, suspicious_homoplasies)
		if (writingFiles) write.nexus(tree, file="Newick_tree_for_XML_file.tree")
	}
txt = c(); tab = c()
for (i in 1:length(tree$tip.label))
	{
		index = which(data[,1]==tree$tip.label[i])
		date = as.character(data[index,"Collection.Data"])
		if (date != "")
			{
				txt = c(txt, paste0(">",tree$tip.label[i]),"NNNN")
				if (grepl("Nextstrain",analysis))
					{
						location = unlist(strsplit(tree$tip.label[i],"\\/"))[1]
					}
				if (grepl("TreeTime",analysis))
					{
						if (grepl("hCoV-19",tree$tip.label[i]))
							{
								location = unlist(strsplit(tree$tip.label[i],"\\/"))[2]
							}	else	{
								location = unlist(strsplit(tree$tip.label[i],"\\/"))[1]
							}
					}
				if (location != "Belgium") location = "other"
				tab = rbind(tab, cbind(tree$tip.label[i],location,date))
			}
	}
colnames(tab) = c("trait","location","collection_date")
if (writingFiles) write.table(tab, paste0(analysis,".txt"), row.names=F, quote=F, sep="\t")
if (writingFiles) write(txt, paste0(analysis,".fasta"))

# 2. Analysing the outputs of the preliminary discrete phylogeographic analysis 

burnIn = 101; computingHPDInterval = FALSE # N.B.: long analysis
if (computingHPDInterval)
	{
		trees = scan(paste0(analysis,".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		indices1 = which(!grepl("tree STATE_",trees)); indices2 = which(grepl("tree STATE_",trees))
		belgianBranches_list = rep(NA,length(trees))
		belgianIntroductions_list = rep(NA,length(trees))
		belgianTipBranches_list = rep(NA,length(trees))
		for (i in (burnIn+1):length(indices2))
			{
				tree1 = trees[c(indices1[1:(length(indices1)-1)],indices2[i],indices1[length(indices1)])]
				write(tree1, paste0("TreeTime_100620_sampled_tree_",i,".tree"))
				tree2 = readAnnotatedNexus(paste0("TreeTime_100620_sampled_tree_",i,".tree"))
				belgianBranches = 0; belgianIntroductions = 0; belgianTipBranches = 0
				for (j in 1:dim(tree2$edge)[1])
					{
						if (tree2$annotations[[j]]$location == "Belgium")
							{
								belgianBranches = belgianBranches + 1
								index = which(tree2$edge[,2]==tree2$edge[j,1])
								if (tree2$annotations[[index]]$location != "Belgium")
									{
										belgianIntroductions = belgianIntroductions + 1
									}
								if (!tree2$edge[j,2]%in%tree2$edge[,1])
									{
										belgianTipBranches = belgianTipBranches + 1
									}
							}
					}
				belgianBranches_list[i] = belgianBranches
				belgianIntroductions_list[i] = belgianIntroductions
				belgianTipBranches_list[i] = belgianTipBranches
				file.remove(paste0("TreeTime_100620_sampled_tree_",i,".tree"))
			}
		quantiles = quantile(belgianIntroductions_list[!is.na(belgianIntroductions_list)],probs=c(0.025,0.975))
		HPD = HDInterval::hdi(belgianIntroductions_list[!is.na(belgianIntroductions_list)])[1:2]
		cat("A minimum number of ",median(belgianIntroductions_list[!is.na(belgianIntroductions_list)])," lineage introductions (95% HPD interval = [",
			HPD[1],"-",HPD[2],"])"," identified from the global phylogenetic analysis of ",belgianTipBranches," SARS-CoV-2 sampled in Belgium (20-04-2020)",sep="")
		# Results for the 1° analysis based on the Nextstrain tree of the 20-04-20 (without removing tip branches associated with suspicious homoplasie):
			# a minimum number of 166 lineage introductions (95% HPD interval = [161-171]) identified from the global phylogenetic analysis of 391 SARS-CoV-2 sampled in Belgium
		# Results for the 2° analysis based on the Nextstrain tree of the 20-04-20 (when dropping tip branches associated with suspicious homoplasie):
			# a minimum number of 157 lineage introductions (95% HPD interval = [151-162]) identified from the global phylogenetic analysis of 370 SARS-CoV-2 sampled in Belgium
		# Results for the 3° analysis based on the Nextstrain tree of the 20-04-20 (based on a lignment without sequences associated with suspicious homoplasie):
			# a minimum number of 144 lineage introductions (95% HPD interval = [138-150]) identified from the global phylogenetic analysis of 370 SARS-CoV-2 sampled in Belgium
		# Results for the 7° analysis based on the Nextstrain tree of the 10-06-20 (based on a lignment without sequences associated with suspicious homoplasie):
			# a minimum number of 331 lineage introductions (95% HPD interval = [315-344]) identified from the global phylogenetic analysis of 740 SARS-CoV-2 sampled in Belgium
	}
tree = readAnnotatedNexus("TreeTime_100620.tree")
if (showingPlots)
	{
		samplingDates = decimal_date(ymd(gsub("\\/","-",tab[,"collection_date"]))); mostRecentSamplingYear = max(samplingDates)
		selectedDates = decimal_date(ymd(c("2019-11-01","2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01")))
		rootHeight = max(nodeHeights(tree)); root_time = mostRecentSamplingYear-rootHeight
		selectedLabels = c("01-11-2019","01-12-2019","01-01-2020","01-02-2020","01-03-2020","01-04-2020","01-05-2020","01-06-2020")
		cols = rep("gray30",dim(tree$edge)[1]); lwds = rep(0.1,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Belgium") & (tree$annotations[[i]]$location=="Belgium"))
							{
								cols[i] = "chartreuse3"; lwds[i] = 0.4
							}
					}
			}
		pdf("Figure_1_NEW.pdf", width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		# dev.new(width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("Belgium",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.3, col="chartreuse3")
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.3, col="gray30", lwd=0.5)
					}
				if (tree$annotations[[i]]$location == "Belgium")
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if (tree$annotations[[index]]$location != "Belgium")
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
		cols = rep("gray50",dim(tree$edge)[1]); lwds = rep(0.05,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Belgium") & (tree$annotations[[i]]$location=="Belgium"))
							{
								cols[i] = "chartreuse3"; lwds[i] = 0.4
							}
					}
			}
		dev.off()
		pdf("Figure_S1_NEW.pdf", width=11, height=8) # dev.new(width=11, height=8)
		plot(tree, show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$annotations[[i]]$location == "Belgium")
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if (tree$annotations[[index]]$location != "Belgium")
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}	else	{
								if (!tree$edge[i,2]%in%tree$edge[,1])
									{	}
							}
					}
			}
		axis(lwd=0.2, at=selectedDates-root_time, labels=selectedLabels, cex.axis=0.65, mgp=c(0,0.1,-0.9), lwd.tick=0.2, 
			 col.lab="gray30", col="gray30", tck=-0.01, side=1)
		dev.off()
	}
communes = shapefile("Shapefile_communes/Shapefile_post_codes.shp")
provinces = spTransform(raster::getData("GADM", country="BEL", level=2), crs(communes))
belgium = spTransform(raster::getData("GADM", country="BEL", level=0), crs(communes))
pop = projectRaster(raster("WorldPop_pop_raster.tif"), crs=crs(communes)); pop[] = log(pop[])
belgianBranches = c(); belgianIntroductions = c()
belgianTipBranches = c(); sampledSequences = c()
for (i in 1:dim(tree$edge)[1])
	{
		if (tree$annotations[[i]]$location == "Belgium")
			{
				belgianBranches = c(belgianBranches,i)
				index = which(tree$edge[,2]==tree$edge[i,1])
				if (tree$annotations[[index]]$location != "Belgium")
					{
						belgianIntroductions = c(belgianIntroductions, i)
					}
				if (!tree$edge[i,2]%in%tree$edge[,1])
					{
						belgianTipBranches = c(belgianTipBranches, i)
						sampledSequences = c(sampledSequences, tree$tip.label[tree$edge[i,2]])
					}
			}
	}
for (i in 1:length(belgianIntroductions))
	{
		if (i == 1) clusters1 = list()
		if (tree$edge[belgianIntroductions[i],2]%in%tree$edge[,1])
			{
				subtree = tree_subset(tree, tree$edge[belgianIntroductions[i],2], levels_back=0)
				clusters1[[i]] = gsub("'","",subtree$tip.label)
			}	else		{
				clusters1[[i]] = gsub("'","",tree$tip.label[tree$edge[belgianIntroductions[i],2]])
			}
	}
sampledSequences = gsub("'","",sampledSequences)
if (!file.exists(paste0("Sampling_Belgium.csv")))
	{
		samplingData = matrix(nrow=length(sampledSequences), ncol=5)
		colnames(samplingData) = c("sequenceID","collectionDate","postCode","longitude","latitude")
		samplingData[,"sequenceID"] = sampledSequences
		for (i in 1:dim(samplingData)[1])
			{
				index = which(data[,"Strain"]==samplingData[i,"sequenceID"])
				date = ymd(gsub("\\/","-",data[index,"Collection.Data"]))
				samplingData[i,"collectionDate"] = decimal_date(date)
				ID = unlist(strsplit(samplingData[i,"sequenceID"],"\\/"))[3]
				index1 = which(grepl(gsub("Rega","rega",ID),data1[,"sequence.name"]))
				if (length(index1) == 1)
					{
						samplingData[i,"postCode"] = data1[index1,"ZIP"]
						indices = which(communes@data[,"nouveau_PO"]==data1[index1,"ZIP"])
						if (length(indices) > 0)
							{
								maxArea = 0; polIndex1 = 0; polIndex2 = 0
								for (j in 1:length(indices))
									{
										for (k in 1:length(communes@polygons[[indices[j]]]@Polygons))
											{
												if (maxArea < communes@polygons[[indices[j]]]@Polygons[[k]]@area)
													{
														maxArea = communes@polygons[[indices[j]]]@Polygons[[k]]@area; polIndex1 = indices[j]; polIndex2 = k
													}
											}
									}
								pol = communes@polygons[[polIndex1]]@Polygons[[polIndex2]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = communes@proj4string
								samplingData[i,c("longitude","latitude")] = coordinates(pol) # to avoid a jitter:
								samplingData[i,c("longitude","latitude")] = spsample(pol, 1, type="random")@coords
							}
					}
				index1 = which(grepl(ID,data2[,"gisaid.virus.name"]))
				if (length(index1) == 1)
					{
						samplingData[i,"postCode"] = as.character(data2[index1,"Postal.code"])
						indices = which(communes@data[,"nouveau_PO"]==data2[index1,"Postal.code"])
						if (length(indices) > 0)
							{
								maxArea = 0; polIndex1 = 0; polIndex2 = 0
								for (j in 1:length(indices))
									{
										for (k in 1:length(communes@polygons[[indices[j]]]@Polygons))
											{
												if (maxArea < communes@polygons[[indices[j]]]@Polygons[[k]]@area)
													{
														maxArea = communes@polygons[[indices[j]]]@Polygons[[k]]@area; polIndex1 = indices[j]; polIndex2 = k
													}
											}
									}
								pol = communes@polygons[[polIndex1]]@Polygons[[polIndex2]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = communes@proj4string
								samplingData[i,c("longitude","latitude")] = coordinates(pol) # to avoid a jitter:
								samplingData[i,c("longitude","latitude")] = spsample(pol, 1, type="random")@coords
							}
					}
			}
		print(samplingData[which(is.na(samplingData[,"postCode"])),"sequenceID"])
		write.csv(samplingData, "Sampling_Belgium.csv", quote=F, row.names=F)
	}	
samplingData = read.csv("Sampling_Belgium.csv", head=T)
for (i in 1:length(belgianIntroductions))
	{
		tab = c()
		if (i == 1)
			{
				clusters2 = list(); centroids = list()
			}
		for (j in 1:length(clusters1[[i]]))
			{
				index = which(samplingData[,"sequenceID"]==clusters1[[i]][j])
				if (length(index) == 1)
					{
						line = cbind(as.numeric(samplingData[index,"collectionDate"]),as.numeric(samplingData[index,"longitude"]),as.numeric(samplingData[index,"latitude"]))
						row.names(line) = clusters1[[i]][j]; tab = rbind(tab, line)
					}
			}
		colnames(tab) = c("collectionDate","longitude","latitude"); clusters2[[i]] = tab
		centroids[[i]] = cbind(mean(tab[!is.na(tab[,"longitude"]),"longitude"]), mean(tab[!is.na(tab[,"latitude"]),"latitude"]))
	}
clusterSizes = rep(NA, length(clusters1))
collectionDates = c()
for (i in 1:length(clusters1))
	{
		clusterSizes[i] = length(clusters1[[i]])
		collectionDates = c(collectionDates, clusters2[[i]][,"collectionDate"])
	}
if (showingPlots)
	{
		collectionDates_filetered = collectionDates
		dev.new(width=3.3, height=8); par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(2,2,1,1), lwd=0.2, col="gray30")
		hist(clusterSizes, breaks=50, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		hist(collectionDates_filetered, breaks=65, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=decimal_date(ymd(c("2020-02-01","2020-03-01","2020-04-01","2020-05-01"))),
			 labels=c("01-02-2020","01-03-2020","01-04-2020","01-05-2020"))
	}
if (showingPlots)
	{
		dev.new(width=7.3, height=6); par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
		cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
		plot(pop, col=cols, axes=F, ann=F, box=F, legend=F)
		plot(provinces, border="white", col=NA, add=T, lwd=1.0)
		plot(belgium, border="gray30", col=NA, add=T, lwd=0.4)
		for (i in 1:length(clusters1))
			{
				if (!is.na(centroids[[i]][,1]))
					{
						if (length(clusters1[[i]]) > 1)
							{
								for (j in 1:dim(clusters2[[i]])[1])
									{
										if (!is.na(clusters2[[i]][j,1]))
											{
												segments(centroids[[i]][,1],centroids[[i]][,2],clusters2[[i]][j,"longitude"],clusters2[[i]][j,"latitude"], lwd=0.5, col="gray30")	
											}
									}
							}
					}
			}
		for (i in 1:length(clusters1))
			{
				if (!is.na(centroids[[i]][,1]))
					{
						for (j in 1:dim(clusters2[[i]])[1])
							{
								points(clusters2[[i]][,"longitude"], clusters2[[i]][,"latitude"], pch=16, cex=0.8, col="chartreuse3")
								points(clusters2[[i]][,"longitude"], clusters2[[i]][,"latitude"], pch=1, cex=0.8, col="gray30", lwd=0.2)
							}
					}
			}
		for (i in 1:length(clusters1))
			{
				if (!is.na(centroids[[i]][,1]))
					{
						if (length(clusters1[[i]]) > 1)
							{
								points(centroids[[i]][,1], centroids[[i]][,2], pch=16, cex=0.6, col="red")
								points(centroids[[i]][,1], centroids[[i]][,2], pch=1, cex=0.6, col="gray30", lwd=0.2)
							}
					}
			}
		legendRast = raster(as.matrix(c(min(pop[],na.rm=T),max(pop[],na.rm=T))))
		mtext("Human population (log-transformed)", col="gray30", cex=0.7, line=-23, at=601000)
		plot(legendRast, legend.only=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.141,0.409,0.18,0.19),
			 alpha=1, horizontal=T, legend.args=list(text="", cex=0.7, line=0.5, col="gray30"), axis.args=list(cex.axis=0.55, lwd=0,
			 lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,-0.05,0)))
	}

# 3. Preparing the continuous phylogeographic analysis (RRW, Cauchy model)

analyses = c()
template = scan("RRW_XMLtemplate.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
for (i in 1:length(clusters2))
	{
		if ((dim(clusters2[[i]])[1] >= 3)&(sum(!is.na(clusters2[[i]][,"longitude"])) >= 3))
			{
				analyses = c(analyses, paste0("Clade_",i)); cluster = clusters2[[i]][which(!is.na(clusters2[[i]][,"longitude"])),]
				xml = gsub("TEMPLATE", paste0("Clade_",i), template)
				tre = tree_subset(tree, tree$edge[belgianIntroductions[i],2], levels_back=0)
				tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%row.names(cluster))]
				if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
				write.tree(tre, paste0("Phylogeographic_runs/Clade_",i,".tre"))
				tre = scan(paste0("Phylogeographic_runs/Clade_",i,".tre"), what="", sep="\n", quiet=T)
				sink(file=paste0("Phylogeographic_runs/Clade_",i,".xml"))
				for (j in 1:length(xml))
					{
						cat(xml[j]); cat("\n")
						if (xml[j]=="\t<taxa id=\"taxa\">")
							{
								for (k in 1:dim(cluster)[1])
									{
										cat(paste0("\t\t<taxon id=\"",row.names(cluster)[k],"\">","\n"))
										cat(paste0("\t\t\t<date value=\"",cluster[k,"collectionDate"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
										cat("\t\t\t<attr name=\"latitude\">\n")
										cat(paste0("\t\t\t\t",cluster[k,"latitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t\t<attr name=\"longitude\">\n")
										cat(paste0("\t\t\t\t",cluster[k,"longitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t\t<attr name=\"coordinates\">\n")
										cat(paste0("\t\t\t\t",cluster[k,"latitude"]," ",cluster[k,"longitude"],"\n"))
										cat("\t\t\t</attr>\n")
										cat("\t\t</taxon>\n")
									}
							}
						if (xml[j]=="\t<alignment id=\"alignment\" dataType=\"nucleotide\">")
							{
								for (k in 1:dim(cluster)[1])
									{
										cat("\t\t<sequence>\n")
										cat(paste0("\t\t\t<taxon idref=\"",row.names(cluster)[k],"\"/>","\n"))
										cat("\t\t\tNNNN\n")
										cat("\t\t</sequence>\n")
									}
							}
						if (xml[j]=="\t<newick id=\"startingTree\">")
							{
								cat(paste0("\t\t",tre,"\n"))
							}
					}
				sink(NULL)
			}
	}

# 4. Running BEAST and building the maximum clade consensus (MCC) tree

source("MCC_tree_extraction.r")
sink(file=paste0("Phylogeographic_runs/Analyses.sh"))
for (i in 1:length(analyses))
	{
		cat(paste0("java -jar beast_1104.jar -overwrite ",analyses[i],".xml\n"))
	}
sink(NULL)
runningNewAnalyses = FALSE
wd = getwd(); setwd(paste0(wd,"/Phylogeographic_runs/"))
if (runningNewAnalyses)
	{
		system("bash Analyses.sh", ignore.stdout=T, ignore.stderr=F)
		burnIns = rep(1001, length(analyses)); burnIns[which(analyses%in%paste0("Clade_",c(30,38,129,157)))] = 5001
		for (i in 1:length(analyses))
			{
system(paste0("BEAST_1104/bin/treeannotator -burninTrees ",burnIns[i]," -heights keep ",analyses[i],".trees ",analyses[i],".tree"), ignore.stdout=F, ignore.stderr=F)
			}
	}
setwd(wd)

# 5. Extracting spatio-temporal information embedded in MCC and posterior trees

prov_WGS = raster::getData("GADM", country="BEL", level=2)
polIndex1 = which(prov_WGS@data[,"NAME_2"]=="Liège")
maxArea = 0; polIndex2 = 0
for (i in 1:length(prov_WGS@polygons[[polIndex1]]@Polygons))
	{
		if (maxArea < prov_WGS@polygons[[polIndex1]]@Polygons[[i]]@area)
			{
				maxArea = prov_WGS@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
			}
	}
pol = prov_WGS@polygons[[polIndex1]]@Polygons[[polIndex2]]
p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
pol = sps; proj4string(pol) = crs(raster("WorldPop_pop_raster.tif"))
pop_liege_WGS = raster::mask(crop(raster("WorldPop_pop_raster.tif"),pol),pol)
wd = getwd(); setwd(paste0(wd,"/Phylogeographic_runs/"))
for (i in 1:length(analyses))
	{
		if (!file.exists(paste0(analyses[i],".csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2])
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
				mcc_tre = readAnnotatedNexus(paste0(analyses[i],".tree"))
				mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
				write.csv(mcc_tab, paste0(analyses[i],".csv"), row.names=F, quote=F)
			}
	}
nberOfTreesToSample = 1000; randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 5
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0(analyses[i],"_ext")
		if (!file.exists(paste0(localTreesDirectory,"/TreeExtractions_1.csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2]); burnIn = burnIns[i]
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collectionDate"])
				allTrees = scan(file=paste0(analyses[i],".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
			}
	}
for (i in 1:length(analyses))
	{
		tab = read.csv(paste0(analyses[i],".csv"), head=T)
		if (i == 1)
			{
				all = tab
			}	else	{
				maxNodeID = max(all[,c("node1","node2")])
				tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
				all = rbind(all, tab)
			}
	}
write.csv(all, "All_clades.csv", row.names=F, quote=F)
dir.create(file.path("All_clades_ext1"), showWarnings=F)
dir.create(file.path("All_clades_ext2"), showWarnings=F)
dir.create(file.path("Bf_180320_ext2"), showWarnings=F)
dir.create(file.path("Af_180320_ext2"), showWarnings=F)
dir.create(file.path("Bf_180320_Liege"), showWarnings=F)
dir.create(file.path("Af_180320_Liege"), showWarnings=F)
nberOfExtractionFiles = nberOfTreesToSample
for (i in 1:nberOfExtractionFiles)
	{
		for (j in 1:length(analyses))
			{
				tab = read.csv(paste0(analyses[j],"_ext/TreeExtractions_",i,".csv"), head=T)
				if (j == 1)
					{
						all = tab
					}	else	{
						maxNodeID = max(all[,c("node1","node2")])
						tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
						all = rbind(all, tab)
					}
			}
		write.csv(all, paste0("All_clades_ext1/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		temp1 = all[,c("startLon","startLat")]; temp2 = all[,c("endLon","endLat")]
		coordinates(temp1) = ~ startLon + startLat; crs(temp1) = crs(communes)
		coordinates(temp2) = ~ endLon + endLat; crs(temp2) = crs(communes)
		temp1 = spTransform(temp1, CRS("+init=epsg:4326"))@coords
		temp2 = spTransform(temp2, CRS("+init=epsg:4326"))@coords
		all[,c("startLon","startLat")] = temp1; all[,c("endLon","endLat")] = temp2
		write.csv(all, paste0("All_clades_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		bf1 = all[which(all[,"endYear"]<decimal_date(dmy("18-03-2020"))),]
		af1 = all[which(all[,"endYear"]>=decimal_date(dmy("18-03-2020"))),]
		write.csv(bf1, paste0("Bf_180320_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		write.csv(af1, paste0("Af_180320_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		vS1 = raster::extract(pop_liege_WGS, bf1[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege_WGS, bf1[,c("endLon","endLat")])
		bf2 = bf1[which((!is.na(vS1))&(!is.na(vS2))),]
		vS1 = raster::extract(pop_liege_WGS, af1[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege_WGS, af1[,c("endLon","endLat")])
		af2 = af1[which((!is.na(vS1))&(!is.na(vS2))),]
		write.csv(bf2, paste0("Bf_180320_Liege/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		write.csv(af2, paste0("Af_180320_Liege/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}
setwd(wd)

# 6. Generating a dispersal history graph (mapped MCC trees, 80% HPD polygons)

localTreesDirectory = paste0("Phylogeographic_runs/All_clades_ext1"); nberOfExtractionFiles = 1000
percentage = 80; prob = percentage/100; precision = 1/(365/7); croppingPolygons = TRUE
mcc = read.csv("Phylogeographic_runs/All_clades.csv", head=T); startDatum = min(mcc[,"startYear"])
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
if (showingPlots)
	{
		colourScale = rev(colorRampPalette(brewer.pal(11,"BrBG"))(141)[16:116])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		startYears_indices = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours = colourScale[startYears_indices]
		endYears_colours = colourScale[endYears_indices]
		polygons_colours = rep(NA, length(polygons))
		cexNode = 0.8; LWD = 1.0
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours[i] = paste0(colourScale[polygon_index],"40")
			}
		firstTimePeriod = TRUE; secondTimePeriod = FALSE
		firstTimePeriod = FALSE; secondTimePeriod = TRUE
		firstTimePeriod = FALSE; secondTimePeriod = FALSE
		if ((firstTimePeriod == TRUE)|(secondTimePeriod == TRUE)) { cexNode = 1.1; LWD = 2.0 }
		pdf("Figure_2_NEW.pdf", width=7.3, height=6) # dev.new(width=7.3, height=6)
		par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
		cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
		cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
		if ((firstTimePeriod == TRUE)|(secondTimePeriod == TRUE))
			{
				plot(pop, col=cols, axes=F, ann=F, box=F, legend=F)
				plot(provinces, border="white", col=NA, add=T, lwd=LWD)
			}	else		{
				plot(provinces, border=NA, col="gray95", lwd=LWD)
			}
		if ((firstTimePeriod == FALSE)&(secondTimePeriod == FALSE))
			{
				for (i in 1:length(polygons))
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(belgium)
						if (croppingPolygons == TRUE) pol = crop(pol, belgium)
						plot(pol, axes=F, col=polygons_colours[i], add=T, border=NA)
					}
			}
		plot(provinces, border="white", col=NA, add=T, lwd=LWD)
		plot(belgium, border="gray30", col=NA, add=T, lwd=0.4); croppingPolygons = TRUE
		selectedBranches = 1:dim(mcc)[1]
		if (firstTimePeriod) selectedBranches = which(mcc[,"endYear"]<decimal_date(ymd("2020-03-18")))
		if (secondTimePeriod) selectedBranches = which(mcc[,"startYear"]>decimal_date(ymd("2020-03-18")))
		for (i in selectedBranches)
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						    arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		if ((firstTimePeriod == FALSE)&(secondTimePeriod == FALSE))
			{
				for (i in 1:length(clusters2))
					{
						if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
							{
								if (sum(!is.na(clusters2[[i]][,"longitude"])) == 2)
									{
										indices = which(!is.na(clusters2[[i]][,"longitude"]))
										if (firstTimePeriod)
											{
												indices = which((!is.na(clusters2[[i]][,"longitude"]))&(clusters2[[i]][,"collectionDate"]<decimal_date(ymd("2020-03-18"))))
											}
										if (secondTimePeriod)
											{
												indices = which((!is.na(clusters2[[i]][,"longitude"]))&(clusters2[[i]][,"collectionDate"]>decimal_date(ymd("2020-03-18"))))
											}
										if (length(indices) == 2)
											{
												curvedarrow(cbind(clusters2[[i]][indices[1],"longitude"],clusters2[[i]][indices[1],"latitude"]),
															cbind(clusters2[[i]][indices[2],"longitude"],clusters2[[i]][indices[2],"latitude"]),
															arr.length=0, arr.width=0, lwd=0.2, lty=2, lcol="gray10", arr.col=NA, arr.pos=F,
															curve=0.1, dr=NA, endhead=F)
											}
									}
							}
					}
			}
		for (i in 1:length(clusters2))
			{
				if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
					{
						for (j in 1:dim(clusters2[[i]])[1])
							{
								if (!is.na(clusters2[[i]][j,"longitude"]))
									{
										plotTheNode = TRUE
										if ((firstTimePeriod==TRUE)&(clusters2[[i]][j,"collectionDate"]>decimal_date(ymd("2020-03-18")))) plotTheNode = FALSE
										if ((secondTimePeriod==TRUE)&(clusters2[[i]][j,"collectionDate"]<decimal_date(ymd("2020-03-18")))) plotTheNode = FALSE
										if (plotTheNode)
											{
												index = (((clusters2[[i]][j,"collectionDate"]-minYear)/(maxYear-minYear))*100)+1
												points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=16, col=colourScale[index], cex=cexNode)
												points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
											}
									}
							}
					}
			}
		for (i in rev(selectedBranches))
			{
				if (!mcc[i,"node1"]%in%mcc[selectedBranches,"node2"])
					{
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
					}
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
			}
		selectedDates = decimal_date(ymd(c("2020-02-03","2020-03-03","2020-03-18","2020-04-03","2020-05-03")))
		selectedLabels = c("03-02-2020","03-03-2020","18-03-2020","03-04-2020","03-05-2020")
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.1,0.5,0.100,0.112),
			 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.6, col.axis="gray30", line=0, mgp=c(0,0.00,0),
		     at=selectedDates, labels=selectedLabels))
		dev.off()
	}
if (showingPlots)
	{
		polIndex1 = which(provinces@data[,"NAME_2"]=="Liège")
		maxArea = 0; polIndex2 = 0; cexNode = 0.85
		for (i in 1:length(provinces@polygons[[polIndex1]]@Polygons))
			{
				if (maxArea < provinces@polygons[[polIndex1]]@Polygons[[i]]@area)
					{
						maxArea = provinces@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
					}
			}
		pol = provinces@polygons[[polIndex1]]@Polygons[[polIndex2]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sps; proj4string(pol) = communes@proj4string
		pop_liege = raster::mask(crop(pop,pol),pol)
		vS1 = raster::extract(pop_liege, mcc[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege, mcc[,c("endLon","endLat")])
		sub = mcc[which((!is.na(vS1))&(!is.na(vS2))),]		
		colourScale = rev(colorRampPalette(brewer.pal(11,"BrBG"))(141)[16:116])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		startYears_indices = (((sub[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices = (((sub[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours = colourScale[startYears_indices]
		endYears_colours = colourScale[endYears_indices]
		pdf("Figure_S2_NEW.pdf", width=11, height=4) # dev.new(width=11, height=4)
		par(mfrow=c(1,2), oma=c(0,2,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		for (i in 1:2)
			{
				cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
				cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
				plot(pop_liege, col=cols, axes=F, ann=F, box=F, legend=F)
				plot(pol, border="gray30", col=NA, add=T, lwd=0.4)
				if (i == 1) selectedBranches = which(sub[,"endYear"]<decimal_date(dmy("18-03-2020")))
				if (i == 2) selectedBranches = which(sub[,"startYear"]>=decimal_date(dmy("18-03-2020")))
				for (j in selectedBranches)
					{
						curvedarrow(cbind(sub[j,"startLon"],sub[j,"startLat"]), cbind(sub[j,"endLon"],sub[j,"endLat"]), arr.length=0,
						    		arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
				for (j in rev(selectedBranches))
					{
						if (!sub[j,"node1"]%in%sub[selectedBranches,"node2"])
							{
								points(sub[j,"startLon"], sub[j,"startLat"], pch=16, col=startYears_colours[j], cex=cexNode)
								points(sub[j,"startLon"], sub[j,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
							}
						points(sub[j,"endLon"], sub[j,"endLat"], pch=16, col=endYears_colours[j], cex=cexNode)
						points(sub[j,"endLon"], sub[j,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
					}
			}
		dev.off()
	}

