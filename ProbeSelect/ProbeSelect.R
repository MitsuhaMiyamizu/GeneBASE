#released 09/26/2007
library(affy)
major_group<-function(grouping) {
	group_id=names(table(grouping))[table(grouping)==max(table(grouping))][1]
	return(names(grouping)[grouping==group_id])

}

split_probe_name<-function(name) {
	return(unlist(strsplit(name,"\\|")))
}

condense_probe_list<-function(probe_list) {
	probes=matrix(ncol=3,nrow=length(probe_list))
	n=0
	for (probe in probe_list) {
		n=n+1
		probes[n,]=split_probe_name(probe)
	}
	probeset_count<-table(probes[,2])
	final_probe<-vector()
	for (k in c(1:length(probe_list))) {
		if (probeset_count[probes[k,2]]>=2) {
			final_probe=c(final_probe,probe_list[k])
		}
	}
	return(final_probe,length(final_probe),sum(probeset_count>=2))		
}

#correlation cutoff
#initial cut 0.2 higher than the correlation cutoff

corr_cutoff=as.numeric(read.table("cutoff.txt"))
cut_height=1-corr_cutoff+0.2

files<-dir("out/",pattern=".txt")

for (n in 1:length(files)) {
	file=files[n]
	len_filename=nchar(file)
	file_prefix=substr(file,0,len_filename-4)

	probe_all<-read.table(paste("out/",file,sep=""),header=T,sep="\t",row.names=1)
	
	#avoid probes with all zero intensities
	probe_all<-probe_all+rnorm(dim(probe_all)[1]*dim(probe_all)[2])
	
	n_core_probe=0
	core_probe=vector()
	probe_name<-rownames(probe_all)
	array_names<-colnames(probe_all)
	
	n_array=length(array_names)
	if (length(probe_name)==0) {
		next
	}


	for (p in 1:length(probe_name)) {
		if (substring(probe_name[p],nchar(probe_name[p])-3,nchar(probe_name[p]))=='core') {
		n_core_probe=n_core_probe+1
		core_probe[n_core_probe]=probe_name[p]
	}
	}
	
	if (length(core_probe)<=4) {
		next
	}
	
	probe<-probe_all[core_probe,]


	probe_corr<-cor(t(probe))
	probe_corr_dist<-as.dist(1-probe_corr)
	probe_tree<-hclust(probe_corr_dist,method='average')

	
	
	cut_tree<-cutree(probe_tree,h=cut_height)
	major_cluster<-major_group(cut_tree)
	group_size=length(major_cluster)
	
	summary<-condense_probe_list(major_cluster)
	probe_number=summary[[2]]
	exon_number=summary[[3]]

	selected_probes=summary[[1]]
	final_probes=vector()


	theta1<-as.integer(fit.li.wong(t(probe))$theta+0.5)


	probe_theta_corr<-t(t(rep(-1,n_core_probe)))
	rownames(probe_theta_corr)<-core_probe

	iter_flag=1
	n_iter=0

	while (iter_flag==1 && n_iter<=50 && length(selected_probes)>=6) {
		theta2<-as.integer(fit.li.wong(t(probe[selected_probes,]))$theta+0.5)
		n_iter<-n_iter+1
		final_probes=vector()
		for (p in 1:length(core_probe)) {
			probe_intensity=as.vector(t(probe[core_probe[p],]))
			
			#if the probe intensity is the same in all samples (e.g. all 0)
			if (var(probe_intensity)<1e-12) {
				probe_theta_corr[p]=0
			}
			else {
				probe_theta_corr[p]=cor(probe_intensity,theta2,method='pearson')
			}
			if (probe_theta_corr[p]>=corr_cutoff) {
				final_probes=c(final_probes,core_probe[p])
			}
		}

		
		if (paste(final_probes,collapse="")==paste(selected_probes,collapse="")) {
			iter_flag=0
			}
		else {
			selected_probes=final_probes
			}
	}

	header=matrix(c("TranscriptCluster",array_names),ncol=n_array+1,nrow=1)

	if (n==1) {
		write.table(header,file='summary.core',quote=F,row.names=F,col.names=F,sep="\t",append=F)
		write.table(header,file='summary.selected',quote=F,row.names=F,col.names=F,sep="\t",append=F)
	}

	if (length(final_probes)<=6) {
		final_probes=vector()
	}
	write.table(matrix(c(file_prefix,t(final_probes)),nrow=1),file='Probes.Selected',quote=F,row.names=F,col.names=F,append=T)
	write.table(matrix(c(file_prefix,n_iter,length(final_probes),n_core_probe,length(final_probes)/n_core_probe),ncol=5,nrow=1),file='Probe_Selection_Summary',quote=F,row.names=F,col.names=F,append=T)
	
	


	
	if (length(final_probes)<=6) {
		write.table(matrix(c(file_prefix,theta1),ncol=n_array+1,nrow=1,byrow=T),file='summary.selected',quote=F,row.names=F,col.names=F,sep="\t",append=T)
		write.table(matrix(c(file_prefix,theta1),ncol=n_array+1,nrow=1,byrow=T),file='summary.core',quote=F,row.names=F,col.names=F,sep="\t",append=T)
		next
	}


	


	write.table(matrix(c(file_prefix,theta2),ncol=n_array+1,nrow=1,byrow=T),file='summary.selected',quote=F,row.names=F,col.names=F,sep="\t",append=T)
	write.table(matrix(c(file_prefix,theta1),ncol=n_array+1,nrow=1,byrow=T),file='summary.core',quote=F,row.names=F,col.names=F,sep="\t",append=T)

}
