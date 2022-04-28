#######################################################################
## Title: Peak calling
## Author: Balaji Santhanam
####### if you find errors or have questions/comments, please email balaji.srinivasan.santhanam@gmail.com

######## R version 3.6.3 (2020-02-29) with required R libraries installed

######## to run this script, use the command 
##       Rscript --vanilla parsePeaks_annotated.R <path_to_ip.bam> <path_to_control.bam> <<path_to_gff_file>


######## to install required libraries
# NOTE: GenomicAlignments depends on Rsamtools
# if (!require("BiocManager", quietly = TRUE)) ; install.packages("BiocManager") ; BiocManager::install("Rsamtools")
# if (!require("BiocManager", quietly = TRUE)) ; install.packages("BiocManager") ; BiocManager::install("GenomicAlignments")

####
args = commandArgs(trailingOnly=TRUE)
if (length(args) <= 2) { stop("Please ensure paths to bam files for IP experiment, matched control as well as the appropriate GFF3 file have been specified.") }

options(echo=TRUE)

bamExpt = args[1] ; bamCtrl = args[2] ; gff_fi = args[3]

startTime = Sys.time()

getNullvec <- function(raw_expt_fragvec, sampling_number=100){
	n = sampling_number * length(raw_expt_fragvec)
	v = sample(raw_expt_fragvec,n, replace=T)
	vs = roll_sum(v,(2*halfwindow_size))
	vs = vs[seq(1,length(vs),(halfwindow_size))]
	idx = sample(length(vs) - (length(vs) %% 2))
	idx1 = idx[1:as.integer(0.5*length(idx))] ; idx2 = idx[(1+as.integer(0.5*length(idx))):length(idx)]
	val_ret = log((vs[idx1]/vs[idx2]),2)
}

tryCatch ( { library('GenomicAlignments') }, error = function(e) { install.packages("BiocManager") ; BiocManager::install("GenomicAlignments") })
tryCatch ( { library('RcppRoll') }, error = function(e) { install.packages('RcppRoll', dependencies=T) })


gff = read.csv(gff_fi, sep='\t', header=F, comment.char = '#')
gff$gene = unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(gff$V9), 'gene='), '[',2)),';'),'[',1))
gff_cds = gff[which(gff$V3 == 'CDS'),]
gff_cds.gr = GRanges(gff_cds$V1, IRanges(gff_cds$V4, gff_cds$V5))

##### index and read bam files and calculate scaling factors
indexBam(bamExpt) ; expt_summary = readGAlignments(bamExpt)
indexBam(bamCtrl) ; ctrl_summary = readGAlignments(bamCtrl)

e_name = gsub('.bam','',gsub('\\[mapped_reads_', '', gsub('\\].bam', '', gsub('_bam', '', bamExpt))))
c_name = gsub('.bam','',gsub('\\[mapped_reads_', '', gsub('\\].bam', '', gsub('_bam', '', bamCtrl))))
read_length = rlen = min(c(as.integer(names(which.max(table(granges(expt_summary)@ranges@width)))), 
							as.integer(names(which.max(table(granges(ctrl_summary)@ranges@width))))))
chr_lengths = expt_summary@seqinfo@seqlengths
names(chr_lengths) = expt_summary@seqinfo@seqnames
chr_list = names(chr_lengths)
expt_scale = 10^6/(length(expt_summary))
ctrl_scale = 10^6/(length(expt_summary))
rm(expt_summary) ; rm(ctrl_summary)


##################################### parameters for peak calling ####################################################
## fragment length
## size of individual fragments that the chromosome 
#  is broken into
## to account for local non-specific binding 
frag_length = 2000 
######################################################################################################################
## half window size
## size of windows within individual fragments to 
#  assess enrichment of reads in ip compared to control
halfwindow_size = 50
######################################################################################################################
## coverage threshold
## coverage indicates number of reads supporting a nucleotide
## goal is to minimize biases from experimental/protocol 
#  artefacts (eg. from PCR)
## value should strike a balance the requirement to minimize protocol-related 
#  bias while making sure enough dynamic range exists for peak detection
cov_thresh = 5
######################################################################################################################


out_filename = paste('ChIPseq_compare_', e_name, '__vs__', c_name, '__', 'fragment_length-', frag_length, '__window-',halfwindow_size,'__coveragethresh-', cov_thresh,'__',Sys.Date(), sep='')

output_matrix = c()
for(ch in chr_list) {
	ch_time1 = Sys.time()
	chr_index = frag_id = start_id = stop_id = enrich_for = enrich_rev = fdr_for = fdr_rev = fdr_for_both = fdr_rev_both = c()
	support_for = support_rev = max_cov_for = max_cov_rev = c()
	ch_l = as.numeric(chr_lengths[ch])
	ef = er = cf = cr = rep(0, ch_l)
	names(ef) = names(er) = names(cf) = names(cr) = seq(1, ch_l)
	param = ScanBamParam(what = c('qname','rname','strand','pos'), which=GRanges(ch,IRanges(c(1),c(ch_l))))
	expt_chr_cov = readGAlignments(bamExpt,use.names=TRUE,param=param)
	ctrl_chr_cov = readGAlignments(bamCtrl,use.names=TRUE,param=param)
	##### '+' strand
	expt_chr_reads = start(ranges(expt_chr_cov[expt_chr_cov@strand == '+']))
	ctrl_chr_reads = start(ranges(ctrl_chr_cov[ctrl_chr_cov@strand == '+']))
	ef[names(table(expt_chr_reads))] = as.vector(table(expt_chr_reads))
	cf[names(table(ctrl_chr_reads))] = as.vector(table(ctrl_chr_reads))
	##### '-' strand
	expt_chr_reads = end(ranges(expt_chr_cov[expt_chr_cov@strand == '-']))
	ctrl_chr_reads = end(ranges(ctrl_chr_cov[ctrl_chr_cov@strand == '-']))
	er[names(table(expt_chr_reads))] = as.vector(table(expt_chr_reads))
	cr[names(table(ctrl_chr_reads))] = as.vector(table(ctrl_chr_reads))
	for(id in seq(1, (as.integer(chr_lengths[ch])), frag_length)){
		frag_start = id ; frag_end = min(c(as.numeric(chr_lengths[ch]), (frag_start + frag_length)))
		frag_name = paste(frag_start, frag_end,sep=':')
		frag_ef = ef[frag_start:frag_end] ; frag_er = er[frag_start:frag_end]
 	   	frag_cf = cf[frag_start:frag_end] ; frag_cr = cr[frag_start:frag_end]
		frag_length_act = length(frag_ef)
	    tmp_v = frag_ef
    	support_f = roll_sum(tmp_v, (2*halfwindow_size))
    	max_cov_f = c()
    	for(ix_cov in seq(1, (frag_length_act - 2*halfwindow_size - 1))) { end_coord = min(c(frag_length_act, (ix_cov+2*halfwindow_size))) ; max_cov_f = c(max_cov_f, id -1 + which.max(tmp_v[ix_cov:end_coord])) }
    	tmp_v = frag_er
    	support_r = roll_sum(tmp_v, (2*halfwindow_size))
    	max_cov_r = c()
    	for(ix_cov in seq(1, (frag_length - 2*halfwindow_size - 1))) { end_coord = min(c(frag_length_act, (ix_cov+2*halfwindow_size))) ; max_cov_r = c(max_cov_r, id -1 + which.max(tmp_v[ix_cov:end_coord])) }
    	frag_ef[frag_ef > cov_thresh] = cov_thresh ; frag_er[frag_er > cov_thresh] = cov_thresh
    	frag_cf[frag_cf > cov_thresh] = cov_thresh ; frag_cr[frag_cr > cov_thresh] = cov_thresh
    	
    	print(paste(ch, frag_name, sep='    '))
    	vec_ef = roll_sum(frag_ef, 2*halfwindow_size) ; vec_er = roll_sum(frag_er, 2*halfwindow_size)
    	vec_cf = roll_sum(frag_cf, 2*halfwindow_size) ; vec_cr = roll_sum(frag_cr, 2*halfwindow_size)
    	ind_ = seq(1,length(vec_ef),halfwindow_size)
    	vec_ef = vec_ef[ind_] ; vec_er = vec_er[ind_] ; vec_cf = vec_cf[ind_] ; vec_cr = vec_cr[ind_]
    	support_f = support_f[ind_] ; support_r = support_r[ind_] ; max_cov_f = max_cov_f[ind_] ; max_cov_r = max_cov_r[ind_]
	    dist_f = log(((expt_scale*vec_ef)/(ctrl_scale*vec_cf)),2) ; dist_r = log(((expt_scale*vec_er)/(ctrl_scale*vec_cr)),2)
    	dist_f[is.na(dist_f)] = dist_r[is.na(dist_r)] = 0

######################################################################################################################
################################################# NULL distributions #################################################
    	null_vecf = getNullvec((expt_scale*(frag_ef))) ; null_vecr = getNullvec((expt_scale*(frag_er)))
    	null_vecf_both = getNullvec(c((expt_scale*(frag_ef)), (ctrl_scale*(frag_cf))))
    	null_vecr_both = getNullvec(c((expt_scale*(frag_er)), (ctrl_scale*(frag_cr))))
    	null_distf = null_vecf ; null_distr = null_vecr
    	null_distf[is.na(null_distf)] = null_distr[is.na(null_distr)] = 0
    	null_distf_both = null_vecf_both ; null_distr_both = null_vecr_both
    	null_distf_both[is.na(null_distf_both)] = null_distr_both[is.na(null_distr_both)] = 0
######################################################################################################################
################################################# FDR calculations ###################################################
	    fdr_f = p.adjust(apply(data.frame(dist_f), 1, function(val,null_vec){sum(null_distf >= val)/length(null_distf)}),'fdr')
    	fdr_r = p.adjust(apply(data.frame(dist_r), 1, function(val,null_vec){sum(null_distr >= val)/length(null_distr)}),'fdr')
	    fdr_fb = p.adjust(apply(data.frame(dist_f), 1, function(val,null_vec){sum(null_distf_both >= val)/length(null_distf_both)}),'fdr')
    	fdr_rb = p.adjust(apply(data.frame(dist_r), 1, function(val,null_vec){sum(null_distr_both >= val)/length(null_distr_both)}),'fdr')

    	dist_f[dist_f == Inf] = median(dist_f[is.finite(dist_f)]) ; dist_f[dist_f == -Inf] = min(dist_f[is.finite(dist_f)]) 
    	dist_r[dist_r == Inf] = median(dist_r[is.finite(dist_r)]) ; dist_r[dist_r == -Inf] = min(dist_r[is.finite(dist_r)])
    	start_i = ind_
    	start_i = start_i + frag_start - 1
    	stop_i = start_i + 2*halfwindow_size
    	max_cov_f = max_cov_f + start_i
    	max_cov_r = max_cov_r + start_i
######################################################################################################################
##################### internal checks for debugging (can be removed)  ################################################
		stopifnot(length(vec_ef) == length(vec_cr))
		stopifnot(length(start_i) == length(stop_i))
		stopifnot(length(start_i) == length(fdr_f))
		stopifnot(length(start_i) == length(fdr_fb))
		stopifnot(length(stop_i) == length(vec_er))
		stopifnot(length(stop_i) == length(fdr_r))
		stopifnot(length(stop_i) == length(fdr_rb))
		stopifnot(length(support_f) == length(dist_r))
		stopifnot(length(support_r) == length(dist_f))
		stopifnot(length(support_f) == length(max_cov_r))
		stopifnot(length(support_r) == length(max_cov_f))
######################################################################################################################
		ll = length(fdr_f)
		chr_index = c(chr_index, rep(ch,ll))
		frag_id = c(frag_id, rep(frag_name, ll))
		start_id = c(start_id,start_i)
		stop_id =c(stop_id, stop_i)
		enrich_for = c(enrich_for, dist_f)
		enrich_rev = c(enrich_rev, dist_r)
		fdr_for = c(fdr_for, fdr_f)
		fdr_rev = c(fdr_rev, fdr_r)
		fdr_for_both = c(fdr_for_both, fdr_fb)
		fdr_rev_both = c(fdr_rev_both, fdr_rb)
		support_for = c(support_for, support_f)
		support_rev = c(support_rev, support_r)
		max_cov_for = c(max_cov_for, max_cov_f)
		max_cov_rev = c(max_cov_rev, max_cov_r)
	}
	output_matrix = as.data.frame(rbind(output_matrix, data.frame(chr = chr_index, gen_frag = frag_id, start = start_id, stop = stop_id, 
	enrich_for = enrich_for, enrich_rev = enrich_rev, fdr_for = fdr_for, fdr_rev = fdr_rev,
	read_support_for = support_for, read_support_rev = support_rev, max_coverage_for = max_cov_for,
	max_coverage_rev = max_cov_rev, fdr_for_both = fdr_for_both, fdr_rev_both = fdr_rev_both)))
	ch_time2 = Sys.time()
	ch_diff = difftime(ch_time2, ch_time1, units='mins')
	print(paste('completed',ch,': ',ch_diff,' mins',sep=''))
}
write.table(output_matrix, out_filename, sep='\t', row.names=F, quote=F)
endTime = Sys.time()
print('Output files written to:')
print(out_filename)
tdiff = difftime(endTime, startTime, units='hours')
print(paste('Analysis completed in', tdiff,' hours',sep=' '))


