library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)




dnaseCompartments <- function(bw.files, chr="chr22", resolution = 100*1000,
 method = c("pearson", "spearman"), normalize = TRUE, gcCorrection = TRUE,
  keep=TRUE, verbose = TRUE){
	method <- match.arg(method)

	# Extraction of the score matrix:
	if (verbose){
		cat("[dnaseCompartments] Extraction of the score matrix \n")
	}
	scoreMatrix <- .getScoreMatrix(bw.files = bw.files, chr = chr,
	 resolution = resolution, normalize = normalize)

	# GC correction:
	if (gcCorrection){
		if (verbose){
			cat("[dnaseCompartments] Correction for GC content \n")
		}
		scoreMatrix <- gcCorrect(scoreMatrix)
	}

	# Correlation matrix:
	if (verbose){
		cat("[dnaseCompartments] Creation of the correlation matrix \n")
	}
	gr <- createCorMatrix(scoreMatrix, method = method)
	rm(scoreMatrix)

	# Extraction of the compartments:
	if (verbose){
		cat("[dnaseCompartments] Extraction of the AB compartments \n")
	}
	gr <- extractAB.dnase(gr = gr, keep = kepp)

	# Centering:
	gr$pc <- gr$pc - median(gr$pc)
	gr
}



getScoreMatrix <- function(bw.files, chr="chr22", 
	resolution=100*1000, version="hg19", normalize=TRUE){
	require(rtracklayer)
	require(BSgenome.Hsapiens.UCSC.hg19)

	chr.lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:22]
	chr.length <- chr.lengths[chr]

		# Normalization for library size:
	.dnaseNormalize <- function(scoreMatrix){
			covs <- colSums(scoreMatrix)/1000000
			t(t(scoreMatrix/covs))
	}
	gr.binned <- tileGenome(seqlengths = chr.length,
                              tilewidth = resolution, cut.last.tile.in.chrom = TRUE)
	gr.selection <- GRanges(seqnames=chr, IRanges(start=1, end=chr.length))

	M <- length(bw.files)
	N <- length(gr.binned)

	scoreMatrix <- do.call(cbind, lapply(1:M, function(m){
		temp <- import(bw.files[m], selection = BigWigSelection(gr.selection))
		ids <- findOverlaps(temp, gr.binned , select="first")
		score <- tapply(temp$score, ids, sum, na.rm=TRUE)
		temp <- rep(NA,N)
		temp[as.numeric(names(score))] <- score
		temp
	}))

	# Some clearning:
	good <- rowSums(is.na(scoreMatrix)) != ncol(scoreMatrix)
	scoreMatrix <- scoreMatrix[good,]
	gr.binned <- gr.binned[good]
	good <- complete.cases(scoreMatrix)
	scoreMatrix <- scoreMatrix[good,]
	gr.binned <- gr.binned[good]

	if (normalize){scoreMatrix <- .dnaseNormalize(scoreMatrix)}
	gr.binned$matrix <- scoreMatrix
	scoreMatrix <- gr.binned
	genome(scoreMatrix) <- "hg19"
	return(scoreMatrix)
}




# Correction for GC content:
gcCorrect <- function(scoreMatrix){
	require(Biostrings)
	require(BSgenome.Hsapiens.UCSC.hg19)

	# Make sure GRanges is a valid one:
	.cleanGr <- function(gr){
		lens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:22]
		chrs <- as.character(seqnames(gr))
		for (i in 1:22){
			indices <- which(chrs==paste0("chr", i))
			if (length(indices)!=0){
				ends <- end(gr[indices])
				aa <- which(ends>=lens[i])
				if (length(aa)!=0){
					end(gr[indices][aa]) <- lens[i]
				}
			}
		}
		start(gr)[start(gr)==0] <- 1
		gr
	}

	.gcCorrect <- function(matrix, gc){
		logMatrix <- log(matrix+1)
		logMatrix <- apply(logMatrix, 2, function(x){
			loess(x~gc)$res
		})
		logMatrix <- exp(logMatrix)-1
		logMatrix
	}


	# Calculating GC content for each bin:
	scoreMatrix <- .cleanGr(scoreMatrix)
	vi <- Views(Hsapiens, scoreMatrix)
	f  <- oligonucleotideFrequency(vi, width=1)
	gc <- rowSums(f[,c("C", "G")])/width(scoreMatrix)

	# Takes time..
	scoreMatrix$matrix <- .gcCorrect(scoreMatrix$matrix, gc)
	scoreMatrix
}

# Creation of the correlation matrix:
createCorMatrix <- function(scoreMatrix, method = c("pearson", "spearman")){
	scoreMatrix$cor.matrix <- cor(t(scoreMatrix$matrix), method = method)
	scoreMatrix$matrix <- NULL
	scoreMatrix
}







extractAB.dnase <- function(gr, keep = TRUE){
    if (! (is(gr, "GRanges") && "cor.matrix" %in% names(mcols(gr)))) {
        stop("'gr' must be an object created by createCorMatrix")
    }

    pc <- .getFirstPC(gr$cor.matrix)
    pc <- .meanSmoother(pc, iter=2)
    pc <- .unitarize(pc)
    ## Fixing sign of eigenvector
    gr$pc <- pc
    if (!keep) {
        gr$cor.matrix <- NULL
    }
    return(gr)
}



.getFirstPC <- function(matrix, ncomp = 1){
    ## Centering the matrix
    center <- rowMeans(matrix, na.rm = TRUE)
    matrix <- sweep(matrix, 1L, center, check.margin = FALSE)
    return(mixOmics::nipals(matrix, ncomp = ncomp)$p[,1])
}


.meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
    meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
        n <- length(x)
        y <- rep(NA,n)
        
        window.mean <- function(x, j, k, na.rm=na.rm){
            if (k>=1){
                return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
            } else {
                return(x[j])
            }    
        }
        
        for (i in (k+1):(n-k)){
            y[i] <- window.mean(x,i,k, na.rm)
        }
        for (i in 1:k){
            y[i] <- window.mean(x,i,i-1, na.rm)
        }
        for (i in (n-k+1):n){
            y[i] <- window.mean(x,i,n-i,na.rm)
        }
        y
    }
    
    for (i in 1:iter){
        x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
    }
    x
}

.unitarize <- function (x, medianCenter = TRUE) {
    if(medianCenter) {
        x <- x - median(x, na.rm = TRUE)
    }
    bad <- is.na(x)
    x[!bad] <- x[!bad] / sqrt(sum(x[!bad]^2))
    n.bad <- sum(bad)
    if (n.bad > 0){
        cat(sprintf("[.unitarize] %i missing values were ignored.\n", n.bad))
    }
    return(x)
}

























