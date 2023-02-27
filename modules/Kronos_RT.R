#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE,
        dplyr.summarise.inform=FALSE,
        warn = 1,
        scipen = 999)

option_list = list(
    make_option(
        c("-s", "--input_s"),
        type = "character",
        default = NULL,
        help = "Long form csv file of cell_id, chr, start, end, copy, state, reads for S-phase cells",
        metavar = "character"
    ),
    make_option(
        c("-g", "--input_g"),
        type = "character",
        default = NULL,
        help = "Long form csv file of cell_id, chr, start, end, copy, state, reads for G1/2-phase cells",
        metavar = "character"
    ),
    make_option(
        c("-m", "--metrics_s"),
        type = "character",
        default = NULL,
        help = "Metrics csv file for S-phase cells",
        metavar = "character"
    ),
    make_option(
        c("-M", "--metrics_g"),
        type = "character",
        default = NULL,
        help = "Metrics csv file for G1/2-phase cells",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = NULL,
        help = "Output file containing replication states for S-phase cells",
        metavar = "character"
    ),
    make_option(
        c("-f", "--output_file_base_name"),
        type = "character",
        default = "out",
        help = "Base name for the output file [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-c", "--cores"),
        type = "integer",
        default = 3,
        help = "Number of parallel jobs to run [default= %default]",
        metavar = "integer"
    ),
    make_option(
        c("--filter_cells"),
        type = "logical",
        default = F,
        action = "store_true",
        help = "Discard S-phase cells with low correlations in the final output [default= %default]",
        metavar = "logical"
    ),
    make_option(
        c("--min_correlation"),
        type = "double",
        default = 0.25,
        help = "Minimum correlation value between one cell and its best correlating cell for this cell to not be discarded [default= %default]",
        metavar = "double"
    )
)

print('parsing arguments')

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

#load libraries
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(doSNOW, quietly = TRUE))


#check inputs
if (!'input_s' %in% names(opt)) {
    stop("S-phase cells not provided in input. See script usage (--help)")
}

if (!'input_g' %in% names(opt)) {
    stop("G1/2-phase cells not provided in input. See script usage (--help)")
}

if (!'out' %in% names(opt)) {
    stop("Output file not provided. See script usage (--help)")
}


# load a tsv file for the S-phase cells
# columns should be Cell, chr, start, end, CN, reads
signal_smoothed <- read_csv(file.path(opt$input_s), col_types = cols(
    cell_id = col_character(),
    chr = col_character(),
    start = col_integer(),
    end = col_integer(),
    copy = col_double(),
    state = col_integer(),
    reads = col_double()
))

# load a tsv file for the G1/G2-phase cells
# columns should be Cell, chr, start, end, CN, reads
G1G2_smoothed <- read_csv(file.path(opt$input_g), col_types = cols(
    cell_id = col_character(),
    chr = col_character(),
    start = col_integer(),
    end = col_integer(),
    copy = col_double(),
    state = col_integer(),
    reads = col_double()
))

# convert column names to be compatible with Kronos
signal_smoothed$Cell = signal_smoothed$cell_id

G1G2_smoothed$Cell = G1G2_smoothed$cell_id
G1G2_smoothed$CN = G1G2_smoothed$copy

# save a copy of the original data
signal_smoothed_original <- signal_smoothed

# load the metrics file for the S-phase cells
# columns should be Cell, mean_ploidy, median_ploidy, Sphase_first_part, Sphase_second_part
metrics_s <- read_csv(file.path(opt$metrics_s), col_types = cols(
    cell_id = col_character(),
    mean_copy = col_double(),
    breakpoints = col_integer(),
    total_mapped_reads_hmmcopy = col_integer()
))

# load the metrics file for the G1/G2-phase cells
# columns should be Cell, mean_ploidy, median_ploidy, Sphase_first_part, Sphase_second_part
metrics_g <- read_csv(file.path(opt$metrics_g), col_types = cols(
    cell_id = col_character(),
    mean_copy = col_double(),
    breakpoints = col_integer(),
    total_mapped_reads_hmmcopy = col_integer()
))


# compute the median of all mean_copy values for G1/G2-phase cells
median_ploidy_G1_G2_cells = median(metrics_g$mean_copy)

# hard-code the default settings for S-phase first and second parts
Sphase_first_part = 0.95
Sphase_second_part = 0.55

# rename mean_copy to mean_ploidy
metrics_s = metrics_s %>%
    rename(mean_ploidy = mean_copy)
metrics_g = metrics_g %>%
    rename(mean_ploidy = mean_copy)

# rename cell_id to Cell
metrics_s = metrics_s %>%
    rename(Cell = cell_id)
metrics_g = metrics_g %>%  
    rename(Cell = cell_id)

# subset both metrics columns to just cell_id, mean_copy, breakpoints, total_mapped_reads_hmmcopy
metrics_s = metrics_s[,c('Cell', 'mean_ploidy', 'breakpoints', 'total_mapped_reads_hmmcopy')]
metrics_g = metrics_g[,c('Cell', 'mean_ploidy', 'breakpoints', 'total_mapped_reads_hmmcopy')]

# merge the metrics_s file with the signal_smoothed file
signal_smoothed = signal_smoothed %>%
    dplyr::inner_join(metrics_s, by = c('Cell'))


# pass in dummy column saying that all s-phase cells are noisy
signal_smoothed$is_noisy = TRUE
G1G2_smoothed$is_noisy = FALSE

print('merged metrics with signal smoothed')
print(head(signal_smoothed))


print('correcting mean ploidy for S-phase cells')
# correct mean ploidy late S phase
# find which cells are early or late S-phase which need correcting
signal_smoothed = signal_smoothed %>%
    mutate(
        mean_ploidy_corrected = ifelse(
            as.logical(is_noisy) == T &
                mean_ploidy < median_ploidy_G1_G2_cells,
            mean_ploidy / Sphase_second_part,
            ifelse(
                as.logical(is_noisy) == T &
                    mean_ploidy > median_ploidy_G1_G2_cells,
                mean_ploidy / Sphase_first_part,
                mean_ploidy
            )
        )
    )

# if the mean_ploidy has been corrected, then correct the copy number
signal_smoothed = signal_smoothed%>%
    dplyr::mutate(
        CN = ifelse(
            mean_ploidy == mean_ploidy_corrected,
            copy,
            copy * mean_ploidy_corrected / mean_ploidy
        )
    )

print(head(signal_smoothed))

# subset to only include columns Cell, chr, start, end, CN, reads
signal_smoothed = signal_smoothed[,c('Cell', 'chr', 'start', 'end', 'CN', 'reads')]
G1G2_smoothed = G1G2_smoothed[,c('Cell', 'chr', 'start', 'end', 'CN', 'reads')]

# create a list of unique values in the chr column
chr_list = unique(signal_smoothed$chr)

# add an index column that represents a 1:N ordering of each unique cell
signal_smoothed = signal_smoothed %>%
    mutate(index = as.numeric(factor(Cell)))

G1G2_smoothed = G1G2_smoothed %>%
    mutate(index = as.numeric(factor(Cell)))


# add dummny columns for basename and group
signal_smoothed$basename = 'A'
signal_smoothed$group = 1
G1G2_smoothed$basename = 'A'
G1G2_smoothed$group = 1


print('computing background G1/2 median CN profile')

#calculate background
# this is the median CN profile across all G1/2 cells
backgroud_smoothed=G1G2_smoothed%>%
    group_by(basename, group, chr, start, end)%>%
    summarise(background=median(CN))

# binarize data into Replicated or not replicated bins
Replication_state = function(Samples, background, chr_list,cores=3){

    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    #merge signal and bg and calculate their ratio
    Samples = Samples %>%
        ungroup() %>%
        mutate(chr = factor(x =  chr, levels = chr_list)) %>%
        inner_join(background,
                   by = c("chr", "start", "end", "basename", 'group')) %>%
        mutate(CN_bg = log2(CN / background)) %>%
        drop_na() %>%
        filter(is.finite(CN_bg))

    # identify threshold that minimizes the difference of the real data with a binary state (1 or 2)
    selecte_th = foreach(
        line = unique(Samples$basename),
        .combine = 'rbind',
        .packages = c('foreach', 'tidyverse')
    ) %dopar% {
        sub_sig = Samples %>%
            filter(basename == line)

        #identify range within looking for a CNV threshold to define replicated and not replicated values
        range = seq(0, 1, 0.01)

        th_temp = foreach(i = range,
                          .combine = 'rbind',
                          .packages = 'tidyverse') %do% {
                              summary = sub_sig %>%
                                  mutate(Rep = ifelse(CN_bg >= i, T, F),
                                         Error = (Rep - CN_bg) ^ 2) %>%
                                  group_by(Cell, index, basename, group)  %>%
                                  summarise(summary = sum(Error))

                              data.frame(
                                  th = i,
                                  basename = summary$basename,
                                  index = summary$index,
                                  sum_error = summary$summary
                              )
                          }
        th_temp
    }

    selecte_th = selecte_th %>%
        group_by(index, basename) %>%
        filter(sum_error == min(sum_error)) %>%
        summarise(th = min(th)) %>%
        ungroup()

    # mark replicated bins
    Samples = Samples %>%
        inner_join(selecte_th, by = c("index", "basename")) %>%
        mutate(Rep = ifelse(CN_bg >= th, T, F))

    #identify new distribution in the S phase based on the amount of replicated bins
    new_index_list = Samples %>%
        group_by(Cell, index, basename, group) %>%
        summarise(PercentageReplication = mean(Rep)) %>%
        ungroup() %>%
        arrange(PercentageReplication) %>%
        group_by(group) %>%
        mutate(newIndex = 1:n()) %>%
        dplyr::select(oldIndex = index, newIndex, Cell, basename, group,PercentageReplication)

    Samples = Samples %>%
        inner_join(new_index_list,
                   by = c('Cell', 'index' = 'oldIndex', 'basename', 'group'))

    stopCluster(cl)


    return(Samples)
}

print('computing replication states')
#merge signal and bg and calculate their ratio
signal_smoothed = Replication_state(Samples =signal_smoothed ,
                                    background = backgroud_smoothed,
                                    chr_list = chr_list,cores = opt$cores)

# remove control track
rm('backgroud_smoothed')

# write results here if we are not filtering cells
if (!opt$filter_cells){

    #write results
    signal_smoothed%>%
        dplyr::select(chr,
                    start,
                    end,
                    CN,
                    background,
                    CN_bg,
                    th,
                    Rep,
                    PercentageReplication,
                    Cell,
                    basename,
                    group,
                    newIndex)%>%
        write_delim(
            file = opt$out,
            delim = '\t',
            col_names = T
        )
} else {
    print('filtering cells')
    #matrix for correlation
    signal_smoothed = signal_smoothed %>%
        ungroup() %>%
        arrange(group, newIndex) %>%
        unite(index, c(group, newIndex), sep = ' _ ') %>%
        mutate(index = factor(index, levels = unique(index)))
    mat = signal_smoothed %>%
        unite(pos, c(chr, start), sep = ':') %>%
        mutate(Rep = as.numeric(Rep)) %>%
        dplyr::select(pos, index, Rep) %>%
        spread(key = index, value = Rep) %>%
        column_to_rownames('pos') %>%
        filter(complete.cases(.)) %>%
        as.matrix()

    #correlation similarity distances
    results = 1-as.matrix(dist.binary(t(mat),method = 2,diag = T,upper = T))
    basenames = str_remove(colnames(mat), ' _ [0-9]{1,10}$')
    Index = colnames(mat)
    basename_n = basenames

    for (i in 1:length(unique(basename_n))) {
        basename_n[basename_n == unique(basename_n)[i]] = i
    }

    #filter out cells that do not correlate
    to_keep = foreach(i = 1:length(unique(basename_n))) %do% {
        sub_mat = results[basename_n == i, basename_n == i]
        diag(sub_mat) = 0
        ! rowQuantiles(x = sub_mat, probs = 0.60,na.rm = T) <= opt$min_correlation
    }

    to_keep = unlist(to_keep)
    results = results[to_keep, to_keep]
    basename_n = basename_n[to_keep]
    Index = Index[to_keep]


    #filter out samples that don't correlate and save
    signal_smoothed = signal_smoothed %>%
        filter(index %in% Index) %>%
        separate(index, c('group', 'index'), sep = ' _ ') %>%
        mutate(group = factor(group, level = unique(opt$groups)))

    rep_percentage = signal_smoothed %>%
        group_by(Cell,basename,group, index) %>%
        summarise(Rep_percentage = mean(Rep))

    #new index
    new_index_list = rep_percentage %>%
        ungroup() %>%
        arrange(Rep_percentage) %>%
        group_by(group) %>%
        mutate(newIndex = 1:n()) %>%
        arrange(group,newIndex) %>%
        dplyr::select(oldIndex=index, newIndex,Cell,basename,group)
    

    signal_smoothed = signal_smoothed %>%
        ungroup() %>%
        inner_join(new_index_list, by = c('Cell','index'='oldIndex', 'basename','group')) %>%
        dplyr::select(chr,
                    start,
                    end,
                    CN,
                    background,
                    CN_bg,
                    th,
                    Rep,
                    PercentageReplication,
                    Cell,
                    basename,
                    group,
                    newIndex)

    write_delim(
        x = signal_smoothed,
        file = opt$out,
        delim = '\t',
        col_names = T
    )

    rm('new_index_list')
}

print('done')