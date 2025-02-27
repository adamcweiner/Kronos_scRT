#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE,
        dplyr.summarise.inform=FALSE,
        warn = 1,
        scipen = 999)

option_list = list(
    make_option(
        c("-F", "--File"),
        type = "character",
        default = NULL,
        help = "Replication timing files separated by a comma. Format: chr <TAB> start <TAB> end <TAB> group",
        metavar = "character"
    ),
    make_option(
        c("-s", "--sort"),
        type = "character",
        default = NULL,
        help = "Group names orders",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "output",
        help = "Output directory [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-f", "--output_file_base_name"),
        type = "character",
        default = "out",
        help = "Base name for the output file [default= %default]",
        metavar = "character"
    )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

#load module
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(GGally, quietly = TRUE))
suppressPackageStartupMessages(library(ggcorrplot, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))

theme_set(theme_bw())


if (!'File' %in% names(opt)) {
    stop("RT files were not provided. See script usage (--help)")
    
} else{
    opt$File = str_split(opt$File, ',')[[1]]
}
if ('sort' %in% names(opt)) {
    opt$sort = str_split(opt$sort, ',')[[1]]
    
}

#create directory
if (!dir.exists(opt$out)) {
    dir.create(opt$out,recursive = T)
}

scRT = foreach(
    i = 1:length(opt$File),
    .packages = 'tidyverse',
    .combine = 'rbind'
) %do% {
    tmp = read_tsv(opt$File[i], col_types = cols(chr = 'c'))
    if ('sort' %in% names(opt)) {
        tmp %>%
            mutate(group = factor(group, levels = opt$sort))
        
    } else{
        tmp
    }
}


scRT = scRT %>% spread(group, RT) %>%
    drop_na() %>%
    dplyr::select(-chr, -start, -end)

plot = ggcorrplot(
    scRT %>%
        cor(method = 'spearman'),
    lab = T,
    lab_col = 'red',
    legend.title = 'Spearman\ncorrelation',
    colors =  c('#BCAF6FFF', '#7C7B78FF', '#00204DFF'),
    ggtheme = ggplot2::theme(aspect.ratio = 1)
)

suppressMessages(ggsave(plot = plot,
                        filename = file.path(
                            opt$out,
                            paste0(opt$output_file_base_name,
                                   '_spearman_correlation.pdf')
                        )))

suppressMessages(ggsave(
    plot = ggpairs(
        scRT,
        diag = list(
            continuous = function(data, mapping, ...) {
                p <- ggplot(data, mapping) +
                    geom_density(aes(y = ..density.. / max(..density..)),
                                 fill = 'grey') +
                    scale_x_continuous(breaks = c(0, 0.5, 1)) +
                    scale_y_continuous(breaks = c(0, 0.5, 1))
                return(p)
            }
        ),
        upper = list(
            continuous = function(data, mapping, ...) {
                Spearman = cor(data[as_label(mapping$x)], data[as_label(mapping$y)], method = 'spearman')
                data = tibble(x = seq(0, 2 * pi, length.out = 200),
                              Corr = Spearman) %>%
                    mutate(y = 0.5 + Corr / 2 * sin(x),
                           x = 0.5 + Corr / 2 * cos(x))
                
                p <- ggplot(data, aes(x = x, y = y, fill = Corr)) +
                    geom_polygon() + theme(axis.title = element_blank(),
                                           panel.grid = element_blank()) +
                    annotate(
                        'text',
                        x = 0.5,
                        y = 0.5,
                        label = paste("Corr:", round(unique(Spearman), 3), sep = '\n'),
                        color = 'red'
                    ) +
                    scale_fill_gradient2(
                        low = '#BCAF6FFF',
                        high = '#00204DFF',
                        mid = '#7C7B78FF',
                        midpoint = 0,
                        limits = c(-1, 1)
                    ) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
                
                return(p)
            }
        ),
        lower = list(
            continuous = function(data, mapping, ...) {
                p <- ggplot(data = data, mapping = mapping) +
                    geom_hex(bins = 50, aes(fill = ..ndensity..)) +
                    scale_fill_gradientn(
                        'Density',
                        colours = c(
                            "#FFEA46FF",
                            "#D3C164FF",
                            "#A69D75FF",
                            "#7C7B78FF",
                            "#575C6DFF",
                            "#233E6CFF",
                            "#00204DFF"
                        )
                    ) +
                    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
                    scale_x_continuous(breaks = c(0, 0.5, 1)) +
                    scale_y_continuous(breaks = c(0, 0.5, 1)) +
                    geom_abline(slope = 1,
                                color = 'black',
                                alpha = 0.5)
                
                return(p)
            }
        ),
        legend = c(2, 1)
    ) + theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)
    ),
    filename = file.path(
        opt$out,
        paste0(opt$output_file_base_name,
               '_pair_scatter_plot_RTs.pdf')
    )
))

print('done')
