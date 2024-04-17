library('leiden')
# library('biomaRt')
library('spatstat')
library('igraph')
library('viridis')
library('ComplexHeatmap')
library('seriation')
library('ggspavis')

celltype_colours = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Definitive endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     "Gut tube" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#D2EAF6",
                     "Splanchnic mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Lateral plate mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     "Mixed mesenchymal mesoderm" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "Erythroid" = "#f79083",
                     "Blood progenitors" = "#f9decf",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A",
                     
                     "Unknown" = "#FFFFFF",
                     "Low quality" = "#e6e6e6",
                     
                     # somitic and paraxial types
                     # colour from T chimera paper Guibentif et al Developmental Cell 2021
                     "Cranial mesoderm" = "#77441B",
                     "Anterior somitic tissues" = "#F90026",
                     "Sclerotome" = "#A10037",
                     "Dermomyotome" = "#DA5921",
                     "Posterior somitic tissues" = "#E1C239",
                     "Presomitic mesoderm" = "#9DD84A"
)

.palette1   = c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#DABEaf", 
    "#FF7F00", "#DA5900", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", 
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", "#f790af", 
    "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",  "#00bfdd",
    "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00", "#808044"
    )


.palette2   = c(
    "#F90026", "#FF7F00", "#FDB462", "#BEAED4", '#7F6874', "#77441B", "#9DD8dd",
    "#7BAFDE", "#EF5A9D", "#B2DF8A",  '#8EC792', "#C3C388", "#9DD84A", "#ffff33",
    "#B17BA6", "#E7298A", "#33A02C", '#354E23', "#f7f79e", "#f79083", "#f7909f",
    "#55A1B1", "#A6761D", "#FACB12", "#DA5921", "#f9decf", '#1e9099', "#aeaecc",
    "#8DD3C7", "#F397C0", "#532C8A",  "#8870ad", "#0F4A9C", "#CDE088", "#aa8299")

.palette3 = c(
    "#666666", "#999999", "#aa8282", "#d4b7b7", 
    "#8600bf", "#ba5ce3", "#808000", "#aeae5c", 
    "#1e90ff", "#00bfff", "#56ff0d", "#ffff00", "#DABE99")


.palette_all = c(.palette1, .palette2, .palette3) 
.palette_all_unique = .palette_all %>% unique

method_ord  = c('true', 'scPotter', 'tangram_cells', 'novosparc', 'scPotter_markers', 'exprs')
method_cols = c("#E1C239", '#33A02C', '#FF7F00', '#882E72', '#B2DF8A', '#666666')
names(method_cols) = method_ord

# mart = useMart('ensembl', dataset='hsapiens_gene_ensembl', host='www.ensembl.org')

# go_cellcycle = getBM(
#     attributes = c('ensembl_gene_id','external_gene_name'), 
#     filters    = 'go', 
#     values     = 'GO:0007049', 
#     mart       = mart)

# go_translation = getBM(
#     attributes = c('ensembl_gene_id','external_gene_name'), 
#     filters    = 'go', 
#     values     = 'GO:0006412', 
#     mart       = mart)

# go_ribosome1 = getBM(
#     attributes = c('ensembl_gene_id','external_gene_name'), 
#     filters    = 'go', 
#     values     = 'GO:0005840', 
#     mart       = mart)

# go_ribosome2 = getBM(
#     attributes = c('ensembl_gene_id','external_gene_name'), 
#     filters    = 'go', 
#     values     = 'GO:0042254', 
#     mart       = mart)

# ex_genes = unique(c(
#     go_cellcycle$external_gene_name, 
#     go_translation$external_gene_name, 
#     go_ribosome1$external_gene_name,
#     go_ribosome2$external_gene_name)) 


nngraph_comm <- function(nngraph, min_cells = 100, res_param=0.1){
    comm_vec = leiden(nngraph, resolution_parameter=res_param)
    comm_dt  = table(comm_vec)
    to.keep  = max(as.numeric(names(comm_dt)[comm_dt > min_cells]))
    comm_vec = ifelse(comm_vec <= to.keep, comm_vec, NA)
    names(comm_vec) = names(V(nngraph))
    comm_vec
}

get_network <- function(data, rho = .1, threshold = .1){
    data     = data.table(data)
    S        = stats::cov.wt(data, method='ML')$cov
    C        = stats::cov2cor(S)
    res      = glasso::glasso(C, rho=rho)
    AM       = abs(res$wi) > threshold
    diag(AM) = F
    rownames(AM) = colnames(AM) = colnames(data)
    g.lasso  = graph_from_adjacency_matrix(AM)
    # names(N(g.lasso)) = colnames(data)
    as(g.lasso, 'igraph')
}

analyze_network <- function(graph_obj, res = 1){
    adj_mat = igraph::as_adj(graph_obj)
    comm  = leiden(adj_mat, resolution_parameter=res)
    isolated_comm = min(which(table(comm) == 1))
    names(comm) = colnames(adj_mat)
    comm = ifelse(comm >= isolated_comm, isolated_comm, comm)

    comm_dt = data.table(ID=names(comm), community=comm) %>%
        setkey(community) %>% 
        .[, color := c(.palette1, .palette2)[community]] %>%
        setkey(ID) %>%
        .[names(V(graph_obj))]

    comm_dt
    # markers = unique(comm) %>% 
    #     map(~names(base::sort(colSums(adj_mat[names(comm)[comm==.x], names(comm)[comm==.x]]), decreasing=T))[1:2]) %>%
    #     unlist

}

graph_col_comm <- function(graph, lay, grp, sz, title=NULL, labels){
    igraph::V(graph)$color <- grp
    v <-  igraph::V(graph)
    # sprintf(comm_out, title) %>% pdf()
    p = plot.igraph(
        graph,
        vertex.size = 4,
        layout = lay,
        vertex.label = labels,
        vertex.frame.color = igraph::V(graph)$color,
        # vertex.label.family = 'Helvetica', 
        vertex.label.family = "sans",
        vertex.label.dist = 0,
        vertex.label.cex = 2,
        # vertex.label.font = 5,
        vertex.label.color = '#695c59',
        main=title,
        edge.arrow.mode='-',
        directed=FALSE,
        edge.arrow.size=0)
    # dev.off()
    p
}

plot_spatial <- function(dim_df, labels, label_cols=.palette1, title='', label_title='label', hide_legend=TRUE, sz=0.5){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2])) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        geom_point(size=sz) +
        theme_void() + 
        theme(axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=label_cols, na.value='gray' , drop = TRUE, limits = unique(labels)) +
        labs(title=title, x='', y='', color=label_title) 
    if(hide_legend)
        dim_plot = dim_plot + theme(legend.position='none')
    dim_plot
}

plot_2d <- function(dim_df, labels, label_cols=.palette1, title='', label_title='label', hide_legend=TRUE, sz=1, alpha=1){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2]), alpha=alpha) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        geom_point(size=sz, alpha=alpha) +
        theme_bw() + 
        theme(axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=label_cols, na.value='gray' , drop = TRUE, limits = unique(labels)) +
        labs(title=title, x='', y='', color=label_title)
    if(hide_legend)
        dim_plot = dim_plot + theme(legend.position='none')
    dim_plot
}


plot_2d_cont <- function(dim_df, labels, label_cols=nice_cols, title='', label_title='label', sz=3){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2])) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        # geom_hex(bins = 30) + 
        geom_point(size=sz) +
        # coord_fixed() +
        scale_color_viridis(na.value='gray') +
        theme_bw() + 
        theme(legend.position = "none",
            axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
        labs(title=title, x='', y='', color=label_title)
    dim_plot
}


hm_col_act <- function(
    mtx, 
    comm_dt,
    col_dt,
    cluster_rows=FALSE,
    cluster_columns=TRUE) {
    mtx = t(mtx)
    stopifnot(dim(mtx)[1] == dim(comm_dt)[1])
    stopifnot(dim(mtx)[2] == dim(col_dt)[1])
    col_dt = col_dt[colnames(mtx)]
    col_cols = col_dt$color
    names(col_cols) = col_dt$group
    # make column annotations
    col_annots  = HeatmapAnnotation(
        group                = col_dt$group, 
        col                  = list(group  = col_cols),
        annotation_name_side = 'left', 
        show_legend          = c(group=FALSE)
    )

    row_cols = comm_dt$color
    names(row_cols) = comm_dt$community
    row_annots  = rowAnnotation(
        gene_module          = comm_dt$community, 
        col                  = list(gene_module  = row_cols),
        show_legend          = c(gene_module=FALSE)
    )

    grp_vals  = seq(min(mtx), max(mtx), length.out=9)
    grp_cols  = circlize::colorRamp2(grp_vals, viridis::viridis(9))
    if(dim(as.matrix(mtx))[2] > 1){
        seriate_obj  = seriate(as.matrix(mtx - min(mtx)), method = "BEA_TSP")
        row_order    = get_order(seriate_obj, 1)
        column_order = get_order(seriate_obj, 2)
    }else{
        row_order    = 1:length(mtx)
        column_order = 1
    }
    
    # do heatmap for these genes
    hm_obj      = Heatmap(
        matrix=mtx, col=grp_cols,
        row_order=row_order, 
        cluster_rows=cluster_rows, cluster_columns=cluster_columns,
        row_names_gp=gpar(fontsize = 8), column_names_gp=gpar(fontsize = 10),
        name="Scaled\nGene Expression", 
        row_names_side="left",
        column_names_side="top",
        row_split=comm_dt$community,
        cluster_column_slices=FALSE,
        top_annotation=col_annots, left_annotation=row_annots
        )

    return(hm_obj)
}



hm_col_act <- function(
    mtx, 
    comm_dt,
    col_dt,
    cluster_rows=FALSE,
    cluster_columns=TRUE) {
    mtx = t(mtx)
    stopifnot(dim(mtx)[1] == dim(comm_dt)[1])
    stopifnot(dim(mtx)[2] == dim(col_dt)[1])
    col_dt = col_dt[colnames(mtx)]
    col_cols = col_dt$color
    names(col_cols) = col_dt$group
    # make column annotations
    col_annots  = HeatmapAnnotation(
        group                = col_dt$group, 
        col                  = list(group  = col_cols),
        annotation_name_side = 'left', 
        show_legend          = c(group=FALSE)
    )

    row_cols = comm_dt$color
    names(row_cols) = comm_dt$community
    row_annots  = rowAnnotation(
        gene_module          = comm_dt$community, 
        col                  = list(gene_module  = row_cols),
        show_legend          = c(gene_module=FALSE)
    )

    grp_vals  = seq(min(mtx), max(mtx), length.out=9)
    grp_cols  = circlize::colorRamp2(grp_vals, viridis::viridis(9))
    if(dim(as.matrix(mtx))[2] > 1){
        seriate_obj  = seriate(as.matrix(mtx - min(mtx)), method = "BEA_TSP")
        row_order    = get_order(seriate_obj, 1)
        column_order = get_order(seriate_obj, 2)
    }else{
        row_order    = 1:length(mtx)
        column_order = 1
    }
    
    # do heatmap for these genes
    hm_obj      = Heatmap(
        matrix=mtx, col=grp_cols,
        row_order=row_order, 
        cluster_rows=cluster_rows, cluster_columns=cluster_columns,
        row_names_gp=gpar(fontsize = 8), column_names_gp=gpar(fontsize = 10),
        name="Scaled\nGene Expression", 
        row_names_side="left",
        column_names_side="top",
        row_split=comm_dt$community,
        cluster_column_slices=FALSE,
        top_annotation=col_annots, left_annotation=row_annots
        )

    return(hm_obj)
}


act_colors <- function(data, var_vec){
    grp    = list()
    for(l in unique(var_vec)){
        idx = which(var_vec == l)
        dn  = data[idx, ]
        if(length(idx) == 1)
            grp[[l]] = dn
        else
            grp[[l]] = colMeans(dn)
    }
    if(length(grp) > 1){
        grp = grp %>% do.call(rbind, .)
        grp = apply(grp, 2, function(x) (x)/(max(x)))
    }else{
        grp = unlist(grp) %>% as.matrix %>% t
        rownames(grp) = unique(var_vec)
    }
    grp[is.infinite(grp)] = 0
    grp
}

graph_col_comm <- function(graph, lay, grp, sz=2, title=NULL, labels){
    igraph::V(graph)$color <- grp
    v <-  igraph::V(graph)
    # sprintf(comm_out, title) %>% pdf()
    p = plot.igraph(
        graph,
        # vertex.size = sz,
        layout = lay,
        vertex.label = labels,
        # vertex.frame.color = igraph::V(graph)$color,
        vertex.label.family = 'Helvetica',
        vertex.label.dist = 0,
        # vertex.label.cex = 1,
        # vertex.label.font = 5,
        vertex.label.color = '#585c59',
        main=title)
    # dev.off()
    p
}

graph_col_act <- function(graph, grp, lay, title){
    # sprintf(graph_out, paste0(title, '_train')) %>% pdf()
    grp_range = c(min(grp)^(1)/sum(grp^(1)),max(grp)^(1)/sum(grp^(1)))
    grp_vals  = seq(grp_range[1],grp_range[2],length.out=9)
    grp_cols  = circlize::colorRamp2(grp_vals, viridis::viridis(9))
    igraph::V(graph)$color = grp_cols(grp^(1)/sum(grp^(1)))
    p = plot.igraph(graph,
        vertex.size = 5,
        layout = lay,
        vertex.frame.color = igraph::V(graph)$color,
        vertex.label = "",
        main=title)
    p
}

plot_2d <- function(dim_df, labels, label_cols=nice_cols, title='', label_title='label', hide_legend=FALSE){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2])) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        geom_point(size=2) +
        theme_bw() + 
        theme(axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
        scale_color_manual(values=label_cols, na.value='gray') +
        labs(title=title, x='', y='', color=label_title)
    if(hide_legend)
        dim_plot = dim_plot + theme(legend.position='none')
    dim_plot
}

plot_2d_cont <- function(dim_df, labels, label_cols=nice_cols, title='', label_title='label'){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2])) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        # geom_hex(bins = 30) + 
        geom_point(size=2) +
        # coord_fixed() +
        scale_color_viridis() +
        theme_bw() + 
        theme(legend.position = "none",
            axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
        labs(title=title, x='', y='', color=label_title)
    dim_plot
}


remove_empty_spots <- function(spe){
    spe = spe[,colSums(counts(spe)) != 0]
    spe
}

remove_isolated_nodes <- function(ig){
    Isolated = which(igraph::degree(ig)==0)
    delete.vertices(ig, Isolated)
}