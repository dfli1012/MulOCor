#' @title correlation plot
#' @description tow omics correlation
#' @details Input two omics data matrix
#' @param data1 one omics matrix
#' @param data2 another omics matrix
#' @param type1 data1 omics type
#' @param type2 data2 omics type
#' @param method correlation method one of "spearman","person",defult is "spearman"
#' @return Correlation analysis of Multi-omics data.
#' @export
#' @import ggplot2 dplyr pheatmap Hmisc
#'
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


CorrPlotTotal <- function(data1 = data1,data2 = data2,type1="omics1",type2="omics2",
                          file_name=file_name,cutoff_p=0.05,labels_num = 1,
                          method="spearman",topn1=50,topn2=50,cutoff_cor=0.6,width = 10,height=8){
  if (identical(colnames(data1),colnames(data2))) {
    #type annotation
    data1$type = type1
    data2$type = type2
    data_com = rbind(data1,data2)
    xi = ncol(data_com)
    res= Hmisc::rcorr(as.matrix(t(data_com[,-xi])),type = method)
    heatmap_data = as.matrix(res$r[1:nrow(data1),(nrow(data1)+1):nrow(data_com)])
    cor_data_p = as.matrix(res$P[1:nrow(data1),(nrow(data1)+1):nrow(data_com)])

    plot_heatmap_data = heatmap_data[c(1:ifelse(nrow(data1)>topn1,topn1,nrow(data1))),
                                     c(1:ifelse(nrow(data2)>topn2,topn2,nrow(data2)))]
    plot_p_data = cor_data_p[c(1:ifelse(nrow(data1)>topn1,topn1,nrow(data1))),
                             c(1:ifelse(nrow(data2)>topn2,topn2,nrow(data2)))]

    display_numbers = matrix(ifelse(plot_p_data < 0.05,
                                    ifelse(plot_p_data < 0.01,
                                           ifelse(plot_p_data < 0.001,"***","**"),"*"),""),
                             nrow = nrow(plot_p_data))
    #plot heatmap

    bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))
    p = pheatmap::pheatmap(plot_heatmap_data,
                cluster_rows=T,cluster_cols=T,
                # annotation_row = annotation_row,
                # annotation_col = annotation_col,
                display_numbers = display_numbers,
                scale="none",
                #color=colorRampPalette(rev(c("#B1001F","orange","white","#1B52D9","#0F177A")))(2000),
                 color=c(colorRampPalette(colors = c("#0F177A","#1B52D9","white"))(round(length(bk)*(abs(min(plot_heatmap_data))/(max(plot_heatmap_data))))),
                         colorRampPalette(colors = c("white","orange","#B1001F"))(round(length(bk)))),
                show_rownames=T,show_colnames=T,
                cellwidth =400/ncol(plot_heatmap_data),
                cellheight=300/nrow(plot_heatmap_data),
                border_color =NA,number_color = "white",
                fontsize_col = 280/ncol(plot_heatmap_data),
                fontsize_row =min(300/nrow(plot_heatmap_data),12),
                fontsize = 14,
                legend_breaks = NA,fontsize_number = 280/ncol(plot_heatmap_data)+2,
                legend = TRUE,angle_col = 45,
                main="")

    ggplot2::ggsave(filename = paste0(file_name,"_相关性热图.png"),plot = p,
                    width = width,height = height,dpi = 900)
    ggplot2::ggsave(filename =  paste0(file_name,"_相关性热图.pdf"),plot = p,
                    width = width,height = height,dpi = 900)

    edge_data = MulOCor::cormatrix(res$r,res$P,data_com,data1)
    #nrow(edge_data)
    #cutoff
    edge_data_clear = edge_data[abs(edge_data$cor) > cutoff_cor & edge_data$p< cutoff_p,]
    #remove inter omics relation
    edge_data_clear = edge_data_clear[(edge_data_clear$rowtype != edge_data_clear$coltype),]

    write.csv(x = edge_data_clear,file =  paste0(file_name,"_相关性网络边文件.csv"),row.names = F)
    write.csv(x = plot_heatmap_data,file =  paste0(file_name,"_相关性数据.csv"),row.names = T)
    write.csv(x = plot_p_data,file =  paste0(file_name,"_相关性p值数据.csv"),row.names = T)



    datab = unique(data.frame(name = edge_data_clear$row,type = edge_data_clear$rowtype))

    # databb = openxlsx::read.xlsx("cor_result_data_cor_0.3.xlsx",1)
    # colnames(databb)[1:2] = c("name","type")
    # datab = left_join(data.frame(name = datab[,-2]),databb,by = "name")

    vertexdata = rbind(datab[,c(1,2)],
                       unique(data.frame(name = edge_data_clear$colum,type = edge_data_clear$coltype)))

    edge_data_clear$edgetype = ifelse(edge_data_clear$cor>0,"positive","negative")
    edge_data_clear$abscor = abs(edge_data_clear$cor)

    dnet <- c_net_from_edgelist(edgelist = edge_data_clear, direct = F,
                                vertex = vertexdata,e_type = "edgetype")

    plot(dnet, coors =nicely())
    # set vertex_class
    dnet <- c_net_set(dnet,vertex_class = c("type"))
    dnet <- c_net_set(dnet,edge_width  = c("abscor"))
    # dnet <- c_net_set(dnet,edge_class  = c("edgetype"))
    #dnet <- c_net_set(dnet,vertexdata,vertex.shapes= c("type"))


    #get_v(dnet) %>% head(5)

    dnet <- c_net_index(dnet)

    dnet <- c_net_set(dnet, vertex_size = "Degree")



    dnet <- c_net_set(dnet, vertex_class  = )

    module_detect(dnet,method = "cluster_walktra") -> dnet
    moduleclass =  get_v(dnet)[, c("name", "module")]


    #
    V(dnet)$shape = ifelse(get_v(dnet)$type== type1,"circle","square")
    # View(get_v(dnet))
    #
    # View((get_e(dnet)))


    c_net_layout(dnet, method = as_star()) -> coors

    pdf(file =  paste0(Sys.Date(),file_name,"_相关性网络图.pdf"),width = 10)
    plot(dnet , labels_num = length(dnet),edge.color= c("lightblue","red"),
         legend_cex = 1.3, coors = coors,group_legend_title = c("type"),
         edge.curved = 0.1,vertex_size_range = c(6,15))
    dev.off()

    pdf(file =  paste0(Sys.Date(),file_name,"_相关性网络图_module.pdf"),width = 10)
    c_net_plot(dnet,plot_module = T,mark_module = T, legend_cex = 1.3,coors = with_fr(),
               labels_num = labels_num*length(dnet), vertex.color = get_cols(11, "col1"), color_legend = F,
               vertex_size_range = c(7,12), vertex.label.color = "black",color_legend_title="x",
               edge_width_range = c(0.5, 3), edge.curved = 0.5,
               legend = T, legend_number = F, #group_legend_order = c("Microbiome", "Metabolome", "Transcriptome"),
               #group_legend_title = c("Phylum", "Metabolome", "Transcriptome"),
               edge_legend_title = "Correlation", edge_legend_order = c("positive", "negative"),
               size_legend = T, #size_legend_title = "Abundance",
               width_legend = T, width_legend_title = "abs(r)",
               lty_legend = F, lty_legend_title = "Omics relationship")
    dev.off()



    }
  else{
    print("两个组学间样本名称不一致或者样本顺序有误，请检查后修改再尝试")
  }

}

