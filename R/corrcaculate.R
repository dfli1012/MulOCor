#' @title correlation
#' @return Correlation analysis of Multi-omics data.
#' @export

#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
cormatrix=function(cor,p,data_com=data_com,data1 = data1){
  ut=upper.tri(cor)
  data.frame(row=rownames(cor)[row(cor)[ut]],
             colum=rownames(cor)[col(cor)[ut]],
             cor=cor[ut],
             p=p[ut],
             rowtype = data_com[rownames(cor)[row(cor)[ut]],ncol(data1)],
             coltype = data_com[rownames(cor)[col(cor)[ut]],ncol(data1)])

}
