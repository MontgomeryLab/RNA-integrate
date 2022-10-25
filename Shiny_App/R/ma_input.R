# MA Plot Function
ma_input = function(res,group1,group2) {
  return(plotMA(res,ylim=c(-4,4),main = paste0(group2,' vs ',group1)))
}