lapply(0:600,
       function(x) {
           print(x)
           load(paste("nmf_result",x,sep=""))
           measures <- nmf_result$measures
           save(measures,file=paste("measures",x,sep=""))
       }
)
