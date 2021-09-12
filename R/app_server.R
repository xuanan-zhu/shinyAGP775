#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  # Your application server logic 
  
  vals <- reactiveValues(
    GAmax=NULL,
    LAmax=NULL,
    AGPmax=NULL
  )
  output$r1<-renderTable({
    if(is.null(vals$GAmax)){
      return()
    }
    else{return(vals$GAmax)}
    
  }) 
  output$r11<-renderTable({
    if(is.null(vals$LAmax)){
      return()
    }
    else{return(vals$LAmax)}
  }) 
  output$r1seq<-renderTable({
    if(input$submit1==0){
      return()
    }
    isolate({
      matrix<-matrix(c(input$m1,input$rm1,input$rm1,input$rm1,input$rm1,input$m1,input$rm1,input$rm1,input$rm1,input$rm1,input$m1,input$rm1,input$rm1,input$rm1,input$rm1,input$m1),nrow=4,ncol=4,dimnames=list(c("A","T","C","G"),c("A","T","C","G")))
      seqdata<-c(input$data1,input$data2)
      seqdata<-as.matrix(seqdata)
      seqdata1<-seqdata[1,]#提取序列1
      seqdata2<-seqdata[2,]#提取序列2
      seqdata1<-strsplit(seqdata1,"",fixed=T)
      seqdata2<-strsplit(seqdata2,"",fixed=T)
      zseqdata1<-as.character(unlist(seqdata1))
      zseqdata2<-as.character(unlist(seqdata2)) #zseqdata1和zseqdata2是转化成单个字符后的序列,去除列表
      #Needleman-Wunsch 算法
      gap=input$gap1#已知gap
      N<-length(zseqdata2)
      M<-length(zseqdata1)
      scorematrix<-matrix(0,N+1,M+1)#构造空矩阵，N+1行,M+1列
      rownames(scorematrix)<-c("-",zseqdata2)
      colnames(scorematrix)<-c("-",zseqdata1)
      dirmatrix<-scorematrix
      #计算第一行第一列
      scorematrix[1,1]=0
      for (i in 0:N+1)
        scorematrix[i,1]=gap*(i-1)
      for (j in 0:M+1)
        scorematrix[1,j]=gap*(j-1)
      #计算剩下的
      for (i in 1:N+1)
        for (j in 1:M+1)
        {
          scorematrix[i,j]=max(c(scorematrix[i-1,j-1]+matrix[rownames(scorematrix)[i],colnames(scorematrix)[j]],
                                 scorematrix[i-1,j]+gap,
                                 scorematrix[i,j-1]+gap))
          if(scorematrix[i,j]==scorematrix[i-1,j-1]+matrix[rownames(scorematrix)[i],colnames(scorematrix)[j]]){
            dirmatrix[i,j]<-"Q"
          }
          if(scorematrix[i,j]==scorematrix[i-1,j]+gap){
            dirmatrix[i,j]<-"W"
          }
          if(scorematrix[i,j]==scorematrix[i,j-1]+gap){
            dirmatrix[i,j]<-"A"
          }
        }
      vals$GAmax<-scorematrix[N+1,M+1]
      result1<-list()
      result2<-list()
      k<-0
      while (N!=0&M!=0)
      {
        ifelse(dirmatrix[N+1,M+1]=="Q",{
          k=k+1
          result1[k]<-colnames(scorematrix)[M+1]
          result2[k]<-rownames(scorematrix)[N+1]
          N=N-1
          M=M-1
        },{
          if(dirmatrix[N+1,M+1]=="A")
          {
            k=k+1
            result1[k]<-colnames(scorematrix)[M+1]
            result2[k]<-"-"
            M=M-1
          }
          else
          {
            k=k+1
            result1[k]<-"-"
            result2[k]<-rownames(scorematrix)[N+1]
            N=N-1
          }
        }
        )
      }
      result1<-rev(as.character(unlist(result1)))
      result2<-rev(as.character(unlist(result2)))
      frame<-cbind(result1,result2)
      colnames(frame)=c("seq1","seq2")
      print(frame)
    })
  })
  output$r11seq<-renderTable({
    if(input$submit3==0){
      return()
    }
    isolate({
      matrix<-matrix(c(input$m1,input$rm1,input$rm1,input$rm1,input$rm1,input$m1,input$rm1,input$rm1,input$rm1,input$rm1,input$m1,input$rm1,input$rm1,input$rm1,input$rm1,input$m1),nrow=4,ncol=4,dimnames=list(c("A","T","C","G"),c("A","T","C","G")))
      seqdata<-c(input$data1,input$data2)
      seqdata<-as.matrix(seqdata)
      seqdata1<-seqdata[1,]#提取序列1
      seqdata2<-seqdata[2,]#提取序列2
      seqdata1<-strsplit(seqdata1,"",fixed=T)
      seqdata2<-strsplit(seqdata2,"",fixed=T)
      zseqdata1<-as.character(unlist(seqdata1))
      zseqdata2<-as.character(unlist(seqdata2)) #zseqdata1和zseqdata2是转化成单个字符后的序列,去除列表
      #Needleman-Wunsch 算法
      gap=input$gap1#已知gap
      N<-length(zseqdata2)
      M<-length(zseqdata1)
      scorematrix<-matrix(0,N+1,M+1)#构造空矩阵，N+1行,M+1列
      rownames(scorematrix)<-c("-",zseqdata2)
      colnames(scorematrix)<-c("-",zseqdata1)
      dirmatrix<-scorematrix
      #计算第一行第一列
      scorematrix[1,1]=0
      for (i in 0:N+1)
        scorematrix[i,1]=0
      for (j in 0:M+1)
        scorematrix[1,j]=0
      #计算剩下的
      for (i in 2:(N+1))
        for (j in 2:(M+1))
        {
          scorematrix[i,j]=max(c(scorematrix[i-1,j-1]+matrix[rownames(scorematrix)[i],colnames(scorematrix)[j]],
                                 scorematrix[i-1,j]+gap,
                                 scorematrix[i,j-1]+gap))
          if(scorematrix[i,j]==scorematrix[i-1,j-1]+matrix[rownames(scorematrix)[i],colnames(scorematrix)[j]]){
            dirmatrix[i,j]<-"Q"
          }
          if(scorematrix[i,j]==scorematrix[i-1,j]+gap){
            dirmatrix[i,j]<-"W"
          }
          if(scorematrix[i,j]==scorematrix[i,j-1]+gap){
            dirmatrix[i,j]<-"A"
            
          }
          scorematrix[i,j]<-max(scorematrix[i,j],0)
        }
      maxt=0
      for (i in 2:(N+1))
        for (j in 2:(M+1))
        {
          if(scorematrix[i,j]>maxt){
            maxt<-scorematrix[i,j]
            maxi<-i
            maxj<-j
          }
        }
      result1<-list()
      result2<-list()
      k<-0
      N<-maxi-1
      M<-maxj-1
      vals$LAmax<-maxt
      while (N!=0&M!=0&&scorematrix[N+1,M+1]!=0)
      {
        ifelse(dirmatrix[N+1,M+1]=="Q",{
          k=k+1
          result1[k]<-colnames(scorematrix)[M+1]
          result2[k]<-rownames(scorematrix)[N+1]
          N=N-1
          M=M-1
        },{
          if(dirmatrix[N+1,M+1]=="A")
          {
            k=k+1
            result1[k]<-colnames(scorematrix)[M+1]
            result2[k]<-"-"
            M=M-1
          }
          else
          {
            k=k+1
            result1[k]<-"-"
            result2[k]<-rownames(scorematrix)[N+1]
            N=N-1
          }
        }
        )
      }
      result1<-rev(as.character(unlist(result1)))
      result2<-rev(as.character(unlist(result2)))
      frame<-cbind(result1,result2)
      colnames(frame)=c("seq1","seq2")
      print(frame)
    })
  })
  output$r2seq<-renderTable({
    if(input$submit2==0){
      
      return()
    }
    isolate({
      matrix<-matrix(c(input$m2,input$rm2,input$rm2,input$rm2,input$rm2,input$m2,input$rm2,input$rm2,input$rm2,input$rm2,input$m2,input$rm2,input$rm2,input$rm2,input$rm2,input$m2),nrow=4,ncol=4,dimnames=list(c("A","T","C","G"),c("A","T","C","G")))
      gap<-input$gap2
      gapp<-input$gapp2
      seqdata<-c(input$data1,input$data2)
      seqdata<-as.matrix(seqdata)
      seqdata1<-seqdata[1,]#提取序列1
      seqdata2<-seqdata[2,]#提取序列2
      seqdata1<-strsplit(seqdata1,"",fixed=T)
      seqdata2<-strsplit(seqdata2,"",fixed=T)
      zseqdata1<-as.character(unlist(seqdata1))
      zseqdata2<-as.character(unlist(seqdata2)) #zseqdata1和zseqdata2是转化成单个字符后的序列,去除列表
      N<-length(zseqdata2)
      M<-length(zseqdata1)
      babymatrix<-matrix(NA,N+1,M+1)#构造空矩阵，N+1行,M+1列
      babymatrix[1,1]=0
      rownames(babymatrix)<-c("-",zseqdata2)
      colnames(babymatrix)<-c("-",zseqdata1)
      Wmatrix<-babymatrix
      for (i in 2:(N+1))
        Wmatrix[i,1]=gap+gapp*(i-2)
      Vmatrix<-babymatrix
      for (i in 2:(M+1))
        Vmatrix[1,i]=gap+gapp*(i-2)
      
      Matrix<-babymatrix
      for (i in 2:(N+1))
        Matrix[i,1]=gap+gapp*(i-2)
      for (i in 2:(M+1))
        Matrix[1,i]=gap+gapp*(i-2)
      #构造方法矩阵
      gw<-Wmatrix
      for (i in 1:(N+1))
        gw[i,1]="F"
      gv<-Vmatrix
      for (i in 1:(M+1))
        gv[1,i]="F"
      #比对算法开始
      dirmatrix<-Matrix
      for (i in 2:(N+1))
        for (j in 2:(M+1))
        {
          Wmatrix[i,j]=max(Matrix[i-1,j]+gap,ifelse(is.na(Wmatrix[i-1,j]),Matrix[i-1,j]+gap,ifelse(gw[i-1,j]=="T",Wmatrix[i-1,j]+gapp,Matrix[i-1,j]+gap)))
          Vmatrix[i,j]=max(Matrix[i,j-1]+gap,ifelse(is.na(Vmatrix[i,j-1]),Matrix[i,j-1]+gap,ifelse(gv[i,j-1]=="T",Vmatrix[i,j-1]+gapp,Matrix[i,j-1]+gap)))
          t=Matrix[i-1,j-1]+matrix[rownames(Matrix)[i],colnames(Matrix)[j]]
          Matrix[i,j]=max(t,Wmatrix[i,j],Vmatrix[i,j])
          if(Matrix[i,j]==t){
            gv[i,j]="F"
            gw[i,j]="F"
            dirmatrix[i,j]="Q"}
          if(Matrix[i,j]==Wmatrix[i,j]){
            gw[i,j]="T"
            gv[i,j]="F"
            dirmatrix[i,j]="W"}
          if(Matrix[i,j]==Vmatrix[i,j]){
            gv[i,j]="T"
            gw[i,j]="F"
            dirmatrix[i,j]="A"}
        }
      vals$AGPmax<-Matrix[N+1,M+1]
      result1<-list()
      result2<-list()
      k<-0
      while (N!=0&M!=0)
      {
        ifelse(dirmatrix[N+1,M+1]=="Q",{
          k=k+1
          result1[k]<-colnames(dirmatrix)[M+1]
          result2[k]<-rownames(dirmatrix)[N+1]
          N=N-1
          M=M-1
        },{
          if(dirmatrix[N+1,M+1]=="A")
          {
            k=k+1
            result1[k]<-colnames(dirmatrix)[M+1]
            result2[k]<-"-"
            M=M-1
          }
          else
          {
            k=k+1
            result1[k]<-"-"
            result2[k]<-rownames(dirmatrix)[N+1]
            N=N-1
          }
        }
        )
      }
      result1<-rev(as.character(unlist(result1)))
      result2<-rev(as.character(unlist(result2)))
      frame<-cbind(result1,result2)
      colnames(frame)=c("seq1","seq2")
      print(frame)
    })
  })
  output$r2<-renderTable({
    if(is.null(vals$AGPmax)){
      return()
    }
    else{return(vals$AGPmax)}
  })
}
