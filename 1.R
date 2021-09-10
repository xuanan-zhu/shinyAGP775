matrix<-matrix(c(10,-1,-3,-4,-1,7,-5,-3,-3,-5,9,0,-4,-3,0,8),nrow=4,ncol=4,dimnames=list(c("A","T","C","G"),c("A","T","C","G")))
seqdata<-c("ACGTC","AATC")
seqdata<-as.matrix(seqdata)
seqdata1<-seqdata[1,]#提取序列1
seqdata2<-seqdata[2,]#提取序列2
seqdata1<-strsplit(seqdata1,"",fixed=T)
seqdata2<-strsplit(seqdata2,"",fixed=T)
zseqdata1<-as.character(unlist(seqdata1))
zseqdata2<-as.character(unlist(seqdata2)) #zseqdata1和zseqdata2是转化成单个字符后的序列,去除列表
#Needleman-Wunsch 算法
gap=-5#已知gap
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
