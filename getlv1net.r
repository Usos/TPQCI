require(Matrix)
require(dplyr)
require(igraph)
require(foreach)


#Calculate corrlation between gene A's copy number and gene B's RNA expression
#... is for corrlation calculate method - pearson, spearman and kendall
getLV1cor<-function(cnv,rna,ppn,start=1,n=0,...)
{
  lv1.adj<-as_adjacency_matrix(ppn)+Diagonal(length(V(ppn)));
  lv1<-graph.adjacency(lv1.adj,mode = "directed");
  nodes<-V(lv1)$name;
  nodesrank<-rank(nodes);
  names(nodesrank)<-nodes;
  
  cnv<- arrange(cnv,Gene) %>% select(-Gene);
  rna<- arrange(rna,Gene) %>% select(-Gene);
  
  cnv<-t(cnv);
  rna<-t(rna);
  
  lv1.edges<-as_edgelist(lv1);
  nEdge<-dim(lv1.edges)[1]
  if(n!=0 && (start+n-1)<nEdge)
    nEdge<-start+n-1
  
  cors<-foreach(i=start:nEdge,.combine = c) %do%
    cor(cnv[,nodesrank[lv1.edges[i,1]]],
        rna[,nodesrank[lv1.edges[i,2]]],...)

  colnames(lv1.edges)<-c("node_c","node_r");
  lv1.edges<-as.data.frame(lv1.edges);
  lv1.edges$node_c<-as.character(lv1.edges$node_c)
  lv1.edges$node_r<-as.character(lv1.edges$node_r)
  cors<-mutate(lv1.edges[start:nEdge,],cor=cors);
  
  #E(lv1)$weight<-cors;
  
  
  #return(lv1);
  return(cors);
  
}

rnaLog2Trans<-function(rna)
{
  genes<-rna$Gene;
  rna<-log2(rna %>% select(-Gene)+1) %>% mutate(Gene=genes) %>% select(Gene,everything())
  
  return(rna)
}

#get percentage of copy number abnormal in all samples
getCNVMass<-function(cnv,threshold=0.3)
{
  mass<-apply(select(cnv,-Gene),1,function(x){sum(abs(x)>threshold)/length(x)})
  names(mass)<-cnv$Gene
  return(mass)
}


#Calculate sigma value of selected pcc threshold
getSigma<-function(pcc.threshold=0.3, dec1=T)
{
  if(dec1)
    return(sqrt(2)/3*(1/pcc.threshold-1))
  return(sqrt(2)/(3*pcc.threshold));
}

#Calculate topology potential
calTopoPot<-function(mass,weight,sigma=0.3, dec1=T)
{
  sigma<-getSigma(sigma)
  topot<-weight %>%
    mutate(dist=1/cor,dist1=1/cor-1,mass=mass[node_c]) %>%
    mutate(field=mass*exp(-(dist/sigma)^2),field1=mass*exp(-(dist1/sigma)^2)) %>%
    group_by(node_c) %>%
    summarise(topot=sum(field),topot1=sum(field1)) 
  if(dec1)
    topot<-topot %>%  
      arrange(desc(topot1)) %>%
      select(node_c,topot=topot1)
  else
    topot<-topot %>%  
      arrange(desc(topot)) %>%
      select(node_c,topot)
      
  topot<- topot %>% mutate(abcR=mass[node_c])
  return(topot)
}

