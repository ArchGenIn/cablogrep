#!/usr/bin/env Rscript
d = commandArgs(trailingOnly = T)
input1 = d[1]
ref = read.table(d[2], stringsAsFactors = F, header = T)                        #read reftab
hgmap = read.table(d[3], stringsAsFactors = F, header = F, sep = ",", quote = "\"")           #read hg map
dname = dirname(input1)                                   #get sample name
input1 = read.table(input1, stringsAsFactors = F, header = F, skip = 5)         #read MUMMER input, exclude first 5 rows
slist = unique(input1[,ncol(input1)])
for(i in 1:length(slist)){
  basename = slist[i]
  input = input1[which(slist[i] == input1[,ncol(input1)]),c(1:3)]
  input$V2 = toupper(input$V2)                                                    #if lowercase, don't
  input$V3 = toupper(input$V3)                                                    #if lowercase, don't
  rl = ref[ref$Position %in% input$V1,]                                           #find informative positions
  #filter out listed positions with unlisted mutation type
  tr1 = c("C", "T")
  tr2 = c("A", "G")
  finfin = as.data.frame(matrix(NA, 1, 5))
  colnames(finfin) = c("Sample", "Best_Hit", "Quality_1", "2nd_Best", "Quality_2")
  finfin$Sample = basename
  if(nrow(rl) > 0){
    for(i in 1:nrow(rl)){
      pair = unname(unlist(input[rl$Position[i] == input$V1,c(2,3)]))
      pair = na.omit(pair)
      if(rl$Mut_type[i] == "X"){                                                    #if not transition
        if(!identical(tr1, sort(pair)) & !identical(tr2, sort(pair))){
          input[which(rl$Position[i] == input$V1),] = NA
        }
      } else if(rl$Mut_type[i] == "D" & pair[2] != "."){                            #if not deleterious
        input[which(rl$Position[i] == input$V1),] = NA
      } else if(rl$Mut_type[i] == "A" | rl$Mut_type[i] == "G" | rl$Mut_type[i] == "T" | rl$Mut_type[i] == "C"){ #if not transversion
        if(identical(tr1, sort(pair)) | identical(tr2, sort(pair)) | any(pair == ".")){
          input[which(rl$Position[i] == input$V1),] = NA
        }
      } else if(rl$Mut_type[i] == "I"){                                             #if not or invalid insertion
        if(pair[1] != "." | pair[2] != rl$Ins_type[i]){
          input[which(rl$Position[i] == input$V1),] = NA
        }
      }
    }
    input = na.omit(input)
    rl = ref[ref$Position %in% input$V1,]
    #get node list seq
    nlist = unique(rl$Haplogroup)
    hmlist = as.list(as.data.frame(t(hgmap)))
    hit = list()
    for(i in 1:length(hmlist)){
      tsub = hmlist[[i]]
      tsub = tsub[tsub != '']
      if(any(tsub[length(tsub)] == nlist)){
        tsub = list(tsub)
        hit = c(hit, tsub)
      }
    }
    #find best sequence hit
    if(length(hit) == 0){
      finfin$Best_Hit = "A1"
      finfin$Quality_1 = "NA"
      finfin$`2nd_Best`= "NA"
      finfin$Quality_2 = "NA"
      print(finfin)
    } else {
      fin = as.data.frame(matrix(NA, length(hit), 2))
      colnames(fin) = c("mtHG", "Q")
      for(i in 1:length(hit)){
        tsub = hit[[i]]
        ttab = as.data.frame(matrix(NA, 3, length(tsub)))
        ttab[1,] = tsub
        for(j in 1:ncol(ttab)){
          ez = ref[ttab[1,j] == ref$Haplogroup,]
          ttab[2,j] = nrow(ez)
          az = rl$Mut_type[ttab[1,j] == rl$Haplogroup]
          ttab[3,j] = length(az) - length(az[az == "back"])
        }
        fin$mtHG[i] = tsub[length(tsub)]
        fin$Q[i] = ((sum(as.numeric(ttab[3,])) / sum(as.numeric(ttab[2,]))) + (sum(as.numeric(ttab[3,])) / nrow(rl))) / 2
      }
      fin2 = fin[order(fin$Q, decreasing = T),]
      fin2$Q = round(fin2$Q, digits = 3)
      write.table(fin2, file = paste(dname, "/", basename, ".HgCandidates", sep = ""), row.names = F, quote = F, sep = "\t")
      finfin$Best_Hit = fin2$mtHG[1]
      finfin$Quality_1 = fin2$Q[1]
      finfin$`2nd_Best`= fin2$mtHG[2]
      finfin$Quality_2 = fin2$Q[2]
      print(finfin)
    }
  } else {
    finfin$Best_Hit = "A1"
    finfin$Quality_1 = "NA"
    finfin$`2nd_Best`= "NA"
    finfin$Quality_2 = "NA"
    print(finfin)
  }
  write.table(finfin, file = paste(dname, "/", basename, ".HG", sep = ""), row.names = F, col.names = F, sep = "\t", quote = F)
}
