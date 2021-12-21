
library(GEOquery)
library(Biobase)
library(limma)
library(stringr)



  experiemtn_id = commandArgs(TRUE)[1]
  geo_id = commandArgs(TRUE)[2]
  pert_ids = commandArgs(TRUE)[3]
  ctrl_ids = commandArgs(TRUE)[4]
  
  experiemtn_id = str_replace(experiemtn_id, ':', '-')

  ctrl_ids = strsplit(ctrl_ids, '|', fixed = TRUE)[[1]]
  pert_ids = strsplit(pert_ids, '|', fixed = TRUE)[[1]]
  
  
  data = getGEO(geo_id, GSEMatrix=TRUE, AnnotGPL = TRUE)
  data = data[[1]]
  
  
  
  expression = exprs(data)
  dim(expression)
  expression = as.data.frame(expression)
  dim(expression)
  
  ### remove negative values
  expression[expression<0] = 0
  expression[is.na(expression)] <- 0
  dim(expression)
  
  if (max(expression)>100){ 
    expression = log2(expression+1)
      }
  
  ### get gene names
  annot = fData(data)$`Gene symbol`
  # apply function (x: matrix or array, Margin 1 or 2: 1 for rows 2 for columns, FUN: function we want to apply)
  expression$MEAN = apply(expression, 1, mean)
  expression$GENE = annot
  
  ### remove genes without gene symbol
  fil = expression$GENE != ''
  expression = expression[fil,]
  
  # sorting in descending order >> (decreasing = true)
  expression = expression[order(expression$MEAN, decreasing = TRUE),]
  ### remove duplicates
  expression = expression[!duplicated(expression$GENE),]
  
  rownames(expression) = expression$GENE
  
  ### removed unnecesarry columns
  n = dim(expression)[2]
  expression = expression[, -c(n, n-1)]
  dim(expression)
  expression = expression[,c(ctrl_ids,pert_ids)]
  colnames(expression)
  ctrl_ids
  
  ### fit limma model
  indicator = rep(0, length(ctrl_ids) + length(pert_ids))
  indicator
  names(indicator) = c(ctrl_ids,pert_ids)
  indicator[pert_ids] = 1
  design=model.matrix(~1 + indicator)
  rownames(design) = names(indicator)
  rownames(design)
  
  dim(design)
  dim (expression)
  
  expression = expression[, rownames(design)]

  fit <- lmFit(expression, design)
  fit = eBayes(fit)
  
  results = topTable(fit, coef = 'indicator', adjust="BH", number = 1000000)
  
  write.csv(results, paste0('../Diseases/', experiemtn_id,'.csv')) 


