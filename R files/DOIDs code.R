
library(DOSE)

data = read.csv('../data/disease_signatures-v1.0.csv', sep=',', header = T, row.names = 1)

doids = unique(data$do_id)
fil = doids!=''
doids = doids[fil]


write.csv(DOSE::doSim(doids,doids),'../Diseases/doidstable.csv')



