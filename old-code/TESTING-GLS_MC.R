library(remotePARTS)

n.pix = 30865 # pixels in AK_ndvi_common-land.csv
parts = sample_partitions(npix = n.pix, npart = 4, partsize = 1000)
data.file = system.file("extdata", "AK_ndvi_common-land.csv",
                        package = "remotePARTS")

## Try 1 ----
df = data.table::fread(data.file)

part_data <- function(part.i, form, df, part_mat, locvars = c("lng", "lat")){
  prt = part_mat[, part.i]
  df = as.data.frame(df)
  df.prt = df[prt, ]
  X = model.matrix(cls.coef ~ 0 + land, data = df.prt)

  return(list(X = as.matrix(X),
              y = as.vector(df.prt$cls.coef),
              coords = as.matrix(df.prt[, locvars])))
}

part_data(1, cls.coef ~ 0 + land, df = df, part_mat = parts)


# GLS with 4 cores
GLS.part.mc = fitGLS.partition.mc(part_f = "part_data", dist_f = "dist_km",
                                  partsize = nrow(parts), npart = ncol(parts),
                                  V.meth = "exponential", spatcor = .5,
                                  # additional arguments passed to part_csv():
                                  ncores = 4,
                                  form = "cls.coef ~ 0 + land", df = df, part_mat = parts)

## Try 2 ----
# partDF = data.frame(part = as.vector(parts))
# df.sql = sqldf::read.csv.sql(file = data.file,
#                          sql = "select file.* from file join partDF on file.rowid = partDF.part",
#                          dbname = tempfile(),
#                          header = TRUE)
# #
# rel.est = data$cls.coef
# part.mat = parts
# mod.mat = model.matrix(cls.coef ~ land, data = df)
# location = data[, c("lng", "lat")]

## Function to get data partitions (will be passed to GLS_partition)
part_func <- function(part.i, part_mat, X.full, y.full, location.full){
  part = part_mat[, part.i]
  return(list(X = as.matrix(X.full[part, ]), y = y.full[part], coords = location.full[part, ]))
}

part_func(1, part_mat = parts, y.full = df$cls.coef, location.full = df[, c("lng", "lat")],
          X.full = model.matrix(cls.coef ~ 0 + land, data = df))

GLS.part.mc2 = fitGLS.partition.mc(part_f = "part_data", dist_f = "dist_km",
                                 partsize = nrow(parts), npart = ncol(parts),
                                 V.meth = "exponential", spatcor = .5,
                                 # additional arguments passed to part_csv():
                                 ncores = 4,
                                 df = df,
                                 part_mat = parts)

