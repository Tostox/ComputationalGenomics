DBs <- matrix(c(
    "dlbcl_new_out.rma.10Kmad.res",  "dlbcl_new_out.cls",        "dlbcl",        "-unlog",
    "duke_breast.3Kmad.res",         "duke_breast.cls",          "duke",         "", 
    "duke_breast.3Kmad.res",         "duke_breast.er.cls",       "duke.er",      "",
    "leukemia.ms.5Kmad.res",         "leukemia.cls",             "leuk",         "",
    "lung.pnas.rma.nomiss.5Kmad.gct","lung.pnas.nomiss.mut.cls", "lung",         "-unlog",
    "prostate.rma.10Kmad.res",       "prostate.cls",             "prostate",     "confound=prostate.study.cls -unlog",
    "rosetta_breast_out.10Kmad.res", "rosetta_breast_out.cls",   "rosetta",      "",
    "rosetta_breast_out.10Kmad.res", "rosetta_breast_out.er.cls","rosetta.er",   "",
    "RosettaNEJM2002.10Kmad.res",    "RosettaNEJM2002.cls",      "rosetta.nejm", "",
    "all_aml_mll.rma.5Kmad.gct",     "all_aml_mll.cls",          "AllAmlMLL",    "-unlog",
    "cns.rma.morpho.3Kmad2.gct",     "cns.morphology.cls",       "cns.morph",    "-unlog",
    "dmap.10Kmad.pop13.gct",         "dmap.pop13.cls",           "dmap13",       "",
    "dlbcl_media_new.rma.5Kmad.res","dlbcl_media_new.cls",      "mlbcl",        "-unlog",
    "EGFR.pnas.rma2.nocm.mut.5Kmad.gct","EGFR.pnas.nocm.mut.cls",      "egfr",  ""),
    nrow=14, ncol=4, byrow=T, dimnames=list(NULL,c("DB","CLS","ID","PARAM")))

stub <- "xchip"
missing <- rep(F,nrow(DBs))

for ( i in 1:nrow(DBs) ) {
  missing[i] <- file.access(paste(stub,DBs[i,"DB"],sep="/"))[1]!=0
}
if (any(missing)) {
  cat( "Files missing:\n", paste(DBs[missing,"DB"],sep="\n"), "\n" )
}
for ( i in 1:nrow(DBs) ) {
  missing[i] <- file.access(paste(stub,DBs[i,"CLS"],sep="/"))[1]!=0
}
if (any(missing)) {
  cat( "Files missing:\n", paste(DBs[missing,"CLS"],sep="\n"), "\n" )
}

