#LIB <- new.env()
#LIB$library.dir <- library.dir
#l <- list.files(path=library.dir)
#l <- l[grep("\\.R",l)]
#tilde <- grep("\\.R~",l)
#if (length(tilde) > 0) l <- l[-tilde]
#l <- l[l != "library.R"]
#s <- sapply(l,FUN=function(x,library.dir, env) {
#  sys.source(file.path(library.dir,x), env)},
#            library.dir, env=LIB)
#rm(s,l,tilde)


l <- list.files(path=library.dir)
l <- l[grep("\\.R",l)]
tilde <- grep("\\.R~",l)
if (length(tilde) > 0) l <- l[-tilde]
l <- l[l != "library.R"]
s <- sapply(l,FUN=function(x,library.dir) {
  sys.source(file.path(library.dir,x), .GlobalEnv)},
            library.dir)
rm(s,l,tilde)
