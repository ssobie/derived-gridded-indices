##Script to clean up temporary directories

tmp.dir <- '/local_temp/ssobie/prism/'


print(list.files(path=tmp.dir,recursive=T))

clean.up <- list.files(path=tmp.dir,recursive=T,full.name=T)
file.remove(clean.up)

