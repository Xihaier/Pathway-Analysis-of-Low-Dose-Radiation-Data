filename=system.file("extdata/gse16873.demo", package = "gage")
demo.data=readExpData(filename, row.names=1)

read.tcsv = function(file, header=TRUE, sep=",", ...) {
    
    n = max(count.fields(file, sep=sep), na.rm=TRUE)
    x = readLines(file)
    
    .splitvar = function(x, sep, n) {
        var = unlist(strsplit(x, split=sep))
        length(var) = n
        return(var)
    }
    
    x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
    x = apply(x, 1, paste, collapse=sep) 
    out = read.csv(text=x, sep=sep, header=header, ...)
    return(out)
    
}

df <- read.tcsv(file = '../GSE6978/data/GSE6978_gsnz.csv')

str(df)
