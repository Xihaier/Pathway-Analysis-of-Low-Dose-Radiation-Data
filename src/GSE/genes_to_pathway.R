
# before running this, make sure you have run 
#  library( gage )
#  data( kegg.gs )
#
# Example:
#
#    genes_to_kegg_pathway
#
#     > genes_to_kegg_pathway( c(335, 336) ) 
#     $`335`
#     [1] "hsa03320 PPAR signaling pathway"          
#     [2] "hsa04975 Fat digestion and absorption"    
#     [3] "hsa04977 Vitamin digestion and absorption"
#     
#     $`336`
#     [1] "hsa03320 PPAR signaling pathway"
#


genes_to_kegg_pathway <- function( genes )
   {
    ngenes <- as.numeric( genes )

    result <- list()    
    nkpaths <- length( kegg.gs )    
    for ( g in ngenes )
       {
        paths <- c()
        for ( j in 1:nkpaths )
           {
            if ( g %in% as.numeric( kegg.gs[[j]] ) )
            paths <- c( paths, names( kegg.gs[j] ) )
           }
        # cat( g, ":", paths, "\n" );
        result <- c( result, list(paths) )
       }
    names( result ) <- genes
    result
   }




