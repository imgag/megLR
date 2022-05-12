# Custom ggplot2 theme
theme_cas <- function(){
    # Change main theme
    theme_bw() %+replace% 
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
}

my_palette <- function(...) {  
    # Colour Brewer Set1
    colours <- c(
        "#E41A1C",
        "#377EB8",
        "#4DAF4A",
        "#984EA3",
        "#FF7F00",
        "#FFFF33",
        "#A65628",
        "#F781BF",
        "#999999"
        )
        
    function(n) {
        if (n == 0) {
            stop("Must request at least one colour from a hue palette.", 
                call. = FALSE)
        } else if (n<10) {
            sel_colours = colours[1:n]
        } else {
            sel_colours = colorRampPalette(colours)(n)
        }
        sel_colours
    }
    
}

scale_fill_custom <- function(..., aesthetics = "fill") {
    
    discrete_scale(aesthetics, "custom", my_palette(), ...)
}


mytheme <- list(
    theme_cas(),
#    scale_fill_custom(),
    scale_fill_distiller()
)
