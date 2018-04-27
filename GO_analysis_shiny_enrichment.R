library(shiny)
library(ggplot2)
library(grid)
library(shiny)
library(gridGraphics)
library(rsconnect)
library(gdata)
make_custom_barplot_for_doris = function(dat, place, set_limit=2, color="dodgerblue", logscale=F, set_size=c(0,1500), xcat= "lfe",remove_underrep = F){
  print(paste("set_limit ",set_limit))
  names(dat) = c("GO_term","total", "uploaded", "expected", "over_under", "fold_enrichment", "raw_P_value", "FDR")
  dat$logFDR = -log(as.numeric(dat$FDR),10)   # this is how it was presented in https://doi.org/10.1016/j.molmet.2015.11.001 the example paper doris gave me
  dat$log10_fold_enrichment = log(as.numeric(dat$fold_enrichment), 10)
  dat$GO_term = substr(x=as.character(dat$GO_term), start = 1, stop = nchar(as.character(dat$GO_term))-13)
  dat$GO_term = gsub(pattern = " \\(", replacement = "", x=as.character(dat$GO_term))
  dat=dat[order(dat$logFDR, decreasing = F),]
  dat=dat[dat$logFDR>set_limit,]
  dat=dat[dat$total<set_size[2] & dat$total>set_size[1],] # if the GO term is too basic
  #dat=drop.levels(dat) #somehow the "up dataset has factors and couses problems downstram. I am going to test that
  print(dat$fold_enrichment)
  if (remove_underrep==T){
    dat = dat[unlist(dat$fold_enrichment, use.names = F) > 1,]
  }
  print(paste("remove_underrep",remove_underrep,"xcat:",xcat,"set_size:",set_size,"logscale:",logscale,"set_limit:",set_limit,"place:",place, sep = ' '))
  dat$GO_term = paste(dat$GO_term, " (", dat$uploaded, "/", dat$total, ")")
  #dat$group = c("log10(fold enrichment", "-log(FDR)")
  title="GO term (genes / full set)"
  
  require(ggplot2)
  title="GO term (genes / full set)"
  #if(logscale==T){
  print(xcat)
  if(xcat=="lfe"){                #######################################################################################################
    #print(paste("dataset", dset))
    dat$GO_term = factor(dat$GO_term, levels = dat$GO_term[order(dat$fold_enrichment, decreasing = F)])      ##### log10(fold enrichment)
    print("B")
    p=ggplot(dat = dat, aes(x = GO_term, width = 0.8, y = log10_fold_enrichment)) +
      ggtitle(title) +
      labs(x="GO terms",y="log10(fold enrichment)")+
      theme(plot.title = element_text(hjust = place),
            axis.title.x = element_text(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain")) +
      geom_bar(stat = "identity", color = "gray49", fill = color) +
      #scale_fill_manual(values=c("+"="dodgerblue","-"= "darkorchid4")) +
      #theme(legend.title=element_blank()) +
      #labs(title = ) +
      coord_flip()
    return(p)
  }
  if(xcat=="lfdr"){             ############################################################################################ -log10FDR
    dat$GO_term = factor(dat$GO_term, levels = dat$GO_term[order(dat$logFDR, decreasing = F)])
    #dat$fold_enrichment = log(as.numeric(dat$fold_enrichment), 2)
    p=ggplot(dat = dat, aes(x = GO_term, width = 0.8, y = logFDR)) +
      ggtitle(title) +
      labs(x="GO terms",y="-log10(FDR)")+
      theme(plot.title = element_text(hjust = place),
            axis.title.x = element_text(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain")) +
      geom_bar(stat = "identity", color = "gray49", fill = color) +
      #scale_fill_manual(values=c("+"="dodgerblue","-"= "darkorchid4")) +
      #theme(legend.title=element_blank()) +
      #labs(title = ) +
      coord_flip()
    return(p)
  }
}
down=read.table("down_005_Stat5_inclusive.txt", header = T, sep ="\t", skip = 10 )
#names(down) = c("GO_term","total", "uploaded", "expected", "over_under", "fold_enrichment", "raw_P_value", "FDR")
#d = make_custom_barplot_for_doris(dat=down, place = -5, set_limit = 0.001, color = "dodgerblue", logscale = T, xcat="lfe", remove_underrep = T);d

################################################## IMPORTANT NOTE ##################################################                                  <<<<<<<<<<<<<<<<<<<#######!!!!!!!!!!!!!!!!!!!
#     the UP file had one fold_change value of "< 0.01". it generated a character factor for dat$fold_change because of the '<' unlike for the DOWN file which is numeric. 
#     this bug was difficult to find! removed manually in file. load file to variable
up=read.table("up_005_Stat5_inclusive.txt", header = T, sep ="\t", skip = 10 )
#names(up) = c("GO_term","total", "uploaded", "expected", "over_under", "fold_enrichment", "raw_P_value", "FDR")
#u = make_custom_barplot_for_doris(dat=up, place = 3, set_limit = 0.001, color="tomato2", logscale = F, set_size = c(0,1500), remove_underrep = T, xcat = "lfdr");u

# Define UI for application that draws a hitogram with GO terms
#################################################################################################################################

ui <- pageWithSidebar(
  headerPanel(textOutput("transformedValue")),
  sidebarPanel(
    selectInput(inputId = "dset", label = "Select set of genes", choices = c("Downregulated" = "down", "Upregulated" = "up")),
    selectInput(inputId = "xcat", label = "Select x-axis category", choices = c("log10(fold enrichment)" = "lfe", "-log10(FDR)" = "lfdr")),
    sliderInput(inputId = "logFDR", label = "set significance cutoff based on -log10(FDR)", min = 1, max = 6,  value = 2, step = 0.1),
    sliderInput(inputId = "set_size", label = "lower and upper limits of GO term set size", min = 0, max = 4000,  value = c(0,1500), step = 5),
    checkboxInput("checkbox_underrepresented", label = "remove underrepresented GO terms", value = TRUE)
  ),
  mainPanel(
    # Use imageOutput to place the image on the page
    imageOutput("myImage")
    
  )
)

server <- function(input, output, session) {
  output$transformedValue = renderText({paste("GO terms with FDR < ", round(10^-input$logFDR, 6))})
  output$myImage <- renderImage({
    # A temp file to save the output.
    # This file will be removed later by renderImage
    outfile <- tempfile(fileext = '.png')
    print(paste("input$dset ", input$dset))
    # Generate the PNG
    #png(outfile, width = 500, height = ifelse(test = input$logFDR >=1, (1/input$logFDR)*800, 800))
    if(input$dset == "up"){
      png(outfile, width = 600, height = ifelse(test = input$logFDR >=1, (1/input$logFDR)*900, 900))
      grid.draw(make_custom_barplot_for_doris(dat=up, place = -5, set_limit = input$logFDR , 
                                              color="tomato2", logscale = T, set_size = c(input$set_size[1], 
                                                                                          input$set_size[2]), xcat  = input$xcat, 
                                              remove_underrep = input$checkbox_underrepresented))
      dev.off()
    }
    if(input$dset == "down"){
      png(outfile, width = 600, height = ifelse(test = input$logFDR >=1, (1/input$logFDR)*900, 900))
      grid.draw(make_custom_barplot_for_doris(dat=down, place = -5, set_limit = input$logFDR ,
                                              color="dodgerblue", logscale = T, set_size = c(input$set_size[1], 
                                                                                             input$set_size[2]), xcat = input$xcat,
                                              remove_underrep = input$checkbox_underrepresented))
      dev.off()
    }
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 600,
         height = ifelse(test = input$logFDR >=1, (1/input$logFDR)*900, 900),
         alt = "This is alternate text")
  }, deleteFile = TRUE)
}
#dev.off()
shinyApp(ui = ui, server = server)