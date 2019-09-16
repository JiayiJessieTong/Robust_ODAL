library(plotly) 

# one variable
one_var <- function(Result_setting, setting){
  R1_plot = cbind(Result_setting[,3],Result_setting[,5],Result_setting[,7], 
                  Result_setting[,9], Result_setting[,11], Result_setting[,13],
                  Result_setting[,15],Result_setting[,17],Result_setting[,22],nn1)
  R1_plot[,1:9] = R1_plot[,1:9] - matrix(beta_true[2], nrow = 10, ncol = 9)
  temp = rep(0,10)
  names(temp) =c("ODAL1 mean #1-9","ODAL1 median #1-9","ODAL1 mean","ODAL1 median", 
                 "ODAL2 mean mean", "ODAL2 mean median", "ODAL2 median mean", "ODAL2 median median",
                 "Meta","n")
  temp = cbind(t(t(temp)),t(R1_plot))
  R1_plot = t(temp)
  R1_plot = R1_plot[2:11,]
  DR1_PE <- as.data.frame(R1_plot[,c(1:10)])
  
  # plot
  p <- plot_ly(DR1_PE, 
               # odal 1 mean # 1-9
               x = DR1_PE$n, y = DR1_PE$`ODAL1 mean #1-9`, name = 'ODAL1 mean #1-9', type = 'scatter', mode = 'lines+markers'
               ,line = list(width = 3), symbols = c(16,15)) %>%
    #odal 1 mean
    add_trace(y = DR1_PE$`ODAL1 mean`, name = 'ODAL1 mean', mode = 'lines+markers',
              line = list( width = 3)) %>%
    #odal 2 mean mean
    add_trace(y = DR1_PE$`ODAL2 mean mean`, name = 'ODAL2 mean mean', mode = 'lines+markers',
              line = list(width = 3)) %>%
    #odal 2 mean median
    add_trace(y = DR1_PE$`ODAL2 mean median`, name = 'ODAL2 mean median', mode = 'lines+markers',
              line = list( width = 3)) %>%
    
    #odal 1 median # 1-9
    add_trace(y = DR1_PE$`ODAL1 median #1-9`, name = 'ODAL1 median #1-9', mode = 'lines+markers',
              line = list(width = 3)) %>%
    #odal 1 median
    add_trace(y = DR1_PE$`ODAL1 median`, name = 'ODAL1 median', mode = 'lines+markers',
              line = list(width = 3)) %>%
    #odal 2 median mean
    add_trace(y = DR1_PE$`ODAL2 median mean`, name = 'ODAL2 median mean', mode = 'lines+markers',
              line = list(width = 3)) %>%
    #odal 2 median median
    add_trace(y = DR1_PE$`ODAL2 median median`, name = 'ODAL2 median median', mode = 'lines+markers',
              line = list(width = 3)) %>%
    
    #Meta
    add_trace(y = DR1_PE$`Meta`, name = 'Meta', mode = 'lines+markers',
              line = list(width = 3, dash = "dash")) %>%
    
    layout(autosize = F, width = 1200, height = 700) %>%
    layout(title = paste(setting,"Bias"))
  return(p)
}

p1 <- one_var(Result_setting1,setting = "1/10 Hetero, setting1")
p2 <- one_var(Result_setting2,setting = "2/10 Hetero, setting2")
p3 <- one_var(Result_setting3,setting = "1/10 Hetero, setting3")
p4 <- one_var(Result_setting4,setting = "2/10 Hetero, setting4")
htmlwidgets::saveWidget(as_widget(p1, p2, p3, p4), "index.html")

