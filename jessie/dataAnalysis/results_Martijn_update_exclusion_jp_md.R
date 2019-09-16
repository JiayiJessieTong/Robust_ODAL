# this is the code to read the real data analysis results from Martijn
N = 64222+1976+59861+69164+62348

setwd("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/dataAnalysis")

# output pooled with all five sites
output1 = readRDS("/Users/jiayito/Dropbox/2_Rui_and_Jessie/Summer_2019/4_ODAL_median/real_data_results/from_Martijn/output_ODAL_median_AMI_ccae.rds")

# output without japan (can use the pooled as gold standard, since japan is an outlier)
output_wo_jp = readRDS("/Users/jiayito/Dropbox/2_Rui_and_Jessie/Summer_2019/4_ODAL_median/real_data_results/from_Martijn_3_correct/output_ODAL_median_AMI_ccae_noJmdc.rds")

# output without japan nor mdcr (can use the pooled as gold standard, since japan and mdcr can both be regarded as outliers)
output_wo_jp_md = readRDS("/Users/jiayito/Dropbox/2_Rui_and_Jessie/Summer_2019/4_ODAL_median/real_data_results/from_Martijn_3_correct/output_ODAL_median_AMI_ccae_noJmdcNorMdcr.rds")



##########################################
odal_result_func <- function(est, hessian, num){
  # this function is used to calculate the variance of ODAL method point estimate with hessian matrix returned from each site
  # there are two methods can be used to do this:
  # method 1: use local hessian hessian matrix only (inverse of negative hessian matrix)
  # method 2: use the mean of all hessian matrix
  # 
  # Arguments:
  #   est: point estimate
  #   hessian: hessian matrix of this point estiamte
  # Output:
  #   se and CI of two methods
  
  ## beta_local 
  beta_odal = est
  
  # # variance method 1: use local 
  # beta_odal_var_1 = diag(solve(matrix(-hessian[1,], ncol = num, nrow = num)))/N
  # beta_odal_se_1 = sqrt(beta_odal_var_1)
  # beta_odal_CI_l_1 = beta_odal - qnorm(0.975)*beta_odal_se_1
  # beta_odal_CI_u_1 = beta_odal + qnorm(0.975)*beta_odal_se_1
  
  # varaince method 2: use all sites
  # variance method 1: use local 
  beta_odal_var_2 = diag(solve(matrix(-apply(hessian, 2, mean), ncol = num, nrow = num)))/N
  beta_odal_se_2 = sqrt(beta_odal_var_2)
  beta_odal_CI_l_2 = beta_odal - qnorm(0.975)*beta_odal_se_2
  beta_odal_CI_u_2 = beta_odal + qnorm(0.975)*beta_odal_se_2
  
  # return the results
  return(list(se_2 = beta_odal_se_2,
              lower_2 = beta_odal_CI_l_2,
              upper_2 = beta_odal_CI_u_2))
}


#####################################
clean_output <- function(output,num){
  
  # beta_pooled
  beta_pooled = output$beta_pooled
  beta_pooled_var = diag(matrix(output$var_matrix_pooled, ncol = num, nrow = num))
  beta_pooled_se = sqrt(beta_pooled_var)
  beta_pooled_CI_l = beta_pooled - qnorm(0.975)*beta_pooled_se
  beta_pooled_CI_u = beta_pooled + qnorm(0.975)*beta_pooled_se
  
  ## beta_meta
  beta_meta = output$beta_meta
  beta_meta_se = sqrt(output$var_meta)
  beta_meta_CI_l = beta_meta - qnorm(0.975)*beta_meta_se
  beta_meta_CI_u = beta_meta + qnorm(0.975)*beta_meta_se
  
  ########### META ###################
  ## beta_ODAL1_meta_mean
  ODAL1_meta_mean = output$ODAL1_meta_mean
  tmp = odal_result_func(output$ODAL1_meta_mean, output$hessian_ODAL1_meta_mean,num)
  # method 1
  ODAL1_meta_mean_se_1 = tmp$se_1
  ODAL1_meta_mean_CI_l_1 = tmp$lower_1
  ODAL1_meta_mean_CI_u_1 = tmp$upper_1
  # method2
  ODAL1_meta_mean_se_2 = tmp$se_2
  ODAL1_meta_mean_CI_l_2 = tmp$lower_2
  ODAL1_meta_mean_CI_u_2 = tmp$upper_2
  
  ## beta_ODAL1_meta_median
  ODAL1_meta_median = output$ODAL1_meta_median
  tmp = odal_result_func(output$ODAL1_meta_median, output$hessian_ODAL1_meta_median,num)
  # method 1
  ODAL1_meta_median_se_1 = tmp$se_1
  ODAL1_meta_median_CI_l_1 = tmp$lower_1
  ODAL1_meta_median_CI_u_1 = tmp$upper_1
  # method2
  ODAL1_meta_median_se_2 = tmp$se_2
  ODAL1_meta_median_CI_l_2 = tmp$lower_2
  ODAL1_meta_median_CI_u_2 = tmp$upper_2
  
  return(rbind(beta_meta,
                 beta_meta_CI_l,
                 beta_meta_CI_u,
                 ODAL1_meta_mean,
                 ODAL1_meta_mean_CI_l_2,
                 ODAL1_meta_mean_CI_u_2,
                 ODAL1_meta_median,
                 ODAL1_meta_median_CI_l_2,
                 ODAL1_meta_median_CI_u_2,
                 beta_pooled,
                 beta_pooled_CI_l,
                 beta_pooled_CI_u
    ))
}


### forest plot
forest.plot <- function(output_matrix, output_matrix_2,output_matrix_3, title, with_label, margin_left){
  ## xlim depends on the value of CI and point_est
  # plot(0,0,type="n", xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", yaxs="i",
  #      ylim=c(1,7),xlim=c(min(CI[,(2*i-1)])-0.5,max(CI[,(2*i)])+0.5), main="")
  # c(bottom, left, top, right)
  par(mar = c(2, margin_left, 2, 2)) # Set the margin on all sides to 2
  plot(0,0,type="n", xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", yaxs="i",
         ylim = c(1,7), xlim = c(-1,2), main="")
 
  axis(1,cex.axis = 1.5)
  mtext(title, outer=FALSE, line=0.4,cex=1.3)
  
  # segments for beta pooled
  segments(output_matrix["beta_pooled_CI_l"], 1.3,output_matrix["beta_pooled_CI_u"], 1.3, lwd=4)
  segments(output_matrix["beta_pooled"], 1.3-0.1,output_matrix["beta_pooled"], 1.3+0.1, lwd=4)
  segments(output_matrix["beta_pooled_CI_l"], 1.3-0.05,output_matrix["beta_pooled_CI_l"], 1.3+0.05, lwd=4)
  segments(output_matrix["beta_pooled_CI_u"], 1.3-0.05,output_matrix["beta_pooled_CI_u"], 1.3+0.05, lwd=4)
  
  abline(v = output_matrix["beta_pooled"], lty = 2, lwd = 1)
  
  # added by Jessie 09.16
  # segments for beta pooled
  # output_matrix_2 is result for w/o japanese
  segments(output_matrix_2["beta_pooled_CI_l"], 2.3,output_matrix_2["beta_pooled_CI_u"], 2.3, lwd=4, col = "red")
  segments(output_matrix_2["beta_pooled"], 2.3-0.1,output_matrix_2["beta_pooled"], 2.3+0.1, lwd=4, col = "red")
  segments(output_matrix_2["beta_pooled_CI_l"], 2.3-0.05,output_matrix_2["beta_pooled_CI_l"], 2.3+0.05, lwd=4, col = "red")
  segments(output_matrix_2["beta_pooled_CI_u"], 2.3-0.05,output_matrix_2["beta_pooled_CI_u"], 2.3+0.05, lwd=4, col = "red")
  
  abline(v = output_matrix_2["beta_pooled"], lty = 3, lwd = 1, col = "red")
  
  # added by Jessie 09.16
  # segments for beta pooled
  # output_matrix_3 is result for w/o japanese nor mdcr
  segments(output_matrix_3["beta_pooled_CI_l"], 3.3,output_matrix_3["beta_pooled_CI_u"], 3.3, lwd=4, col = "blue")
  segments(output_matrix_3["beta_pooled"], 3.3-0.1,output_matrix_3["beta_pooled"], 3.3+0.1, lwd=4,  col = "blue")
  segments(output_matrix_3["beta_pooled_CI_l"], 3.3-0.05,output_matrix_3["beta_pooled_CI_l"], 3.3+0.05, lwd=4,  col = "blue")
  segments(output_matrix_3["beta_pooled_CI_u"], 3.3-0.05,output_matrix_3["beta_pooled_CI_u"], 3.3+0.05, lwd=4,  col = "blue")
  
  abline(v = output_matrix_3["beta_pooled"], lty = 3, lwd = 1, col = "blue")
  
  # color: 
  color_list = c("#eba833","#8cc56d","#95cbdb")
  # segments for beta meta 
  segments(output_matrix["beta_meta_CI_l"], 4.3,output_matrix["beta_meta_CI_u"], 4.3, lwd=4, col=color_list[3])
  segments(output_matrix["beta_meta"], 4.3-0.1,output_matrix["beta_meta"], 4.3+0.1, lwd=4, col=color_list[3])
  segments(output_matrix["beta_meta_CI_l"], 4.3-0.05,output_matrix["beta_meta_CI_l"], 4.3+0.05, lwd=4, col=color_list[3])
  segments(output_matrix["beta_meta_CI_u"], 4.3-0.05,output_matrix["beta_meta_CI_u"], 4.3+0.05, lwd=4, col=color_list[3])
  
  # segments for beta odal meta median
  segments(output_matrix["ODAL1_meta_median_CI_l_2"], 5.3,output_matrix["ODAL1_meta_median_CI_u_2"], 5.3, lwd=4, col = color_list[2])
  segments(output_matrix["ODAL1_meta_median"], 5.3-0.1,output_matrix["ODAL1_meta_median"], 5.3+0.1, lwd=4, col = color_list[2])
  segments(output_matrix["ODAL1_meta_median_CI_l_2"], 5.3-0.05,output_matrix["ODAL1_meta_median_CI_l_2"], 5.3+0.05, lwd=4, col = color_list[2])
  segments(output_matrix["ODAL1_meta_median_CI_u_2"], 5.3-0.05,output_matrix["ODAL1_meta_median_CI_u_2"], 5.3+0.05, lwd=4, col = color_list[2])
  
  # segments for beta odal meta mean
  segments(output_matrix["ODAL1_meta_mean_CI_l_2"], 6.3,output_matrix["ODAL1_meta_mean_CI_u_2"], 6.3, lwd=4, col = color_list[1])
  segments(output_matrix["ODAL1_meta_mean"], 6.3-0.1,output_matrix["ODAL1_meta_mean"], 6.3+0.1, lwd=4, col = color_list[1])
  segments(output_matrix["ODAL1_meta_mean_CI_l_2"], 6.3-0.05,output_matrix["ODAL1_meta_mean_CI_l_2"], 6.3+0.05, lwd=4, col = color_list[1])
  segments(output_matrix["ODAL1_meta_mean_CI_u_2"], 6.3-0.05,output_matrix["ODAL1_meta_mean_CI_u_2"], 6.3+0.05, lwd=4, col = color_list[1])
  
  if (with_label == TRUE){
    u = par("usr")
    text(u[1]-0.01,1.3,'Pooled\nAll sites',xpd=T,adj=1,cex=1.7)
    text(u[1]-0.01,2.3,'Pooled\nno JMDC',xpd=T,adj=1,cex=1.7)
    text(u[1]-0.01,3.3,'Pooled\nno JMDC\nnor MDCR',xpd=T,adj=1,cex=1.7)
    text(u[1]-0.01,4.3,'Meta\nAnalysis',xpd=T,adj=1,cex=1.7)
    text(u[1]-0.01,5.3,'Robust\nODAL',xpd=T,adj=1,cex=1.7)
    text(u[1]-0.01,6.3,'ODAL',xpd=T,adj=1,cex=1.7)
  }
}

forest.plot_2 <- function(output_matrix, output_matrix_2, output_matrix_3, subtitle){
  # pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/dataAnalysis/AMI_CCAE_update_1.pdf", width = 12, height = 3)
  par(mfrow = c(1,3))
  # forest.plot(output_matrix[,1], title = paste("Intercept",subtitle), 
  #             with_label = TRUE, margin_left = 8)
  forest.plot(output_matrix[,2], output_matrix_2[,2], output_matrix_3[,2],
              title = paste("Obesity",subtitle), 
              with_label = TRUE, margin_left = 6)
  forest.plot(output_matrix[,3], output_matrix_2[,3], output_matrix_3[,3],
              title = paste("Alcohol dependence",subtitle), 
              with_label = FALSE, margin_left = 2)
  forest.plot(output_matrix[,4], output_matrix_2[,4], output_matrix_3[,4],
              title = paste("Hypertensive disorder",subtitle), 
              with_label = FALSE, margin_left = 2)
  # dev.off()
  
  # pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/dataAnalysis/AMI_CCAE_update_2.pdf", width = 12, height = 3)
  par(mfrow = c(1,3))
  forest.plot(output_matrix[,5], output_matrix_2[,5], output_matrix_3[,5],
              title = paste("Major depressive disorder",subtitle), 
              with_label = TRUE, margin_left = 6)
  forest.plot(output_matrix[,6], output_matrix_2[,6], output_matrix_3[,6],
              title = paste("Type 2 diabetes mellitus",subtitle), 
              with_label = FALSE, margin_left = 2)
  forest.plot(output_matrix[,7], output_matrix_2[,7], output_matrix_3[,7],
              title = paste("Hyperlipidemia",subtitle), 
              with_label = FALSE, margin_left = 2)
  # dev.off()
}

output_matrix = clean_output(output1,n=7)
output_matrix_2 = clean_output(output_wo_jp,n=7)
output_matrix_3 = clean_output(output_wo_jp_md,n=7)

# output1
pdf("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/dataAnalysis/AMI_CCAE_update_0916.pdf", width = 12, height = 6)
forest.plot_2(output_matrix, output_matrix_2, output_matrix_3, subtitle = "")
dev.off()


