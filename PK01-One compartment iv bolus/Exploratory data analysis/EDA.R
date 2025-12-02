# Author: Vesper
# Date: 2025/12/02
# Purpose: Build the dataset for PK1 One-compartment iv bolus

# 环境准备
rm(list=ls()) # 清除所有对象

# 导入包
library(ggplot2)
library(dplyr)
library(tidyr)
library(rxode2)
library(readr)
library(ncappc)

# 读取数据
simulation_nonmem_data <- read.csv("results/simulation_data.csv",na.strings = ".")

# 检查数据类型
str(simulation_nonmem_data)

# 绘图查看数据的基本特征
plot = ggplot(simulation_data, aes(x = TIME, y = DV, color = ID)) +
  stat_summary(fun=mean, geom='line') +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 160, by = 20)) +
  theme_bw()
plot

# 初始参数估计--NCA

nca_result = ncappc(
    obsFile   = simulation_nonmem_data,
    adminType = 'iv-bolus',
    # LambdaTimeRange = c(120,160), # 手动设置进行末端回归的时间范围
    extrapolate = TRUE, # 此处不设置为TRUE则不进行外推，CL和V就计算不出来
    method    = 'linearup-logdown',
    onlyNCA   = TRUE,
    noPlot    = TRUE,
    printOut  = FALSE)

nca_result$ncaOutput$Rsq
nca_result$ncaOutput$Rsq_adjusted
nca_result$ncaOutput$No_points_Lambda_z
# 查看关键结果（应包含 CL/V）
print(dplyr::select(nca_result$ncaOutput, ID, AUClast, AUCINF_obs, HL_Lambda_z, Cl_obs, Vz_obs, No_points_Lambda_z))

mean_nca_cl = mean(nca_result$ncaOutput$Cl_obs)
mean_nca_vz = mean(nca_result$ncaOutput$Vz_obs)


