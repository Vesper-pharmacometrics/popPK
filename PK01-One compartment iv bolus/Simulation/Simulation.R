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

# 使用rxode2进行仿真
# 定义一室模型，2男，2女
model <- rxode2({
  # 使用群体典型值定义个体的参数
  CL = TVCL * exp(eta.CL)
  V = TVV * exp(eta.V)

  # 定义ODE模型, 注意此处是药量的微分方程
  d/dt(central) = -(CL/V) * central 

  # 定义输出变量：浓度
  Cpred = central/V
  Cobs = Cpred * exp(eps.prop)
})

# 定义参数
theta = c(TVCL = 0.2, TVV = 10) # CL的单位是L/h，V的单位是L
omega <- lotri({ eta.CL + eta.V ~ c(0.1, 0.03, 0.1) }) # 顺序var(CL)， cor, var(V)
sigma = lotri(eps.prop ~ 0.04)

# 定义初始值
ini_iv = c(central = 0)

# 定义给药方案
ev_bolus <- eventTable()
ev_bolus$add.dosing(dose = 10, start.time = 0, cmt='central') 
ev_bolus$add.sampling(seq(from=0, to=160, by=20))

# 求解ODE方程
set.seed(123)
output = rxSolve(object=model, 
                 params=theta, omega=omega, sigma=sigma, inits=ini_iv, 
                 events = ev_bolus, nSub = 10)
output

# 查看结果
out_df = as_tibble(output) %>% 
  rename(ID = sim.id, TIME= time, DV = Cobs)


out_plot = ggplot(out_df, aes(x=TIME, y=Cpred)) +
    stat_summary(fun=mean, geom='line') +
    geom_point() +
    scale_x_continuous(breaks = seq(0, 160, by = 20)) +
    theme_bw()
out_plot
# 保存原始数据
write.csv(out_df, "results/simulation_data_raw.csv",na = ".")



# 转换成NONMEM格式
## 1) 观测记录 ---------------------------------------------------------------
obs_record <- out_df %>%
  filter(TIME != 0) %>%
  mutate(
    AMT  = NA_real_,
    EVID = 0, # 表示观测值
    MDV  = 0, # 不缺失观测值
    CMT  = 1, # 观测在中央室
    DV   = round(DV, 2)
  ) %>%
  filter(!is.na(DV)) %>%
  select(ID, TIME, DV, MDV, EVID, CMT, AMT)

## 2) 给药事件 ---------------------------------------------------------------
admin_event <- data.frame(
  ID   = 1:length(unique(out_df$ID)),
  TIME = 0,
  AMT  = 10,
  DV   = NA_real_,
  MDV  = 1,
  EVID = 1,
  CMT  = 1
)


## 3) 合并、排序、加 C 列、NA→"."  -------------------------------
dataset_nonmem <- bind_rows(obs_record, admin_event) %>%
  arrange(ID, TIME, desc(EVID), desc(MDV)) %>%  # 同时刻先给药(EVID=1)后观测(EVID=0)
  mutate(C = NA_real_) %>%                 # 数值型占位符，避免字符型干扰分析
  select(C, everything())        # C 列放最前

str(dataset_nonmem)
write.csv(dataset_nonmem, "results/simulation_data.csv",
            na = ".",row.names = FALSE)

