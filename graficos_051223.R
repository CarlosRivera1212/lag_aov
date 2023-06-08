library(tidyverse)

setwd('~/CARLOS/UNAL/D/articulos/rezago_espacial')

path_graf = 'graf/060723/'
# alfas_df = readRDS('alfas_df.rds')
# alfas_df = readRDS('alfas_df_trt3blq2rep10sim1000.rds')
data = readRDS('sim_trt3blq2rep10sim1000_052023.rds')
alfas_df = readRDS('alfas_df_trt3blq2rep10sim1000_052023.rds')
colnames(alfas_df)
var_graf = 'OT'

alfas_df = alfas_df %>% 
  mutate(name = str_replace(name, 'ord', 'ord_'))

tibble(alfas_df)

x = alfas_df %>% 
  filter(name == 'ord_12') %>% 
  pull(OT)
quantile(x, 0.95)
qchisq(0.95, 1, 6.2)

pc = 6.2
xq = qchisq(seq(0,1, l=1000), df = 1, ncp = pc)
yq = dchisq(xq, df = 1, ncp = pc)

hist(x, probability = T, breaks = 100)
abline(v = quantile(x, 0.95))
lines(xq, yq, col='red')

int_sel = c(0, 12, 24, 36, 48, 60)

alfas_df_graf = alfas_df %>% 
  filter(name %in% paste0('ord_', int_sel))
unique(alfas_df$name)


names(data) = paste0('ord_', seq(0,60,2))
data_df = plyr::ldply(data) %>% 
  mutate(txt = paste0('T[', str_sub(trt, 2), ']~',
                      'B[', str_sub(blq, 2), ']'))
tibble(data_df)

# FIGURA 1 - SPATIAL GRID
fig1 = data_df %>% 
  # filter(.id %in% paste0('ord_', c(0, 36, 60))) %>% 
  filter(.id %in% paste0('ord_', c(0, 60))) %>%
  ggplot()+
  aes(x,y,fill=sim1, label=txt)+
  geom_tile()+
  geom_text(parse = TRUE, color='white')+
  # coord_equal()+
  facet_wrap(~.id)+
  # scale_fill_viridis_b()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

fig1
ggsave('fig1.tif', fig1, 'tiff', path_graf,
       0.9, 16, 9, 'cm')



# FIGURA 2 - MORAN INDEX
fig2 = alfas_df_graf %>% 
  ggplot()+
  aes(pmoran, pmoran_aov)+
  geom_point(size=0.1)+
  facet_wrap(~name)+
  labs(x = 'p-value Moran index in Response',
       y = 'p-value Moran index in Residuals')+
  geom_vline(xintercept = 0.05, color='red')+
  geom_hline(yintercept = 0.05, color='red')


fig2 = alfas_df_graf %>% 
  mutate(pvmoran_y = ifelse(pmoran<0.05, 'Non-Independent', 'Independent'),
         pvmoran_r = ifelse(pmoran_aov<0.05, 'Non-Independent', 'Independent')) %>% 
  mutate(pvmoran_y = factor(pvmoran_y, c('Non-Independent', 'Independent'), ordered = T),
         pvmoran_r = factor(pvmoran_r, c('Non-Independent', 'Independent'), ordered = T)) %>% 
  group_by(name, pvmoran_y, pvmoran_r) %>% 
  summarise(n=n(),
            p = 100*n/1000,
            p_txt = paste0(p, ' %')) %>% 
  ggplot()+
  aes(pvmoran_y, pvmoran_r, label=p_txt, fill=p)+
  geom_tile()+
  geom_text(color='white')+
  geom_vline(xintercept = 1.5)+
  geom_hline(yintercept = 1.5)+
  facet_wrap(~name, ncol=2)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  labs(x = 'p-value Moran index in Response',
       y = 'p-value Moran index in Residuals')+
  theme(legend.position = 'none',
        panel.grid = element_blank())

fig2
ggsave('fig2.tif', fig2, 'tiff', path_graf,
       0.9, 16, 9, 'cm')


# ALFA

fig3 = alfas_df_graf %>% 
  ggplot()+
  aes(alfa, fill=name)+
  geom_density(alpha=0.3)+
  labs(x = 'Overlap Spatial Coefficient Estimated',
       y = 'Density')+
  scale_fill_viridis_d()

fig3
ggsave('fig3.tif', fig3, 'tiff', path_graf,
       0.9, 16, 9, 'cm')

# OT
fig4 = alfas_df_graf %>% 
  ggplot()+
  aes(OT, fill=name)+
  geom_density()+
  geom_text(data = alfas_df_graf %>%
               group_by(name) %>%
               summarise(xOT = max(OT)*0.9,
                         yOT = max(density(OT)$y)*0.9,
                         sOT = sum(OT > qchisq(0.95, 1))),
             aes(x=xOT, y=yOT, label = sOT))+
  facet_wrap(~name, scales = 'free')+
  geom_vline(xintercept = qchisq(0.95, 1), col='red')+
  labs(x = 'Spatial Overlap Statistic (OT)',
       y = 'Density')+
  scale_fill_viridis_d()+
  theme(legend.position = 'none')

fig4
ggsave('fig4.tif', fig4, 'tiff', path_graf,
       0.9, 16, 9, 'cm')

# ALFA - OT
library(ggExtra)

prefig5 = alfas_df_graf %>% 
  ggplot()+
  aes(alfa, OT, color=name, fill=name)+
  geom_point()+
  geom_hline(yintercept = qchisq(0.95, 1))+
  # scale_x_log10()+
  # scale_y_log10()+
  labs(x='Overlap Spatial Coefficient Estimated',
       y='Spatial Overlap Statistic (OT)')+
  scale_color_viridis_d()+
  theme(legend.position = 'bottom')+
  guides(color=guide_legend(nrow=1))

prefig5

ggMarginal(prefig5, groupFill = TRUE, type='violin')
fig5 = ggMarginal(prefig5, groupFill = TRUE,
                  type='boxplot', size = 2)

fig5
ggsave('fig5.tif', fig5, 'tiff', path_graf,
       1.3, 16, 9, 'cm')


# p5a = alfas_df_graf %>% 
#   ggplot()+
#   aes(alfa, name, fill=name)+
#   geom_boxplot()+
#   theme_void()+
#   theme(legend.position = 'none')
# p5b = alfas_df_graf %>% 
#   ggplot()+
#   aes(name, OT, fill=name)+
#   geom_boxplot()+
#   theme_void()+
#   theme(legend.position = 'none')
# 
# gridExtra::grid.arrange(p5a, p5, p5b)
#   

# # # # # # # # # # # # # 
# saveRDS(Z, 'Z_matrix.rds')
# saveRDS(X, 'X_matrix.rds')
# writexl::write_xlsx(data_df, 'datos_trt3blq2_1000sim_30iter_052023.xlsx')
