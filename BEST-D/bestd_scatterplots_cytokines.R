########################
# Get variables of interest:
colnames(all_data)
all_data_melt_0 <- melt(all_data, measure.vars = c(29, 88:92))
all_data_melt_0
head(all_data_melt_0)[1:5, 1:5]
dim(all_data_melt_0)
colnames(all_data_melt_0)
plyr::count(all_data_melt_0$variable)
summary(all_data_melt_0$value)
summary(all_data_melt_0[which(all_data_melt_0$variable == 'Ln_TNFalpha0'), 'value'])
all_data_melt_0[1:10, c('value', 'variable')]
all_data_melt_0$variable <- factor(all_data_melt_0$variable, 
                                 levels = c('vitd0',
                                            'Ln_IFNgamma0',
                                            'Ln_IL10_0',
                                            'Ln_IL6_0',
                                            'Ln_IL8_0',
                                            'Ln_TNFalpha0'
                                            ),
                                 labels = c('25OHD baseline',
                                            'IFNg baseline',
                                            'IL10 baseline',
                                            'IL6 baseline',
                                            'IL8 baseline',
                                            'TNFa baseline'
                                            )
                                 )
plyr::count(all_data_melt_0$variable)
plyr::count(all_data_melt_0$arm2)
group <- factor(all_data_melt_0$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group)

# ggplot(data = all_data, aes(x = vitd0, Ln_IFNgamma0, colour = group)) +

ggplot(data = as.data.frame(all_data_melt_0),
         aes(x = value, y = group, colour = group)) +
  facet_wrap(~variable, scale='free_x') +
  geom_point(shape = 1)  +# Hollow circles
  scale_colour_hue(l = 50) + # darker palette
  geom_smooth(method = lm,   # regression line
  #             se = FALSE    # exclude confidence region
  ) + # fullrange = TRUE # Extend regression line
  labs(title = '', x = '', y = '') +
  theme_classic() +
  theme(text = element_text(size = 14), 
        legend.title=element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.position="bottom",
        strip.background = element_blank()
        )
########################

########################
# Get variables of interest:
colnames(all_data)
all_data_melt_12 <- melt(all_data, measure.vars = c(42, 93:97))
all_data_melt_12
head(all_data_melt_12)[1:5, 1:5]
dim(all_data_melt_12)
colnames(all_data_melt_12)
plyr::count(all_data_melt_12$variable)
summary(all_data_melt_12$value)
summary(all_data_melt_12[which(all_data_melt_12$variable == 'Ln_TNFalpha12'), 'value'])
all_data_melt_12[1:10, c('value', 'variable')]
all_data_melt_12$variable <- factor(all_data_melt_12$variable, 
                                   levels = c('vitd12',
                                              'Ln_IFNgamma12',
                                              'Ln_IL10_12',
                                              'Ln_IL6_12',
                                              'Ln_IL8_12',
                                              'Ln_TNFalpha12'
                                   ),
                                   labels = c('25OHD 12m',
                                              'IFNg 12m',
                                              'IL10 12m',
                                              'IL6 12m',
                                              'IL8 12m',
                                              'TNFa 12m'
                                              )
                                   )
plyr::count(all_data_melt_12$variable)
plyr::count(all_data_melt_12$arm2)
group <- factor(all_data_melt_12$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group)

ggplot(data = as.data.frame(all_data_melt_12),
       aes(x = value, y = group, colour = group)) +
  facet_wrap(~variable, scale='free_x') +
  geom_point(shape = 1)  +# Hollow circles
  scale_colour_hue(l = 50) + # darker palette
  geom_smooth(method = lm,   # regression line
              #             se = FALSE    # exclude confidence region
  ) + # fullrange = TRUE # Extend regression line
  labs(title = '', x = '', y = '') +
  theme_classic() +
  theme(text = element_text(size = 14), 
        legend.title=element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.position="bottom",
        strip.background = element_blank()
  )
########################

########################
# Get variables of interest:
colnames(all_data)
all_data_melt_scat <- melt(data = all_data, id.vars = c('vitd0', 'arm2'), 
                           measure.vars = c(88:92))
all_data_melt_scat
head(all_data_melt_scat)
dim(all_data_melt_scat)
colnames(all_data_melt_scat)
plyr::count(all_data_melt_scat$variable)
summary(all_data_melt_scat$value)
summary(all_data_melt_scat[which(all_data_melt_scat$variable == 'Ln_TNFalpha0'), 'value'])
all_data_melt_scat[1:10, c('value', 'variable')]
all_data_melt_scat$variable <- factor(all_data_melt_scat$variable, 
                                     levels = c('vitd0',
                                                'Ln_IFNgamma0',
                                                'Ln_IL10_0',
                                                'Ln_IL6_0',
                                                'Ln_IL8_0',
                                                'Ln_TNFalpha0'
                                     ),
                                     labels = c('25OHD baseline',
                                                'IFNg baseline',
                                                'IL10 baseline',
                                                'IL6 baseline',
                                                'IL8 baseline',
                                                'TNFa baseline'
                                     )
                                     )
plyr::count(all_data_melt_scat$variable)
plyr::count(all_data_melt_scat$arm2)
group <- factor(all_data_melt_scat$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group)

ggplot(data = as.data.frame(all_data_melt_scat),
       aes(x = vitd0, y = variable, colour = group)) +
  facet_wrap(~variable, scale = 'free_x') +
  geom_point(shape = 1)  +# Hollow circles
  scale_colour_hue(l = 50) + # darker palette
  geom_smooth(method = lm   # regression line
              #             se = FALSE    # exclude confidence region
  ) + # fullrange = TRUE # Extend regression line
  labs(title = '', x = '', y = '') +
  theme_classic() +
  theme(text = element_text(size = 14), 
        legend.title=element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.position="bottom",
        strip.background = element_blank()
  )

all_data_melt_scat <- reshape(all_data,  
                              direction = "long",  
                              varying = list(c(29, 88),
                                             c(29, 89),
                                             c(29, 90),
                                             c(29, 91),
                                             c(29, 92)),
                              # sep = "",
                              # v.names = c('Ln_IFNgamma0',
                              #             'Ln_IL10_0',
                              #             'Ln_IL6_0',
                              #             'Ln_IL8_0',
                              #             'Ln_TNFalpha0'
                              #             ),
                              idvar = 'arm2',
                              ids = rownames(all_data)
                              )

head(all_data_melt_scat)










