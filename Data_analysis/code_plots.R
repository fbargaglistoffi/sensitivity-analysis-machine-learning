rm(list = ls())
library(party)

# Classification Tree Analysis

# Below we will check, through a Classification tree (Breiman et al., 1984) to see which are the variable that explain better the low financial literacy score status.


low$dummy <- rep(1, nrow(low))
data$dummy <- rep(0, nrow(data))
data_low <- rbind(low, data)
data_low <- data_low[!duplicated(data_low[,-ncol(data_low)]), ]
table(data_low$dummy)



mf <- model.frame(formula, data=data_low)
mt <- attr(mf, "terms")
predvarnames <- attr(mt, "term.labels")
formula_cart <- as.formula(paste("dummy ~", paste(predvarnames, collapse="+")))
ctree <- ctree(formula_cart, 
               data = data_low, 
               control = ctree_control(mincriterion=0.99, 
                                       minsplit=0, 
                                       minbucket=0,
                                       maxdepth = 3))




plot(ctree, terminal_panel=node_barplot2,
     tp_args = list(ylines = c(2, 4)))


### Classification Tree for Flanders.


ctree_flanders <- ctree(formula_cart, 
                        data = data_low[which(!is.na(data_low$PV1FLIT)),], 
                        control = ctree_control(mincriterion=0.99, 
                                                minsplit=0, 
                                                minbucket=0,
                                                maxdepth = 3))



plot(ctree_flanders, terminal_panel=node_barplot2,
     tp_args = list(ylines = c(2, 4)))


### Classification tree for Wallonia.


ctree_wallonia <- ctree(formula_cart, 
                        data = data_low[which(is.na(data_low$PV1FLIT)),], 
                        control = ctree_control(mincriterion=0.99, 
                                                minsplit=0, 
                                                minbucket=0,
                                                maxdepth = 3))



plot(ctree_wallonia, terminal_panel=node_barplot2,
     tp_args = list(ylines = c(2, 4)))


### Classification tree by grade.

Here, we split our analyisis by grade. We group grade 7-9, and 10-12. Indeed, in Flanders, after grade 9 some things change (teacher qualification requirements etc.). Hence, it is interesting to see which are the drivers of heterogeneous effects by grade.


ctree_10_12_grades <- ctree(formula_cart, 
                            data = data_low[which(data_low$ST001D01T=="Grade 10" |
                                                    data_low$ST001D01T=="Grade 11" |
                                                    data_low$ST001D01T=="Grade 12"),], 
                            control = ctree_control(mincriterion=0.99, 
                                                    minsplit=0, 
                                                    minbucket=0,
                                                    maxdepth = 3))



plot(ctree_10_12_grades, terminal_panel=node_barplot2,
     tp_args = list(ylines = c(2, 4)))



ctree_7_9_grades <- ctree(formula_cart, 
                          data = data_low[which(data_low$ST001D01T=="Grade 7" |
                                                  data_low$ST001D01T=="Grade 8" |
                                                  data_low$ST001D01T=="Grade 9"),], 
                          control = ctree_control(mincriterion=0.99, 
                                                  minsplit=0, 
                                                  minbucket=0,
                                                  maxdepth = 3))



plot(ctree_7_9_grades, terminal_panel=node_barplot2,
     tp_args = list(ylines = c(2, 4)))


### Classification tree by track.

Regarding the track, it is interesting to split by general and vocational education.


ctree_general <- ctree(formula_cart, 
                       data = data_low[which(data_low$ISCEDO=="General"),], 
                       control = ctree_control(mincriterion=0.99, 
                                               minsplit=0, 
                                               minbucket=0,
                                               maxdepth = 3))



plot(ctree_general, terminal_panel=node_barplot2,
     tp_args = list(ylines = c(2, 4)))



ctree_vocational <- ctree(formula_cart, 
                          data = data_low[which(data_low$ISCEDO=="Vocational"),], 
                          control = ctree_control(mincriterion=0.99, 
                                                  minsplit=0, 
                                                  minbucket=0,
                                                  maxdepth = 3))



plot(ctree_vocational, terminal_panel=node_barplot2,
     tp_args = list(ylines = c(2, 4)))