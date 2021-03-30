## ----setup, include=FALSE-------------------------------------------------------------
library(tibble)
library(ggplot2)
library(tidyr)
library(dplyr)
mytheme <- function(){
  theme(axis.title=element_text(size=15), axis.text=element_text(size=12),
                           text=element_text(size=12), title=element_text(size=15))
}

uwogray <- "#807F83"
purples <- c("#4F2683", "#8463AD", "#26064E", "#643F94", "#3B166A")
pal2 <- c("#4F2683","#807F83")
# pal3 <- c("#4F2683", "#811D7C", "302E87")
pal3 <- c(rgb(79,38,131, maxColorValue = 255), rgb(129,29,124, maxColorValue = 255), rgb(48,46,135, maxColorValue = 255))
# pal4 <- c("#4F2683", "#C1582C", "#1E845A", "#C1B52C")
pal4 <- c(rgb(79,38,131, maxColorValue = 255), rgb(193,88,44, maxColorValue = 255), rgb(30,132,90, maxColorValue = 255), rgb(193,181,44, maxColorValue = 255))




## ----echo=FALSE, fig.align="center", out.width="70%", fig.width=10, fig.height=6------
library(tibble)
library(recipes)
library(dplyr)
library(tidyr)
set.seed(454303)
x <- -rchisq(500, 3)
data <- tibble(skl = -rchisq(1000, 3), 
               skr = rchisq(1000, 3), 
               bm = c(rchisq(500, 3), x - min(x)+.1), 
               norm = rnorm(1000, 0,1)) 
data2 <- recipe(skl + skr + bm + norm ~ 1, data=data)  %>% 
  step_YeoJohnson(all_numeric()) %>% 
  prep(., data=data) %>% 
  bake(., new_data=data)
  
data <- add_column(data, transform=factor(1, 
                                          levels=1:2, 
                                          labels=c("Raw", "Transformed")))

data2 <- add_column(data2, transform=factor(2, 
                                          levels=1:2, 
                                          labels=c("Raw", "Transformed")))

data <- bind_rows(data, data2)
data <- data %>% 
  pivot_longer(cols = -transform) %>% 
  mutate(name = factor(name, 
                            levels=c("skl", "bm", "skr", "norm"), 
                            labels=c("Skew-Left", "Bimodal", "Skew-Right", "Normal")))

library(ggplot2)
ggplot(data, aes(x=value)) + 
  geom_histogram() + 
  facet_wrap(transform ~ name, scales="free", ncol=4) + 
  theme_bw()



## -------------------------------------------------------------------------------------
data(Ericksen, package="carData")
library(recipes)
blueprint <- recipe(undercount ~ ., data=Ericksen) %>% 
  step_YeoJohnson(all_outcomes()) 


## -------------------------------------------------------------------------------------
prepare <- prep(blueprint, training=Ericksen)


## -------------------------------------------------------------------------------------
train_data <- bake(prepare, new_data=Ericksen)


## ----echo=FALSE, fig.align="center", out.width="85%"----------------------------------
knitr::include_graphics("subset_alg.png")


## ---- echo=F, fig.height=7, fig.width=7, out.height="80%", fig.align="center"---------
library(tibble)
n <- 1000
res <- NULL
for(i in 1:2500){
    X <- matrix(rnorm(n*5), ncol=5)
    b <- runif(5, .25,5)
    yhat <- X %*% b
    y <- yhat + rnorm(n, 0, runif(1, .25, 4)*sd(yhat))
    tmp <- lm(y ~ X)
    st <- summary(tmp)
    res <- rbind(res, c(st$r.squared, st$adj.r.squared, AIC(tmp), BIC(tmp), MuMIn::Cp(tmp))) 
    }
colnames(res) <- c("R2", "R2_adj", "AIC", "BIC", "Cp")
res <- as.data.frame(res)
library(GGally)
ggpairs(res) + theme_bw()


## ----modsel, echo=T-------------------------------------------------------------------
## load leaps package, where the regsubsets
## function comes from
library(leaps)
library(car)
data(Ericksen, package="carData")
library(recipes)
data <- recipe(undercount ~ ., data=Ericksen) %>% 
  step_YeoJohnson(all_numeric()) %>%
  step_scale(all_numeric()) %>% 
  step_dummy(all_nominal(), one_hot = FALSE) %>% 
  prep(., data=Ericksen) %>% 
  bake(., new_data=Ericksen)

X <- data %>% 
  select(-undercount) %>% 
  as.matrix()

y <- data %>% 
  select(undercount) %>% 
  pull

rmods <- regsubsets(x=X, y=y, method="exhaustive",
    all.best=TRUE, nbest=10)


## ----subsetfig, echo=T, eval=FALSE----------------------------------------------------
## Load the car package to have access to
## the subsets function
library(car)
## subsets plots the all-subsets results
## from the regsubsets function
## cp is smaller for better models and
## has a lower bound of p (number of parameters)
subsets(rmods, statistic="rsq", legend=F)


## ----echo=FALSE, results='hide', fig.align="center", out.width="85%", fig.height=7, fig.width=15----
s <- subsets(rmods, statistic="rsq", legend=F)


## -------------------------------------------------------------------------------------
library(tibble)
library(dplyr)
mods <- as_tibble(summary(rmods)$which, rownames = "size")
mods <- mods %>% add_column(rsq=summary(rmods)$rsq)
mods <- mods %>% group_by(size) %>% filter(rsq == max(rsq))
n <- names(mods)[3:10]
forms <- apply(mods[,3:10], 1, function(x)
  paste("undercount ~ ", paste(n[which(x)], collapse="+")))
forms <- gsub("city_state", "city", forms)
best_mods <- lapply(forms, function(f)glm(f, data=Ericksen))


## -------------------------------------------------------------------------------------
library(rsample)
library(purrr)
cv_fun <- function(split, ...){
  tmp <- lapply(best_mods, function(x)update(x, .~., data=analysis(split)))
  preds <- sapply(tmp, function(x)predict(x, newdata=assessment(split)))
  mse <- colMeans((c(assessment(split)$undercount)-preds)^2)
  tibble(
    err = mse, 
    size=1:8)
}
set.seed(5390423)
v <- vfold_cv(Ericksen, v=10, repeats = 10) %>%
      mutate(err = map(splits, cv_fun)) %>% 
      unnest(err) %>% 
      group_by(size) %>% 
      summarise(err= mean(err))
v


## -------------------------------------------------------------------------------------
v2 <- vfold_cv(Ericksen, v=10, repeats = 10) %>%
      mutate(err = map(splits, cv_fun)) %>% 
      unnest(err) %>% 
      group_by(size, id) %>% 
      summarise(err= mean(err)) %>% 
      ungroup() %>% 
      group_by(size) %>% 
      summarise(se=sd(err)/sqrt(n()),
                err= mean(err)) %>% 
      mutate(thresh = min(err) +
               se[which(err==min(err))]) %>% 
      filter(err < thresh) %>% 
      slice(1)
v2



## -------------------------------------------------------------------------------------
library(caret)
train.control <- trainControl(method = "cv", number = 10)
back.model <- train(undercount ~., data = data,
                    method = "leapBackward", 
                    nested=TRUE, 
                    tuneGrid = data.frame(nvmax = 1:8),
                    trControl = train.control
                    )
back.model$results
back.model$bestTune


## -------------------------------------------------------------------------------------
for.model <- train(undercount ~., data = data,
                    method = "leapForward", 
                    nested=TRUE, 
                    tuneGrid = data.frame(nvmax = 1:8),
                    trControl = train.control
                    )
for.model$results
for.model$bestTune


## -------------------------------------------------------------------------------------
pick1se <- function(x){
  res <- x$results
  res <- res %>% 
    mutate(across(ends_with("SD"), function(z)z/sqrt(x$control$number), .names="se_{.col}"), 
           RMSE_thresh = min(RMSE) + se_RMSESD[which(RMSE == min(RMSE))], 
           Rsquared_thresh = max(Rsquared) - se_RsquaredSD[which(Rsquared == max(Rsquared))], 
           MAE_thresh = min(MAE) + se_MAESD[which(MAE == min(MAE))])

f1 <- filter(res, RMSE < RMSE_thresh) %>% slice(1)
f1 <- add_column(f1, measure="RMSE", .before="nvmax")
f2 <- filter(res, Rsquared > Rsquared_thresh) %>% slice(1)  
f2 <- add_column(f2, measure="Rsquared", .before="nvmax")
f3 <- filter(res, MAE < MAE_thresh) %>% slice(1)  
f3 <- add_column(f3, measure="MAE", .before="nvmax")
bind_rows(f1,f2,f3) %>% select(measure, nvmax, RMSE, Rsquared, MAE)
}


## -------------------------------------------------------------------------------------
pick1se(back.model)  
pick1se(for.model)  


## -------------------------------------------------------------------------------------
mods %>% filter(size == 5) %>% 
  select(`(Intercept)`:poverty, conventional)
which(summary(back.model)$which[3,])
which(summary(for.model)$which[4,])
which(summary(for.model)$which[3,])


## ---- rcv1, echo=TRUE, eval=FALSE, fig.height=6, fig.width=6, out.width="65%", fig.align="center"----
library(glmnet)
library(ggplot2)
library(rio)
library(tidyr)
banks99 <- import(
  "http://quantoid.net/files/reg3/banks99.dta")
## standardize quantitative variables in banks data
banks99 <- banks99 %>%
  select(-GWNo, -Year, -country) %>%
  mutate(parl_resp = as.factor(parl_resp),
         eff_leg = as.factor(eff_leg))

b99 <- recipe(gdppc_mp ~ ., data=banks99) %>%
  step_knnimpute(all_predictors(), all_outcomes()) %>%
  step_YeoJohnson(all_numeric()) %>%
  step_scale(all_numeric()) %>%
  step_dummy(all_nominal(), one_hot=FALSE) %>%
  prep(., data=banks99) %>%
  bake(., new_data=banks99)

library(glmnet)
loglam <- seq(6.8, -5, length=100)
X <- b99 %>% select(-gdppc_mp) %>% as.matrix()
y <- b99 %>% select(gdppc_mp) %>% pull

g1 <- glmnet(X, y, alpha=0)

rcv <- cv.glmnet(X, y, alpha=0, lambda=exp(loglam))
plot(rcv)



## ----rcv, echo=T, eval=FALSE----------------------------------------------------------
r <- glmnet(X, y, alpha=0, lambda=exp(loglam))
ridge.mod <- glmnet(X,y, alpha=0, lambda=rcv$lambda.min)
mod <- lm(y ~ X)

l2o <- sqrt(sum(coef(mod)^2))
l2r <- apply(r$beta, 2, function(x)sqrt(sum(x^2)))
br <- r$beta %>% as.matrix %>% t %>% as.data.frame
br$ratio <- l2r/l2o
br <- br %>% pivot_longer(under5_mort:all_veh_pc, names_to="variable", values_to="coef")
ggplot(br, aes(x=ratio, y=coef)) +
  geom_line() +
  geom_vline(xintercept=(l2r/l2o)[87], lty=2) +
  facet_wrap(~variable) +
  theme_bw() +
  mytheme() +
  labs(x="Ratio of L2(ridge)/L2(OLS)", y="Coefficient")


## -------------------------------------------------------------------------------------
set.seed(1234)
Sig <- diag(5)
Sig[3:5,3:5] <- .99
diag(Sig) <- 1
X <- MASS::mvrnorm(500,rep(0,5), Sig)
b <- c(1,1,1,0,0)
ystar <- X %*% b
y <- ystar + rnorm(500, 0, 2)


## -------------------------------------------------------------------------------------
summary(m1 <- lm(y~ X))


## ---- echo=FALSE, eval=TRUE, fig.width=7, fig.height=7, out.height="80%", fig.align="center"----
rcv2 <- cv.glmnet(scale(X), y, alpha=0)
plot(rcv2)


## ---- echo=FALSE, eval=TRUE, fig.width=7, fig.height=7, out.height="80%", fig.align="center"----
r2 <- glmnet(scale(X), y, alpha=0, lambda=rcv2$lambda.1se)
coefs <- tibble(
  b = c(as.vector(m1$coef), as.vector(coef(r2))), 
  model = factor(rep(1:2, each=length(coef(m1))), 
                 labels=c("LM", "Ridge")), 
  variable = rep(names(m1$coef), 2))
ggplot(coefs, aes(x=b, y=variable, colour=model)) + 
  geom_point() + 
  theme_bw() + 
  scale_colour_manual(values=pal2) + 
  geom_vline(xintercept=0, lty=3)+ 
  mytheme()


## -------------------------------------------------------------------------------------
library(boot)
df <- cbind(y, X)
boot.ridge <- function(data, inds, ...){
  tmp <- data[inds,]
  y <- tmp[,1]
  X <- tmp[,-1]
  out <- glmnet(X,y,alpha=0, lambda=.8736)
  as.vector(coef(out))
}
br <- boot(statistic=boot.ridge, 
           data=df, R=100)
v <- var(br$t)
pred.vars <- diag(cbind(1, X) %*% 
                    v %*% t(cbind(1,X)))
lm.vars <- diag(cbind(1, X) %*% 
                  vcov(m1) %*% t(cbind(1,X)))


## ---- echo=FALSE, fig.width=7, fig.height=7, out.height="80%", fig.align="center"-----
df <- data.frame(x=lm.vars/pred.vars)
ggplot(df, aes(x=x)) + 
  geom_histogram() + 
  theme_bw() + 
  mytheme() + 
  labs(x="OLS Vars/Ridge Vars")


## -------------------------------------------------------------------------------------
summary(pred.vars/lm.vars)
ridge.preds <- cbind(1, X) %*% coef(r2) 
lm.preds <- cbind(1,X) %*% coef(m1)
cor(as.vector(lm.preds), as.vector(ridge.preds))


## ----lasso1, echo=T-------------------------------------------------------------------
X <- b99 %>% select(-gdppc_mp) %>% as.matrix()
y <- b99 %>% select(gdppc_mp) %>% pull
m1 <- lm(y ~ X)
cvg <- cv.glmnet(X,y, lambda=exp(loglam))
## estimate the LASSO for the Banks data
g <- glmnet(X, y, lambda=cvg$lambda.min)
## compare coefficients with the linear model


## -------------------------------------------------------------------------------------
round(cbind(coef(cvg), coef(mod)), 4)


## ----lcv, echo=T, eval=FALSE----------------------------------------------------------
r <- glmnet(X, y, alpha=0, lambda=exp(loglam))
g <- glmnet(X, y, alpha=1, lambda=exp(loglam))
mod <- lm(y ~ X)


br1 <- r$beta %>% as.matrix %>% t %>% as.data.frame
br2 <- g$beta %>% as.matrix %>% t %>% as.data.frame
br1$lambda <- br2$lambda <- loglam
br1 <- br1 %>% pivot_longer(under5_mort:all_veh_pc, names_to="variable", values_to="coef")
br2 <- br2 %>% pivot_longer(under5_mort:all_veh_pc, names_to="variable", values_to="coef")
br1$model <- factor(1, levels=c(1,2), labels=c("Ridge", "LASSO"))
br2$model <- factor(2, levels=c(1,2), labels=c("Ridge", "LASSO"))

br <- bind_rows(br1, br2)

ggplot(br, aes(x=lambda, y=coef, colour=model)) +
  geom_line() +
  facet_wrap(~variable, scales="free_y") +
  geom_vline(xintercept=log(rcv$lambda.min), col=pal2[1], lty=3) +
  geom_vline(xintercept=log(cvg$lambda.min), col=pal2[2], lty=3) +
  scale_colour_manual(values=pal2) +
  theme_bw() +
  mytheme() +
  labs(x="log(Lambda)", y="Coefficient")


## -------------------------------------------------------------------------------------
r1 <- glmnet(X, y, alpha=0, lambda=rcv$lambda.min)
g1 <- glmnet(X, y, alpha=1, lambda=cvg$lambda.min)
yhat <- mod$fitted.values
names(yhat) <- NULL

preds <- tibble(
  ols=yhat, 
  ridge = as.vector(predict(r1, newx=X)), 
  lasso= as.vector(predict(g1, newx=X))
)

## ---- echo=FALSE, fig.height=8 ,fig.width=8, out.height="80%", fig.align="center"-----
ggpairs(preds) + mytheme() + theme_bw()


## ---- echo=FALSE----------------------------------------------------------------------
set.seed(1234)
Sig <- diag(5)
Sig[3:5,3:5] <- .99
diag(Sig) <- 1
coll <- list()
coll$X <- MASS::mvrnorm(500,rep(0,5), Sig)
b <- c(1,1,1,0,0)
coll$ystar <- coll$X %*% b
coll$y <- coll$ystar + rnorm(500, 0, 2)


## -------------------------------------------------------------------------------------
cvg2 <- cv.glmnet(scale(coll$X), coll$y, alpha=1)
rcv2 <- cv.glmnet(scale(coll$X), coll$y, alpha=0)
m1 <- lm(coll$y ~ scale(coll$X))
g2 <- glmnet(scale(coll$X), coll$y, 
             alpha=1, lambda=cvg2$lambda.min)
r2 <- glmnet(scale(coll$X), coll$y, 
             alpha=0, lambda=rcv2$lambda.1se)
coefs <- tibble(
  b = c(as.vector(m1$coef), as.vector(coef(r2)), 
        as.vector(coef(g2))), 
  model = factor(rep(1:3, each=length(coef(m1))), 
          labels=c("LM", "Ridge", "LASSO")), 
  variable = rep(names(m1$coef), 3))
p1 <- ggplot(coefs, aes(x=b, y=variable, 
                  colour=model, shape=model)) + 
  geom_point() + 
  theme_bw() + 
  scale_colour_manual(values=pal3) + 
  geom_vline(xintercept=0, lty=3)+ 
  mytheme()



## ----echo=FALSE, fig.height=6, fig.width=7, out.height="80%", fig.align="center"------
p1


## ----echo=T---------------------------------------------------------------------------
X <- b99 %>% select(-gdppc_mp) %>% as.matrix()
y <- b99 %>% select(gdppc_mp) %>% pull
cv.enet <- list()
s <- seq(0.01, .99, length=25)
for(i in 1:length(s)){
	cv.enet[[i]] <- cv.glmnet(X, y, alpha = s[i])
}
cv.err <- sapply(cv.enet, function(x)min(x$cvm))
s[which.min(cv.err)]


## ----enet, echo=T---------------------------------------------------------------------
b <- sapply(cv.enet[c(1,7,25)], 
	function(x)as.matrix(coef(x)))

plot.dat <- data.frame(
    b = c(b), 
    group = as.factor(rep(c(.010,.255,.990), 
        each = 24)),
    var = factor(rep(
        rownames(coef(cv.enet[[1]])), 3)))
plot.dat <- filter(plot.dat, var != "(Intercept)")
library(ggplot2)
g <- ggplot(plot.dat, aes(x=b, 
    y=reorder(var, b, mean))) + 
    geom_point() + 
    scale_colour_manual(values=pal3) +
    theme_bw() + 
    facet_wrap(~group, nrow=1) + 
    mytheme() + 
    ylab("")
g


## ---- echo=FALSE, fig.height=6, fig.width=7, out.height="80%", fig.align="center"-----
cv2.enet <- list()
for(i in 1:length(s)){
	cv2.enet[[i]] <- cv.glmnet(scale(coll$X), coll$y, alpha = s[i])
}
cv2.err <- sapply(cv2.enet, function(x)min(x$cvm))
alph <- s[which.min(cv2.err)]
cvg3 <- cv.glmnet(scale(coll$X), coll$y, alpha=alph)
g3 <- glmnet(scale(coll$X), coll$y, 
             alpha=alph, lambda=cvg3$lambda.min)

coefs <- tibble(
  b = c(as.vector(m1$coef), as.vector(coef(r2)), 
        as.vector(coef(g2)), as.vector(coef(g3))), 
  model = factor(rep(1:4, each=length(coef(m1))), 
          labels=c("LM", "Ridge", 
                   "LASSO", "E-Net")), 
  variable = rep(names(m1$coef), 4))
ggplot(coefs, aes(x=b, y=variable, 
                  colour=model, shape=model)) + 
  geom_point() + 
  theme_bw() + 
  scale_colour_manual(values=pal4) + 
  geom_vline(xintercept=0, lty=3)+ 
  mytheme()


## -------------------------------------------------------------------------------------
banks_tc <- trainControl(method = "cv", number = 10)
enet_res <- train(gdppc_mp ~., data = b99,
                    method = "glmnet", 
                    nested=TRUE, 
                    trControl = banks_tc, 
                    tuneLength=10
                    )

enet_res$results %>%
  filter(alpha == enet_res$bestTune$alpha &  
         lambda == enet_res$bestTune$lambda)


## ----gge, eval=FALSE------------------------------------------------------------------
ggplot(enet_res) + theme_bw()

## ----alasso, echo=T-------------------------------------------------------------------
set.seed(592043)
# estimate initial ridge regression and save coefficients
b.ridge <- coef(cv.glmnet(X,y, alpha=0))
b.lm <- coef(lm(y ~ X))
b.lasso <- coef(cv.glmnet(X,y, alpha=1))
b.lasso <- ifelse(b.lasso == 0, .01, b.lasso)
b.enet <- coef(glmnet(X, y, alpha=.3, lambda=0.066))
b.enet <- ifelse(b.enet == 0, .01, b.enet)
# calculate weights
gamma <- 1
w.ridge <- 1/(abs(b.ridge)^gamma)
w.lm <- 1/(abs(b.lm)^gamma)
w.lasso <- 1/(abs(b.lasso)^gamma)
w.enet <- 1/(abs(b.enet)^gamma)
# estimate the LASSO with the weights
cval.ridge <- cv.glmnet(X,y, penalty.factor=w.ridge)
cval.lasso <- cv.glmnet(X,y, penalty.factor=w.lasso)
cval.enet <- cv.glmnet(X,y, penalty.factor=w.enet)
cval.lm <- cv.glmnet(X,y, penalty.factor=w.lm)


## ----echo=FALSE-----------------------------------------------------------------------
coefs <- cbind(coef(cval.lm), coef(cval.ridge), 
      coef(cval.enet), coef(cval.lasso))
colnames(coefs) <- c("lm", "ridge", "enet", "lasso")
round(coefs,3)


## -------------------------------------------------------------------------------------
enet_tc <- trainControl(method = "cv", number = 10)
coll.dat <- cbind(coll$y, coll$X)
colnames(coll.dat) <- c("y", paste0("X", 1:5))
coll.dat <- as.data.frame(coll.dat)
enet_res <- train(y ~ ., data = coll.dat,
                    method = "glmnet", 
                    nested=TRUE, 
                    trControl = enet_tc, 
                    tuneLength=10
                    )

b.enet <- coef(glmnet(X, y, 
    alpha=enet_res$bestTune$alpha,
    lambda=enet_res$bestTune$lambda))
b.enet <- ifelse(b.enet == 0, .01, b.enet)
# calculate weights
gamma <- 1
w <- 1/(abs(b.enet)^gamma)
# estimate the LASSO with the weights
cval <- cv.glmnet(scale(coll$X),coll$y, 
                  penalty.factor=w)



## ---- echo=FALSE, fig.height=6, fig.width=7, out.height="80%", fig.align="center"-----
coefs <- tibble(
  b = c(as.vector(m1$coef), as.vector(coef(r2)), 
        as.vector(coef(g2)), as.vector(coef(cval))), 
  model = factor(rep(1:4, each=length(coef(m1))), 
          labels=c("LM", "Ridge", 
                   "LASSO", "A-LASSO")), 
  variable = rep(names(m1$coef), 4))
ggplot(coefs, aes(x=b, y=variable, 
                  colour=model, shape=model)) + 
  geom_point() + 
  theme_bw() + 
  scale_colour_manual(values=pal4) + 
  geom_vline(xintercept=0, lty=3)+ 
  mytheme()


## -------------------------------------------------------------------------------------
r1 <- glmnet(X, y, alpha=0, lambda=rcv$lambda.min)
g1 <- glmnet(X, y, alpha=1, lambda=cvg$lambda.min)
g3 <- glmnet(X, y, alpha=s[which.min(cv.err)], 
             lambda=cvg$lambda.min)
g4 <- glmnet(X, y, alpha=1, penalty.factor=w, 
             lambda=cval$lambda.min)
yhat <- mod$fitted.values
names(yhat) <- NULL

preds <- tibble(
  ols=yhat, 
  ridge = as.vector(predict(r1, newx=X)), 
  lasso= as.vector(predict(g1, newx=X)), 
  enet = as.vector(predict(g3, newx=X)),
  alasso = as.vector(predict(g4, newx=X)))

## ---- echo=FALSE, fig.height=8 ,fig.width=8, out.height="80%", fig.align="center"-----
ggpairs(preds) + mytheme() + theme_bw()

