plot_catdes <- function (x, show = "all", level = 0.01, sort = NULL, col.upper = "indianred2", 
                         col.lower = "royalblue1", numchar = 10, barplot = FALSE, 
                         ...) 
{
  if (!attr(x, "class")[1] == "catdes") {
    stop("'x' must be of type catdes")
  }
  if (show != "all" && show != "quanti" && show != "quali") {
    stop("Invalid value of show")
  }
  if (!barplot) {
    if (!is.null(x$quanti)) {
      rows <- names(x$quanti)
    }
    else {
      rows <- names(x$category)
    }
    if (is.null(rows)) {
      stop("Invalid value of x")
    }
    nb_cluster = length(rows)
    if (!is.null(sort)) {
      if (is.character(sort)) {
        sort = which(rows == sort)
        if (length(sort) == 0) {
          stop("Invalid value of sort")
        }
      }
      else if (!is.double(sort) || sort != round(sort) || 
               sort > nb_cluster || sort < 1) {
        stop("Invalid value of sort")
      }
    }
    if (!col.upper %in% colours() || !col.lower %in% colours()) {
      stop("Invalid colour name. Please use colours() to see the colours available.")
    }
    if (level < 0) {
      stop("The level must be positive")
    }
    col.upper.rgb <- t(col2rgb(col.upper)/255)
    col.lower.rgb <- t(col2rgb(col.lower)/255)
    pval = c()
    if (is.null(sort)) {
      if (show == "quanti" || show == "all") {
        quanti = c()
        for (q in x$quanti) {
          quanti = union(quanti, rownames(q))
          pval = c(pval, q[, "p.value"])
        }
      }
      if (show == "quali" || show == "all") {
        quali = c()
        for (q in x$category) {
          quali = union(quali, rownames(q))
          pval = c(pval, q[, "p.value"])
        }
      }
      if (show == "all") {
        columns = c(quanti, quali)
      }
      else if (show == "quali") {
        columns = quali
      }
      else if (show == "quanti") {
        columns = quanti
      }
    }
    else {
      quanti_all = c()
      pval = c()
      for (q in x$quanti) {
        quanti_all = union(quanti_all, rownames(q))
        pval = c(pval, q[, "p.value"])
      }
      quali_all = c()
      for (q in x$category) {
        quali_all = union(quali_all, rownames(q))
        pval = c(pval, q[, "p.value"])
      }
      k <- sort
      if ((show == "quanti" || show == "all") && !(is.null(x$quanti))) {
        pvals_k_quanti <- as.data.frame(x$quanti[[k]][, 
                                                      "p.value"])
        colnames(pvals_k_quanti) <- "pvals"
        names <- rownames(pvals_k_quanti)
        rownames(pvals_k_quanti) <- NULL
        pvals_k_quanti <- cbind(names, pvals_k_quanti)
        quanti_all <- as.data.frame(cbind(quanti_all, 
                                          rep(0, length(quanti_all))))
        colnames(quanti_all) <- c("names", "pvals")
        quanti_all <- merge(x = quanti_all, y = pvals_k_quanti, 
                            by = "names", all.x = TRUE)[, c("names", "pvals.y")]
        colnames(quanti_all) <- c("names", "pvals")
      }
      if ((show == "quali" || show == "all") && !(is.null(x$category))) {
        pvals_k_quali <- as.data.frame(x$category[[k]][, 
                                                       "p.value"])
        colnames(pvals_k_quali) <- "pvals"
        names <- rownames(pvals_k_quali)
        rownames(pvals_k_quali) <- NULL
        pvals_k_quali <- cbind(names, pvals_k_quali)
        quali_all <- as.data.frame(cbind(quali_all, rep(0, 
                                                        length(quali_all))))
        colnames(quali_all) <- c("names", "pvals")
        quali_all <- merge(x = quali_all, y = pvals_k_quali, 
                           by = "names", all.x = TRUE)[, c("names", "pvals.y")]
        colnames(quali_all) <- c("names", "pvals")
      }
      if (show == "all") {
        columns = rbind(quali_all, quanti_all)
      }
      if (show == "quali") {
        columns = rbind(quali_all)
      }
      if (show == "quanti") {
        columns = rbind(quanti_all)
      }
      columns <- columns[order(columns$pvals), ]
      columns <- columns$names
    }
    dim1 <- length(rows)
    dim2 <- length(columns)
    plot.new()
    for (i in 1:dim1) {
      ybottom = 1 - i * (1/(dim1 + 1))
      xright = 1/(dim2 + 1)
      ytop = 1 - (i + 1) * (1/(dim1 + 1))
      xleft = 0
      center <- c(mean(c(xleft, xright)), mean(c(ytop, 
                                                 ybottom)))
      text(center[1], center[2], rows[i], cex = 0.5, font = 2)
      for (j in 1:(dim2)) {
        xleft = j * (1/(dim2 + 1))
        ybottom = 1 - i * (1/(dim1 + 1))
        xright = (j + 1) * (1/(dim2 + 1))
        ytop = 1 - (i + 1) * (1/(dim1 + 1))
        rect(xleft, ybottom = ybottom, xright, ytop)
      }
    }
    for (j in 1:(dim2)) {
      xleft = j * (1/(dim2 + 1))
      ybottom = 1
      xright = (j + 1) * (1/(dim2 + 1))
      ytop = 1 - (1/(dim1 + 1))
      center <- c(mean(c(xleft, xright)), mean(c(ytop, 
                                                 ybottom)))
      text(center[1], 1 - (1/(dim1 + 1)), columns[j], cex = 0.6, 
           srt = 45, adj = c(0, 0))
    }
    max <- level
    min <- min(pval)
    range <- max - min
    for (i in 1:dim1) {
      quanti = x$quanti
      quali = x$category
      if (show == "quanti" || show == "all") {
        for (name_row in rownames(quanti[[i]])) {
          j <- which(columns == name_row)
          p <- quanti[[i]][name_row, "p.value"]
          if (p < level) {
            col <- pmin(1, pmax(0, 0.9 * (p - min)/(min - 
                                                      max) + 1))
            xleft = j * (1/(dim2 + 1))
            ybottom = 1 - i * (1/(dim1 + 1))
            xright = (j + 1) * (1/(dim2 + 1))
            ytop = 1 - (i + 1) * (1/(dim1 + 1))
            if (quanti[[i]][name_row, "v.test"] >= 0) {
              rect(xleft, ybottom = ybottom, xright, 
                   ytop, col = rgb(col.upper.rgb, alpha = col), 
                   border = NULL)
              tmp.center <- c(mean(c(xleft, xright)), mean(c(ybottom, ytop)))
              text(tmp.center[1], tmp.center[2], labels = round(quanti[[i]][name_row,"Mean in category"],2),cex = 0.6, col="black")
            }
            else {
              rect(xleft, ybottom = ybottom, xright, 
                   ytop, col = rgb(col.lower.rgb, alpha = col), 
                   border = NULL)
              tmp.center <- c(mean(c(xleft, xright)), mean(c(ybottom, ytop)))
              text(tmp.center[1], tmp.center[2], labels = round(quanti[[i]][name_row,"Mean in category"],2),cex = 0.6, col="black")
            }
          }
          else {
            col <- pmin(1, pmax(0, 0.9 * (p - min)/(min - 
                                                      max) + 1))
            xleft = j * (1/(dim2 + 1))
            ybottom = 1 - i * (1/(dim1 + 1))
            xright = (j + 1) * (1/(dim2 + 1))
            ytop = 1 - (i + 1) * (1/(dim1 + 1))
            tmp.center <- c(mean(c(xleft, xright)), mean(c(ybottom, ytop)))
            text(tmp.center[1], tmp.center[2], labels = round(quanti[[i]][name_row,"Mean in category"],2),cex = 0.6, col="black")
          }
        }
      }
      if (show == "quali" || show == "all") {
        for (name_row in rownames(quali[[i]])) {
          j <- which(columns == name_row)
          p <- quali[[i]][name_row, "p.value"]
          if (p < level) {
            col <- pmin(1, pmax(0, 0.9 * (p - min)/(min - 
                                                      max) + 1))
            xleft = j * (1/(dim2 + 1))
            ybottom = 1 - i * (1/(dim1 + 1))
            xright = (j + 1) * (1/(dim2 + 1))
            ytop = 1 - (i + 1) * (1/(dim1 + 1))
            if (quali[[i]][name_row, "v.test"] >= 0) {
              rect(xleft, ybottom = ybottom, xright, 
                   ytop, col = rgb(col.upper.rgb, alpha = col), 
                   border = NULL)
              tmp.center <- c(mean(c(xleft, xright)), mean(c(ybottom, ytop)))
              text(tmp.center[1], tmp.center[2], labels = paste0(round(quali[[i]][name_row,"Mod/Cla"],2)," %"),cex = 0.6, col="black")
            }
            else {
              rect(xleft, ybottom = ybottom, xright, 
                   ytop, col = rgb(col.lower.rgb, alpha = col), 
                   border = NULL)
              tmp.center <- c(mean(c(xleft, xright)), mean(c(ybottom, ytop)))
              text(tmp.center[1], tmp.center[2], labels = paste0(round(quali[[i]][name_row,"Mod/Cla"],2)," %"),cex = 0.6, col="black")
            }
          }
          else  {
            col <- pmin(1, pmax(0, 0.9 * (p - min)/(min - 
                                                      max) + 1))
            xleft = j * (1/(dim2 + 1))
            ybottom = 1 - i * (1/(dim1 + 1))
            xright = (j + 1) * (1/(dim2 + 1))
            ytop = 1 - (i + 1) * (1/(dim1 + 1))
            tmp.center <- c(mean(c(xleft, xright)), mean(c(ybottom, ytop)))
            text(tmp.center[1], tmp.center[2], labels = paste0(round(quali[[i]][name_row,"Mod/Cla"],2)," %"),cex = 0.6, col="black")
          }
        }
      }
    }
  }
  else {
    lengthX = max(length(x$quanti), length(x$category))
    long = rep(0, lengthX)
    list.catdes = list(long)
    minimum = 0
    maximum = 0
    count = 0
    for (i in 1:lengthX) {
      if (!is.null(x$quanti[[i]])) {
        quanti = as.data.frame(x$quanti[[i]])
        quanti.catdes = as.vector(quanti[, 1])
        names(quanti.catdes) = rownames(quanti)
      }
      else quanti.catdes = NULL
      if (!is.null(x$category[[i]])) {
        category = as.data.frame(x$category[[i]])
        category.catdes = as.vector(category[, 5])
        names(category.catdes) = rownames(category)
      }
      else category.catdes = NULL
      if (show == "all") 
        catdes.aux = c(quanti.catdes, category.catdes)
      if (show == "quanti") 
        catdes.aux = quanti.catdes
      if (show == "quali") 
        catdes.aux = category.catdes
      if (!is.null(catdes.aux)) {
        count = count + 1
        long[i] = length(catdes.aux)
        minimum = min(catdes.aux, minimum)
        maximum = max(catdes.aux, maximum)
      }
      else long[i] = 0
      list.catdes[[i]] = catdes.aux
    }
    if (count != 0) {
      if (count <= 4) {
        numc = count
        numr = 1
      }
      else {
        numc = 4
        numr = count%/%4 + 1
      }
      plot.new()
      par(las = 3)
      par(mfrow = c(numr, numc))
      for (i in 1:lengthX) {
        catdes.aux <- list.catdes[[i]]
        if (!is.null(catdes.aux)) {
          catdes.aux = sort(catdes.aux, decreasing = FALSE)
          coul <- rep(col.upper, length(catdes.aux))
          coul[catdes.aux < 0] <- col.lower
          if (is.null(x$category)) 
            titre <- names(x$quanti)[i]
          else titre <- names(x$category)[i]
          barplot(catdes.aux, width = c(1, 1), col = coul, 
                  border = "black", ylim = c(minimum - 1, maximum + 
                                               1), xlim = c(0, max(long) + 1), main = titre, 
                  cex.names = 1, ylab = "v.test", names.arg = substr(names(catdes.aux), 
                                                                     1, numchar))
        }
      }
    }
    par(las = 0)
  }
}

catdes2 <- function (donnee, num.var, proba = 0.05, row.w = NULL, all.mods=F) 
{
  moy.p <- function(V, fac = NULL, poids, na.rm = TRUE) {
    poids[is.na(V)] <- 0
    if (is.null(fac)) {
      res <- sum(V * poids, na.rm = na.rm)/sum(poids)
    }
    else {
      res = NULL
      for (i in 1:nlevels(fac)) res = c(res, sum(V[fac == 
                                                     levels(fac)[i]] * poids[fac == levels(fac)[i]], 
                                                 na.rm = na.rm)/sum(poids[fac == levels(fac)[i]]))
    }
    return(res)
  }
  ec <- function(V, fac = NULL, poids, na.rm = TRUE) {
    poids[is.na(V)] <- 0
    if (is.null(fac)) {
      V <- V - moy.p(V, fac = NULL, poids, na.rm)
      res <- sum(V^2 * poids, na.rm = na.rm)/sum(poids)
    }
    else {
      moy.par.mod = moy.p(V, fac = fac, poids, na.rm)
      res = NULL
      for (i in 1:nlevels(fac)) {
        res = c(res, sum((V[fac == levels(fac)[i]] - 
                            moy.par.mod[i])^2 * poids[fac == levels(fac)[i]], 
                         na.rm = na.rm)/sum(poids[fac == levels(fac)[i]]))
      }
    }
    return(sqrt(res))
  }
  donnee <- as.data.frame(donnee)
  if (is.numeric(donnee[, num.var])) 
    stop(paste("The variable", num.var, "must be qualitative"))
  is.quali <- which(!unlist(lapply(donnee, is.numeric)))
  donnee[, is.quali] <- lapply(donnee[, is.quali, drop = FALSE], 
                               as.factor)
  donnee <- droplevels(donnee)
  if (is.null(row.w)) 
    row.w = rep(1, nrow(donnee))
  lab.sauv <- lab <- colnames(donnee)
  quali = NULL
  for (i in 1:length(lab)) {
    lab[i] = gsub(" ", ".", lab[i])
    if (is.factor(donnee[, i])) {
      if (any(is.na(donnee[, i]))) {
        levels(donnee[, i]) <- c(levels(donnee[, i]), 
                                 "NA")
        donnee[, i][is.na(donnee[, i])] <- "NA"
      }
      if (levels(donnee[, i])[1] == "") 
        levels(donnee[, i])[1] = "NA"
      if (i != num.var) 
        quali = c(quali, i)
    }
  }
  quanti = (1:ncol(donnee))[-c(quali, num.var)]
  if (length(quanti) == 0) 
    quanti = NULL
  colnames(donnee) = lab
  res = list()
  nb.modalite <- nlevels(donnee[, num.var])
  nb.quali = length(quali)
  old.warn = options("warn")
  if (nb.quali > 0) {
    options(warn = -1)
    Test.chi = matrix(NA, nrow = nb.quali, ncol = 2)
    marge.li = apply(sweep(tab.disjonctif(donnee[, num.var]), 
                           1, row.w, FUN = "*"), 2, sum)
    nom = tri = structure(vector(mode = "list", length = nb.modalite), 
                          names = levels(donnee[, num.var]))
    indicateur.quali <- 0
    for (i in 1:nb.quali) {
      Table <- t(tab.disjonctif(donnee[, num.var])) %*% 
        sweep(tab.disjonctif(donnee[, quali[i]]), 1, 
              row.w, FUN = "*")
      marge.col = apply(Table, 2, sum)
      Test <- chisq.test(Table, correct = FALSE)
      Test.chi[i, 1] <- Test$p.value
      Test.chi[i, 2] <- Test$parameter
      for (j in 1:nlevels(donnee[, num.var])) {
        for (k in 1:nlevels(donnee[, quali[i]])) {
          aux2 = Table[j, k]/marge.li[j]
          aux3 = marge.col[k]/sum(marge.col)
          aux4 <- min(phyper(round(Table[j, k], 0) - 
                               1, round(marge.li[j], 0), round(sum(marge.li), 
                                                               0) - round(marge.li[j], 0), round(marge.col[k], 
                                                                                                 0)) * 2 + dhyper(round(Table[j, k], 0), round(marge.li[j], 
                                                                                                                                               0), round(sum(marge.li), 0) - round(marge.li[j], 
                                                                                                                                                                                   0), round(marge.col[k], 0)), phyper(round(Table[j, 
                                                                                                                                                                                                                                   k], 0), round(marge.li[j], 0), round(sum(marge.li), 
                                                                                                                                                                                                                                                                        0) - round(marge.li[j], 0), round(marge.col[k], 
                                                                                                                                                                                                                                                                                                          0), lower.tail = FALSE) * 2 + dhyper(round(Table[j, 
                                                                                                                                                                                                                                                                                                                                                           k], 0), round(marge.li[j], 0), round(sum(marge.li), 
                                                                                                                                                                                                                                                                                                                                                                                                0) - round(marge.li[j], 0), round(marge.col[k], 
                                                                                                                                                                                                                                                                                                                                                                                                                                  0)))
          if (all.mods | aux4 <= proba) {
            aux5 = (1 - 2 * as.integer(aux2 > aux3)) * 
              qnorm(aux4/2)
            aux1 = Table[j, k]/marge.col[k]
            tri[[j]] = rbind(tri[[j]], c(aux1 * 100, 
                                         aux2 * 100, aux3 * 100, aux4, aux5))
            nom[[j]] = rbind(nom[[j]], c(levels(donnee[, 
                                                       quali[i]])[k], colnames(donnee)[quali[i]]))
          }
        }
      }
      rownames(Test.chi) = colnames(donnee)[quali]
    }
    if (nrow(matrix(Test.chi, ncol = 2)) > 1) {
      if (sum(Test.chi[, 1] <= proba) == 1) {
        nomaux = rownames(Test.chi[order(Test.chi[, 1]), 
        ])[1]
        Test.chi = matrix(Test.chi[Test.chi[, 1] <= proba, 
        ], ncol = 2)
        rownames(Test.chi) = nomaux
      }
      else Test.chi = Test.chi[Test.chi[, 1] <= proba, 
      ]
    }
    else if (Test.chi[, 1] > proba) 
      Test.chi = NULL
    if (!is.null(Test.chi)) {
      if (nrow(matrix(Test.chi, ncol = 2)) > 1) {
        oo = order(Test.chi[, 1])
        Test.chi = Test.chi[oo, ]
      }
      colnames(Test.chi) = c("p.value", "df")
      res$test.chi2 = Test.chi
    }
    for (j in 1:nb.modalite) {
      if (!is.null(tri[[j]])) {
        indicateur.quali <- 1
        oo = rev(order(tri[[j]][, 5]))
        tri[[j]] = tri[[j]][oo, ]
        nom[[j]] = nom[[j]][oo, ]
        if (nrow(matrix(tri[[j]], ncol = 5)) > 1) 
          rownames(tri[[j]]) = paste(nom[[j]][, 2], nom[[j]][, 
                                                             1], sep = "=")
        else {
          tri[[j]] = matrix(tri[[j]], ncol = 5)
          rownames(tri[[j]]) = paste(nom[[j]][2], nom[[j]][1], 
                                     sep = "=")
        }
        colnames(tri[[j]]) = c("Cla/Mod", "Mod/Cla", 
                               "Global", "p.value", "v.test")
      }
    }
    if (indicateur.quali > 0) 
      res$category = tri
  }
  if (!is.null(quanti)) {
    nom = result = structure(vector(mode = "list", length = nb.modalite), 
                             names = levels(donnee[, num.var]))
    tabF <- matrix(0, length(quanti), 2)
    for (i in 1:length(quanti)) {
      res.aov <- summary(aov(donnee[, quanti[i]] ~ donnee[, 
                                                          num.var], na.action = na.exclude, weights = row.w))[[1]]
      tabF[i, 1] <- res.aov[1, 2]/sum(res.aov[, 2])
      tabF[i, 2] <- res.aov[1, 5]
      moy.mod = moy.p(donnee[, quanti[i]], fac = donnee[, 
                                                        num.var], poids = row.w)
      n.mod = apply(sweep(tab.disjonctif(donnee[, num.var]), 
                          1, row.w, FUN = "*"), 2, sum)
      sd.mod = ec(donnee[, quanti[i]], fac = donnee[, num.var], 
                  poids = row.w)
      moy = moy.p(donnee[, quanti[i]], poids = row.w)
      et = ec(donnee[, quanti[i]], poids = row.w)
      for (j in 1:nb.modalite) {
        v.test = (moy.mod[j] - moy)/et * sqrt(n.mod[j])/sqrt((sum(n.mod) - 
                                                                n.mod[j])/(sum(n.mod) - 1))
        p.value = pnorm(abs(v.test), lower.tail = FALSE) * 
          2
        if (!is.na(v.test)) {
          if (all.mods | p.value <= proba) {
            result[[j]] = rbind(result[[j]], c(v.test, 
                                               moy.mod[j], moy, sd.mod[j], et, p.value))
            nom[[j]] = c(nom[[j]], colnames(donnee)[quanti[i]])
          }
        }
      }
    }
    dimnames(tabF) <- list(colnames(donnee)[quanti], c("Eta2", 
                                                       "P-value"))
    auxF <- tabF[order(tabF[, 2]), , drop = FALSE]
    select1 <- (1:nrow(auxF))[auxF[, 2, drop = FALSE] <= 
                                proba]
    if (length(select1) > 0) 
      resF <- auxF[select1, , drop = FALSE]
    for (j in 1:nb.modalite) {
      if (!is.null(result[[j]])) {
        oo = rev(order(result[[j]][, 1]))
        result[[j]] = result[[j]][oo, , drop = FALSE]
        nom[[j]] = nom[[j]][oo]
        result[[j]] = matrix(result[[j]], ncol = 6)
        rownames(result[[j]]) = nom[[j]]
        colnames(result[[j]]) = c("v.test", "Mean in category", 
                                  "Overall mean", "sd in category", "Overall sd", 
                                  "p.value")
      }
    }
    if (length(select1) > 0) {
      res$quanti.var = resF
      res$quanti = result
    }
  }
  res$call = list(num.var = num.var, proba = proba, row.w = row.w, 
                  X = donnee)
  options(old.warn)
  class(res) <- c("catdes", "list ")
  return(res)
}
