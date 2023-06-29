#################### R wrapper of mtool ####################

# mtool.estimate

mtool.estimate <- function(X, Y, weight, N_covariates,
                           loss = 'ls', regularization = 'l1',
                           transpose = F,
                           lambda1, lambda2 = 0, lambda3 = 0,
                           learning_rate = 1e-4, tolerance = 1e-3,
                           niter_inner_mtplyr = 5, maxit = 100, ncores = -1,
                           group_id, group_weights,
                           groups, groups_var,
                           own_variables, N_own_variables) {
    ## Dimensions and checks
    nx <- nrow(X)
    ny <- nrow(Y)
    if (nx == ny) {
        n <- nx
    } else {
        stop('X and Y have different number of observations.')
    }
    
    if (missing(weight)) {
        weight = rep(1, nx)
    }
    
    if (!length(weight) == nx) {
        stop(" 'weight' should be a vector with length equal to the sample size.")
    }
    
    p <- ncol(X)
    
    K <- ncol(Y)
    
    ## loss
    lossfns <- c('ls', 'logistic', 'cox')
    
    if (loss %in% lossfns) {
        lossfunc <- which(lossfns == loss)
    } else {
        stop('The provided loss function is not supported.')
    }
    
    ## regularization
    pen1 <- c("l0", "l1", "l2", "linf", "l2-not-squared",
              "elastic-net", "fused-lasso",
              "group-lasso-l2", "group-lasso-linf",
              "sparse-group-lasso-l2", "sparse-group-lasso-linf",
              "l1l2", "l1linf", "l1l2+l1", "l1linf+l1", "l1linf-row-column",
              "trace-norm", "trace-norm-vec", "rank", "rank-vec", "none")
    pen2 <- c("graph", "graph-ridge", "graph-l2", "multi-task-graph")
    pen3 <- c("tree-l0", "tree-l2", "tree-linf", "multi-task-tree")
    
    if (regularization %in% pen1) { penalty <- 1 }
    if (regularization %in% pen2) { penalty <- 2 }
    if (regularization %in% pen3) { penalty <- 3 }
    if (! regularization %in% c(pen1, pen2, pen3)) {
        stop('The provided regularization is not supported.')
    }
    
    ### check regularization-specific inputs
    #### penalty = 1, call proximal(Flat), requires `group_id` in integer vector
    if (penalty == 1) {
        if (missing(group_id)) { group_id <- rep(0L, p) }
        group_weights <- vector(mode = 'double')
        groups <- matrix(NA)
        groups_var <- matrix(NA)
        own_variables <- vector(mode = 'integer')
        N_own_variables <- vector(mode = 'integer')
    }
    
    #### penalty = 2, call proximalGraph
    #### requires `groups` and `groups_var` in integer matrices and `group_weights` in double vector
    if (penalty == 2) {
        if (missing(groups)) { stop('Required input `groups` is missing.') }
        if (missing(groups_var)) { stop('Required input `groups_var` is missing.') }
        if (missing(group_weights)) { stop('Required input `group_weights` is missing.') }
        group_id <- rep(0L, p)
        own_variables <- vector(mode = 'integer')
        N_own_variables <- vector(mode = 'integer')
    }
    
    #### penalty = 3, call proximalGraph
    #### requires `own_variables` and `N_own_variables` in integer vectors, `group_weights` in double vector
    #### and `groups` in integer matrix
    if (penalty == 3) {
        if (missing(groups)) { stop('Required input `groups` is missing.') }
        if (missing(own_variables)) { stop('Required input `own_variables` is missing.') }
        if (missing(N_own_variables)) { stop('Required input `N_own_variables` is missing.') }
        if (missing(group_weights)) { stop('Required input `group_weights` is missing.') }
        group_id <- rep(0L, p)
        groups_var <- matrix(NA)
    }
    
    ## extract number of observations per task and indices of the observations
    obs.flag <- !is.na(Y)
    n_k <- apply(obs.flag, MARGIN = 2, FUN = sum)
    ### extract the indices of non-missing observations for each task
    task_rowid <- apply(obs.flag, MARGIN = 2, FUN = which)
    #### if missingness pattern is identical across tasks, task_rowid will be a matrix instead of a list
    if (is.list(task_rowid)) {
        task_rowid <- unlist(task_rowid)
    } else {
        if (is.matrix(task_rowid)) {
            task_rowid <- as.vector(task_rowid)
        }
    }
    
    ## call mtool main function
    result <- mtool(X = X, Y = Y, wt = weight, K = K, nk_vec = n_k,
                    task_rowid = task_rowid, reg_p = p - N_covariates,
                    loss = lossfunc, penalty = penalty, regul = regularization,
                    transpose = transpose,
                    grp_id = group_id, etaG = group_weights,
                    grp = groups, grpV = groups_var,
                    own_var = own_variables, N_own_var = N_own_variables,
                    lam1 = lambda1, lam2 = lambda2, lam3 = lambda3,
                    learning_rate = learning_rate, tolerance = tolerance,
                    niter_inner = niter_inner_mtplyr * nx, maxit = maxit,
                    ncores = ncores)
    nzc <- length(result$`Sparse Estimates`@i)
    return(list(coefficients = result$`Sparse Estimates`,
                no_non_zero = nzc))
}






########## include a given random fold id assignment

# cross-validation
mtool.cv_error <- function(X, Y, weight, nfold, fold_seed,
                           CV_criterion,
                           N_covariates,
                           loss = 'ls', regularization = 'l1',
                           transpose = F,
                           lambda1, lambda2 = 0, lambda3 = 0,
                           learning_rate = 1e-4, tolerance = 1e-3,
                           niter_inner_mtplyr = 5, maxit = 100, ncores = 1,
                           group_id, group_weights,
                           groups, groups_var,
                           own_variables, N_own_variables) {
    ## Dimensions and checks
    nx <- nrow(X)
    ny <- nrow(Y)
    if (nx == ny) {
        n <- nx
    } else {
        stop('X and Y have different number of observations.')
    }
    
    if (ncores != 1) {
        cat('More than 1 core is being used for the computation of the proximal operator. Be careful if the cross validation procedure is running on parallel.')
    }
    
    if (missing(weight)) {
        weight = rep(1, nx)
    }
    
    if (!length(weight) == nx) {
        stop(" 'weight' should be a vector with length equal to the sample size.")
    }
    
    p <- ncol(X)
    
    K <- ncol(Y)
    
    ## loss
    lossfns <- c('ls', 'logistic', 'cox')
    
    if (loss %in% lossfns) {
        lossfunc <- which(lossfns == loss)
    } else {
        stop('The provided loss function is not supported.')
    }
    
    ## regularization
    pen1 <- c("l0", "l1", "l2", "linf", "l2-not-squared",
              "elastic-net", "fused-lasso",
              "group-lasso-l2", "group-lasso-linf",
              "sparse-group-lasso-l2", "sparse-group-lasso-linf",
              "l1l2", "l1linf", "l1l2+l1", "l1linf+l1", "l1linf-row-column",
              "trace-norm", "trace-norm-vec", "rank", "rank-vec", "none")
    pen2 <- c("graph", "graph-ridge", "graph-l2", "multi-task-graph")
    pen3 <- c("tree-l0", "tree-l2", "tree-linf", "multi-task-tree")
    
    if (regularization %in% pen1) { penalty <- 1 }
    if (regularization %in% pen2) { penalty <- 2 }
    if (regularization %in% pen3) { penalty <- 3 }
    if (! regularization %in% c(pen1, pen2, pen3)) {
        stop('The provided regularization is not supported.')
    }
    
    ### check regularization-specific inputs
    #### penalty = 1, call proximal(Flat), requires `group_id` in integer vector
    if (penalty == 1) {
        if (missing(group_id)) { group_id <- rep(0L, p) }
        group_weights <- vector(mode = 'double')
        groups <- matrix(NA)
        groups_var <- matrix(NA)
        own_variables <- vector(mode = 'integer')
        N_own_variables <- vector(mode = 'integer')
    }
    
    #### penalty = 2, call proximalGraph
    #### requires `groups` and `groups_var` in integer matrices and `group_weights` in double vector
    if (penalty == 2) {
        if (missing(groups)) { stop('Required input `groups` is missing.') }
        if (missing(groups_var)) { stop('Required input `groups_var` is missing.') }
        if (missing(group_weights)) { stop('Required input `group_weights` is missing.') }
        group_id <- rep(0L, p)
        own_variables <- vector(mode = 'integer')
        N_own_variables <- vector(mode = 'integer')
    }
    
    #### penalty = 3, call proximalGraph
    #### requires `own_variables` and `N_own_variables` in integer vectors, `group_weights` in double vector
    #### and `groups` in integer matrix
    if (penalty == 3) {
        if (missing(groups)) { stop('Required input `groups` is missing.') }
        if (missing(own_variables)) { stop('Required input `own_variables` is missing.') }
        if (missing(N_own_variables)) { stop('Required input `N_own_variables` is missing.') }
        if (missing(group_weights)) { stop('Required input `group_weights` is missing.') }
        group_id <- rep(0L, p)
        groups_var <- matrix(NA)
    }
    
    ## extract number of observations per task and indices of the observations
    obs.flag <- !is.na(Y)
    n_k <- apply(obs.flag, MARGIN = 2, FUN = sum)
    task_rowid <- apply(obs.flag, MARGIN = 2, FUN = which)
    #### if missingness pattern is identical across tasks, task_rowid will be a matrix instead of a list
    if (is.matrix(task_rowid)) {
        task_rowid <- as.list(as.data.frame(task_rowid))
    }
    set.seed(fold_seed)
    foldid <- lapply(task_rowid, FUN = function(x) {sample(rep(seq(nfold), length = length(x)))})
    set.seed(NULL)
    
    # cat("fold_ids", utils::head(foldid[[5]],10), "random", rnorm(1))
    
    ## cross validation
    ntrain <- nx * (1 - 1/nfold)
    pred.mat <- matrix(NA, nrow = nx, ncol = K)
    err.fold <- numeric(nfold)
    
    for (i in seq(nfold)) {
        cat("\n fold", i)
        
        ## mapply() works on lists and arrays (vectors, matrices, etc.)
        task_rowid_train <- mapply(task_rowid, foldid, FUN = function(x, y) {x[y != i ]}, SIMPLIFY = F)
        task_rowid_test <- mapply(task_rowid, foldid, FUN = function(x, y) {x[y == i ]}, SIMPLIFY = F)
        
        n_k <- unlist(lapply(task_rowid_train, FUN = length))
        
        fit <-  mtool(X = X, Y = Y, wt = weight, K = K, nk_vec = n_k,
                      task_rowid = unlist(task_rowid_train), reg_p = p - N_covariates,
                      loss = lossfunc, penalty = penalty, regul = regularization,
                      transpose = transpose,
                      grp_id = group_id, etaG = group_weights,
                      grp = groups, grpV = groups_var,
                      own_var = own_variables, N_own_var = N_own_variables,
                      lam1 = lambda1, lam2 = lambda2, lam3 = lambda3,
                      learning_rate = learning_rate, tolerance = tolerance,
                      niter_inner = niter_inner_mtplyr * ntrain,
                      maxit = maxit, ncores = ncores)
        beta.list <- split(fit$`Sparse Estimates`, f = rep(1:K, each = p))
        
        pred <- mapply(task_rowid_test, beta.list, FUN = function(x, y) {X[x,] %*% y}, SIMPLIFY = F)
        
        pred.mat.fold <- matrix(NA, nrow = nx, ncol = K)
        
        for (k in seq(K)) {
            pred.mat[task_rowid_test[[k]], k] <- pred[[k]]
            pred.mat.fold[task_rowid_test[[k]], k] <- pred[[k]]
        }
        
        ### cv error for each fold to calculate var(cv.error)
        if (CV_criterion == "MSE") {
            error.mat.fold <- (pred.mat.fold - Y)^2
            err.fold[i] <- sum(apply(error.mat.fold, MARGIN = 2, FUN = mean, na.rm = T))
        }
        
        if (CV_criterion == "Pearson") {
            err.fold[i] <- mean(diag(stats::cor(pred.mat, Y, use = "pairwise.complete.obs", method = "pearson")))
        }
    }
    
    # confirm missingness patterns in pred.mat and y are the same
    ## Somehow there can be false alarms.
    missing.check <- identical(is.na(Y), is.na(pred.mat))
    if (!missing.check) {
        warning("Something may be wrong with cross-validation. Missingness in CV prediction and in original Y may be different.")
    }
    
    if (CV_criterion == "MSE") {
        error.mat <- (pred.mat - Y)^2
        # print(sum(apply(error.mat, MARGIN = 2, FUN = mean, na.rm = T)))
        error <- sum(apply(error.mat, MARGIN = 2, FUN = mean, na.rm = T))
    }
    
    if (CV_criterion == "Pearson") {
        error <- mean(diag(stats::cor(pred.mat, Y, use = "pairwise.complete.obs", method = "pearson")))
    }
    
    var.error <- 1/K * stats::var(err.fold)
    
    return(list(CV_prediction = pred.mat,
                CV_error = error,
                CV_error_var = var.error))
}


# Automatic tree builder with manual input
## input: a list of tree nodes
## equally spaced node height
## output: a list of own_variables, N_own_variables, groups, eta_g, as required by mtool.estimate and others.


tree_builder <- function(tree) {
    n.nodes <- length(tree)
    
    tree <- lapply(tree, sort)
    
    reorder1 <- order(sapply(tree, FUN = length), decreasing = T)
    
    tree <- tree[reorder1]
    
    reorder2 <- order(sapply(tree, FUN = utils::head, n = 1))
    
    tree <- tree[reorder2]
    
    # pre-allocate own_variable, N_own_variable, group weights, groups matrix
    n_own_var <- integer(n.nodes)
    eta_g <- numeric(n.nodes)
    g <- matrix(0L, nrow = n.nodes, ncol = n.nodes)
    
    # own_variables
    own_var <- sapply(tree, FUN = utils::head, n = 1) - 1L
    
    # N_own_variables
    for (i in seq(n.nodes)) {
        nodes.set <- tree[-(1:i)]
        node.tmp <- tree[[i]]
        nov.tmp <- 0L
        for (j in seq(length(node.tmp))) {
            if (!node.tmp[j] %in% unlist(nodes.set)) {
                nov.tmp <- nov.tmp + 1L
            }
        }
        n_own_var[i] <- nov.tmp
    }
    
    # groups
    for (i in 2:n.nodes) {
        for (j in (i-1):1) {
            if ( setequal(intersect(tree[[i]], tree[[j]]),  tree[[i]]) ) {
                g[i, j] <- 1L
                break
            }
        }
    }
    
    # group weights: eta_g
    
    w <- pass.down <- numeric(n.nodes)
    pass.down[1] <- 1
    w[1] <- 0
    
    sub.tree.old <- tree
    
    for (i in 2:n.nodes) {
        node <- tree[[i]]
        
        ancester <- which(g[i,] == 1L)
        
        sub.tree <- tree[sapply(tree, FUN = function(x) {setequal(x, intersect(x, node))})]
        sub.tree.ancester <- sub.tree
        sub.tree.ancester[[length(sub.tree.ancester) + 1]] <- tree[[ancester]]
        
        ancester.height <- max(table(unlist(sub.tree.ancester)))
        node.height <- max(table(unlist(sub.tree)))
        pass.down[i] <- (node.height - 1) / (ancester.height - 1) * pass.down[ancester]
        
        if (ancester == 1) {
            w[i] <- 1 - pass.down[i]
        } else {
            w[i] <- w[ancester] / (1 - pass.down[ancester]) * pass.down[ancester] * (1 - pass.down[i])
        }
    }
    return(list(own_variables = own_var,
                N_own_variables = n_own_var,
                groups = g,
                eta_g = w))
    
}








# Automatic tree builder with hclust() input
## input: an hclust object
## height generated by hclust()
## output: a list of own_variables, N_own_variables, groups, eta_g, as required by mtool.estimate and others.

tree_builder_hclust <- function(cluster) {
    merge <- cluster$merge
    n.nonleafnodes <- nrow(merge)
    n.leafnodes <- - min(merge)
    n.nodes <- n.nonleafnodes + n.leafnodes
    
    tree <- vector(mode = 'list', length = n.nodes)
    height <- numeric(n.nodes)
    
    for (i in seq(n.nonleafnodes)) {
        node <- merge[i,]
        if (node[1] < 0 & node[2] < 0) {
            tree[[i]] <- - node
        } else {
            if (node[1] <0 | node[2] < 0) {
                leaf <- - min(node)
                prev.node <- max(node)
                tree[[i]] <- c(leaf, tree[[prev.node]])
            } else {
                tree[[i]] <- c(tree[[ node[1] ]], tree[[ node[2] ]])
            }
        }
    }
    
    leaves.indices <- n.nonleafnodes + 1:n.leafnodes
    
    for (i in leaves.indices) {
        tree[[i]] <- i - n.nonleafnodes
    }
    
    height <- c(cluster$height, rep(0, n.leafnodes))
    
    tree <- lapply(tree, sort)
    
    reorder1 <- order(sapply(tree, FUN = length), decreasing = T)
    
    tree <- tree[reorder1]
    height <- height[reorder1]
    
    reorder2 <- order(sapply(tree, FUN = utils::head, n = 1))
    
    tree <- tree[reorder2]
    height <- height[reorder2]
    height <- height / max(height)
    
    # pre-allocate own_variable, N_own_variable, group weights, groups matrix
    n_own_var <- integer(n.nodes)
    eta_g <- numeric(n.nodes)
    g <- matrix(0L, nrow = n.nodes, ncol = n.nodes)
    
    # own_variables
    own_var <- sapply(tree, FUN = utils::head, n = 1) - 1L
    
    # N_own_variables
    for (i in seq(n.nodes)) {
        nodes.set <- tree[-(1:i)]
        node.tmp <- tree[[i]]
        nov.tmp <- 0L
        for (j in seq(length(node.tmp))) {
            if (!node.tmp[j] %in% unlist(nodes.set)) {
                nov.tmp <- nov.tmp + 1L
            }
        }
        n_own_var[i] <- nov.tmp
    }
    
    # groups
    for (i in 2:n.nodes) {
        for (j in (i-1):1) {
            if ( setequal(intersect(tree[[i]], tree[[j]]),  tree[[i]]) ) {
                g[i, j] <- 1L
                break
            }
        }
    }
    
    # group weights: eta_g
    
    w <- numeric(n.nodes)
    w[1] <- 0
    
    sub.tree.old <- tree
    
    for (i in 2:n.nodes) {
        
        ancestor <- which(g[i,] == 1L)
        
        # sub.tree <- nodes[sapply(nodes, FUN = function(x) {setequal(x, intersect(x, node))})]
        # sub.tree.ancestor <- sub.tree
        # sub.tree.ancestor[[length(sub.tree.ancestor) + 1]] <- nodes[[ancestor]]
        # 
        # ancestor.height <- max(table(unlist(sub.tree.ancestor)))
        # node.height <- max(table(unlist(sub.tree)))
        # pass.down[i] <- (node.height - 1) / (ancestor.height - 1) * pass.down[ancestor]
        
        if (ancestor == 1) {
            w[i] <- 1 - height[i]
        } else {
            w[i] <- w[ancestor] / (1 - height[ancestor]) * height[ancestor] * (1 - height[i])
        }
    }
    
    return(list(own_variables = own_var,
                N_own_variables = n_own_var,
                groups = g,
                eta_g = w,
                tree = tree))
}






# Penalized Multinomial Logistic Regression
mtool.MNlogistic <- function(X, Y, offset, N_covariates,
                             regularization = 'l1', transpose = F,
                             lambda1, lambda2 = 0, lambda3 = 0,
                             learning_rate = 1e-4, tolerance = 1e-4,
                             niter_inner_mtplyr = 7, maxit = 100, ncores = -1,
                             group_id, group_weights,
                             groups, groups_var,
                             own_variables, N_own_variables) {
    ## Dimensions and checks
    nx <- nrow(X)

    if (!is.vector(Y)) {Y <- as.vector(Y)}
    ny <- length(Y)

    if (!is.vector(offset)) {offset <- as.vector(offset)}
    noff <- length(offset)

    if (nx == ny & nx == noff) {
        n <- nx
    } else {
        stop('X, Y and offset have different number of observations.')
    }

    p <- ncol(X)

    K <- length(unique(Y)) - 1

    ## regularization
    pen1 <- c("l0", "l1", "l2", "linf", "l2-not-squared",
              "elastic-net", "fused-lasso",
              "group-lasso-l2", "group-lasso-linf",
              "sparse-group-lasso-l2", "sparse-group-lasso-linf",
              "l1l2", "l1linf", "l1l2+l1", "l1linf+l1", "l1linf-row-column",
              "trace-norm", "trace-norm-vec", "rank", "rank-vec", "none")
    pen2 <- c("graph", "graph-ridge", "graph-l2", "multi-task-graph")
    pen3 <- c("tree-l0", "tree-l2", "tree-linf", "multi-task-tree")

    if (regularization %in% pen1) { penalty <- 1 }
    if (regularization %in% pen2) { penalty <- 2 }
    if (regularization %in% pen3) { penalty <- 3 }
    if (! regularization %in% c(pen1, pen2, pen3)) {
        stop('The provided regularization is not supported.')
    }

    ### check regularization-specific inputs
    #### penalty = 1, call proximal(Flat), requires `group_id` in integer vector
    if (penalty == 1) {
        if (missing(group_id)) { group_id <- rep(0L, p) }
        group_weights <- vector(mode = 'double')
        groups <- matrix(NA)
        groups_var <- matrix(NA)
        own_variables <- vector(mode = 'integer')
        N_own_variables <- vector(mode = 'integer')
    }

    #### penalty = 2, call proximalGraph
    #### requires `groups` and `groups_var` in integer matrices and `group_weights` in double vector
    if (penalty == 2) {
        if (missing(groups)) { stop('Required input `groups` is missing.') }
        if (missing(groups_var)) { stop('Required input `groups_var` is missing.') }
        if (missing(group_weights)) { stop('Required input `group_weights` is missing.') }
        group_id <- rep(0L, p)
        own_variables <- vector(mode = 'integer')
        N_own_variables <- vector(mode = 'integer')
    }

    #### penalty = 3, call proximalGraph
    #### requires `own_variables` and `N_own_variables` in integer vectors, `group_weights` in double vector
    #### and `groups` in integer matrix
    if (penalty == 3) {
        if (missing(groups)) { stop('Required input `groups` is missing.') }
        if (missing(own_variables)) { stop('Required input `own_variables` is missing.') }
        if (missing(N_own_variables)) { stop('Required input `N_own_variables` is missing.') }
        if (missing(group_weights)) { stop('Required input `group_weights` is missing.') }
        group_id <- rep(0L, p)
        groups_var <- matrix(NA)
    }

    ## call mtool main function
    result <- MultinomLogistic(X = X, Y = Y, offset = offset, K = K, reg_p = p - N_covariates,
                               penalty = penalty, regul = regularization, transpose = transpose,
                               grp_id = group_id, etaG = group_weights,
                               grp = groups, grpV = groups_var,
                               own_var = own_variables, N_own_var = N_own_variables,
                               lam1 = lambda1, lam2 = lambda2, lam3 = lambda3,
                               learning_rate = learning_rate, tolerance = tolerance,
                               niter_inner = niter_inner_mtplyr * nx, maxit = maxit,
                               ncores = ncores)
    nzc <- length(result$`Sparse Estimates`@i)
    return(list(coefficients = result$`Sparse Estimates`,
                no_non_zero = nzc))
}



