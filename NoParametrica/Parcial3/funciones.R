monregcool <- function(x, y, topred = NULL) {
    y <- y[order(x)]; x <- sort(x); n <- length(x)
    rx <- rank(x); ry <- rank(y)
    df <- data.frame(x, y, rx, ry)
    b2 <- (sum(rx * ry) - n * (n + 1)^2 / 4) / (sum(rx^2) - n * (n + 1)^2 / 4)
    a2 <- (1 - b2) * (n + 1) / 2
    ry_hat <- a2 + b2 * rx
    rx_hat <- (ry - a2) / b2
    df$rx_hat <- rx_hat; df$ry_hat <- ry_hat
    y_hat <- c()
    for (i in 1:n) {
        if (ry_hat[i] > max(ry)) {
            y_hat[i] <- max(y)
        } else if (ry_hat[i] < min(ry)) {
            y_hat[i] <- min(y)
        } else if (ry_hat[i] %in% ry) {
            y_hat[i] <- y[which.max(ry_hat[i] == ry)]
        } else {
            temp <- ry - ry_hat[i]
            ryi <- max(temp[temp < 0]) + ry_hat[i]
            ryj <- min(temp[temp > 0]) + ry_hat[i]
            yi <- y[which.max(ryi == ry)]
            yj <- y[which.max(ryj == ry)]
            y_hat[i] <- yi + (yj - yi) * (ry_hat[i] - ryi) / (ryj - ryi)
        }
    }
    df$y_hat <- y_hat

    x_hat <- c()
    for (i in 1:n) {
        if (rx_hat[i] < min(rx) || rx_hat[i] > max(rx)) {
            x_hat[i] <- NA
        } else if (rx_hat[i] %in% rx) {
            x_hat[i] <- x[which.max(rx_hat[i] == rx)]
        } else {
            temp <- rx - rx_hat[i]
            rxi <- max(temp[temp < 0]) + rx_hat[i]
            rxj <- min(temp[temp > 0]) + rx_hat[i]
            xi <- x[which.max(rxi == rx)]
            xj <- x[which.max(rxj == rx)]
            x_hat[i] <- xi + (xj - xi) * (rx_hat[i] - rxi) / (rxj - rxi)
        }
    }
    df$x_hat <- x_hat

    if (! is.null(topred)) {
        npred <- length(topred)
        temp <- c()
        preds <- c()
        for (i in 1:npred) {
            if (topred[i] < min(x) || rx_hat[i] > max(x)) {
                temp[i] <- NA
                preds[i] <- NA
                break
            } else if (topred[i] %in% x) {
                temp[i] <- rx[which.max(topred[i] == x)]
            } else {
                dists <- x - topred[i]
                xi <- max(dists[dists < 0]) + topred[i]
                xj <- min(dists[dists > 0]) + topred[i]
                rxi <- rx[which.max(xi == x)]
                rxj <- rx[which.max(xj == x)]
                temp[i] <- rxi + (rxj - rxi) * (topred[i] - xi) / (xj - xi)
            }
            tempy <- a2 + b2 * temp[i]
            if (tempy > max(ry)) {
                preds[i] <- max(y)
            } else if (tempy < min(ry)) {
                preds[i] <- min(y)
            } else if (tempy %in% ry) {
                preds[i] <- y[which.max(tempy == ry)]
            } else {
            distsy <- ry - tempy
            ryi <- max(distsy[distsy < 0]) + tempy
            ryj <- min(distsy[distsy > 0]) + tempy
            yi <- y[which.max(ryi == ry)]
            yj <- y[which.max(ryj == ry)]
            preds[i] <- yi + (yj - yi) * (tempy - ryi) / (ryj - ryi)
            }
        }
        df2 <- data.frame(x = topred, yhat = preds)
        return(list(main = df, predictions = df2))
    }
    return(list(main = df))
}