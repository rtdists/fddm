devtools::load_all(recompile = TRUE)
devtools::load_all()
library("fddm")

data <- read.csv("unused_stuff/ticket_price_data.csv")
data <- data[ !is.na(data[["response"]]), ]
data <- data[ data[["rt"]] >= 0.25 & data[["rt"]] <= 3, ]
data$is_first <- factor(ifelse(data$pos == 1, "first", "other"),
                        c("first", "other"))
data$response = data$response + 1
data$price_p <- pnorm(data$prices, 180, 20)
data$price_q <- (data$prices - 180) / 20

m1 <- ddm(drift = rt + response ~ price_q,
          boundary = ~pos:is_first,
          bias = ~pos:is_first,
          ndt = ~pos:is_first,
          data = data)
