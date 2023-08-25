# Title: Simulations/generate_data.R
# Author: Caleb Leedy
# Date Created: April 25, 2023
# Purpose: This file contains the code to generate data for different simulation
# studies. Instead of writing the functions in Simulations/nonmonotone.R and
# then copying them to other files, I thought that it would be a good idea to
    # p2 = Pr(R2 = 1) # NOTE: This is NOT true because R2 = 1 if p12 works.
    # p12 = Pr(R2 = 1 | R1 == 1)
    # p21 = Pr(R1 = 1 | R2 == 1)
    # Instead, p0, p1, and p2 are Pr(i in Group(r)) which determines the order
    # which we fill missing values. Thus,
    # Pr(R1 = 1) = p1 + p21 * p2
    # Pr(R2 = 1) = p2 + p12 * p1

    if (miss_resp) {

      p0 <- abs(x_vec^2 - 1)
      p1 <- abs(x_vec^2 - 2)
      p2 <- abs(x_vec^2)

      tot <- p0 + p1 + p2
      p0 <- p0 / tot
      p1 <- p1 / tot
      p2 <- p2 / tot

      p12 <- expit(y1^2)
      p21 <- expit(y2^2)


    } else {

      p0 <- 0.2
      # p1 <- expit(2 * x_vec)
      p1 <- 0.4
      # p2 <- expit(-2 * x_vec)
      p2 <- 0.4
      tot <- p0 + p1 + p2
      p0 <- p0 / tot
      p1 <- p1 / tot
      p2 <- p2 / tot

      p12 <- expit(y1)
      p21 <- expit(y2)

    }

    rand1 <- runif(n)
    rand2 <- runif(n)

  }

  # 3. Reveal hidden variables.
  tibble(x = x_vec, y1 = y1, y2 = y2,
         p0 = p0, p1 = p1, p2 = p2,
         p12 = p12, p21 = p21,
         rand1 = rand1, rand2 = rand2) %>%
  mutate(out = case_when(
                rand1 < p0 ~ "00",
                rand1 < p0 + p1 & rand2 < 1 - p12 ~ "10",
                rand1 < p0 + p1 & rand2 > 1 - p12 ~ "11a",
                rand1 < p0 + p1 + p2 & rand2 < 1 - p21 ~ "01",
                rand1 < p0 + p1 + p2 & rand2 > 1 - p21 ~ "11b")) %>%
  mutate(r1 = ifelse(out %in% c("10", "11a", "11b"), 1, 0)) %>%
  mutate(r2 = ifelse(out %in% c("01", "11a", "11b"), 1, 0)) %>%
  dplyr::select(-rand1, -rand2)

}

