library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("cpp_files/onlineblock.cpp")

### some other data ####
fb_wall <- read_tsv(gzfile("data/facebook-wall.txt.gz"))

colnames(fb_wall) <- c("send","rec","time_stamp")

head(fb_wall)

fb_wall %>% mutate(time_stamp = lubridate::as_datetime(time_stamp)) %>%
  group_by(send,rec) %>% 
  tally() %>%
  arrange(desc(n))

### could compare this clustering with known adjacency matrix, but is that
### sort of cheating?


yt <- read_tsv(gzfile("data/youtube-links.txt.gz"))




#### Columbia Email Data ####

per_info <- read_lines(gzfile("C:/Users/owenw/Google Drive/Tian/EmailNetwork/Email-network CU/02/data0002.gz"),
                     skip_empty_rows = FALSE)

per_info <- as.data.frame(per_info)

nrow(per_info)/19

vars <- c("uni",
          "is_faculty",
          "is_student",
          "is_officer",
          "is_staff",
          "is_affiliate",
          "birth_year",
          "campus",
          "aca_field",
          "aca_department",
          "adm_department",
          "division",
          "dorm_building",
          "from_us",
          "gender",
          "postal_code",
          "school",
          "school_year",
          "student_status")


per_data <- tibble(value = per_info$per_info,variable = rep(vars,nrow(per_info)/19))

per_data %>% tail(n=20)

# this is correctly capturing everything up to here at least
unis <- per_data %>% filter(variable == "uni") %>% dplyr::select(value)

per_data$id <- rep(as.vector(unis$value),each = 19)

tidy_info <- per_data %>% pivot_wider(names_from = variable,values_from = value,
                                      values_fill = NA)

tidy_info %>% 
  select(-id) %>%
  filter(is_faculty == 1 | is_officer == 1 | is_staff == 1) %>%
  group_by(adm_department) %>% 
  count() %>% 
  arrange(desc(n))


faculty_staff <- tidy_info %>%
  select(-id) %>%
  filter(is_faculty == 1 | is_officer == 1 | is_staff == 1) %>%
  select(uni,adm_department)# %>% # could have a category for those with
  ## no dept also
  # filter(adm_department != "")


faculty_staff
## we could use these known administrative departments as truth ##
## then extract some of the emails

emails <- read_lines(gzfile("C:/Users/owenw/Google Drive/Tian/EmailNetwork/Email-network CU/02/day001.gz"),
                     skip_empty_rows = FALSE)

emails <- emails %>% str_split(pattern = ":")

# now a list entry for each element

emails_tidy <- emails %>% purrr::map(function(.x) tibble(time =.x[1],
                           send = .x[2],
                           rec = .x[3],
                           size = .x[4]))

do.call(emails_tidy,bind_rows)

em <- data.frame(emails_tidy)

em <- do.call(rbind,emails_tidy)

# need to remember how to split these

em %>%
  head() %>%
  stringr::str_split(rec,pattern= ",")

stringr::str_split(em$rec,pattern = ",")


# then need to back this into the data we have