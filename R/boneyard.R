


# as <- unique(c(Ns$a - 0.5, Ns$a + 0.5))
# ts <- unique(c(Ns$t - 0.5, Ns$t + 0.5))
# ggplot(Ns, aes(x = t, y = a)) +
#   geom_tile(aes(fill = Nt)) +
#   geom_text(
#     aes(label = TeX(N, output = "character")),
#     parse = TRUE,
#     size = 7/.pt
#   ) +
#   geom_hline(yintercept = as) +
#   geom_vline(xintercept = ts) +
#   scale_fill_viridis_c(option = "H", alpha = 0.6) +
#   theme_void() +
#   coord_cartesian(expand = FALSE)




# ggplot(Rs, aes(x = t, y = a)) +
#   geom_tile(aes(fill = R)) +
#   geom_text(
#     aes(label = TeX(Rlab, output = "character")),
#     parse = TRUE,
#     size = 7/.pt
#   ) +
#   geom_hline(yintercept = as) +
#   geom_vline(xintercept = ts) +
#   scale_fill_viridis_c(option = "H", alpha = 0.6) +
#   theme_void() +
#   coord_cartesian(expand = FALSE)



# ggplot(Es, aes(x = t, y = a)) +
#   geom_tile(aes(fill = Eft)) +
#   geom_text(
#     aes(label = TeX(Ef_str, output = "character")),
#     parse = TRUE,
#     size = 7/.pt
#   ) +
#   geom_hline(yintercept = as) +
#   geom_vline(xintercept = ts) +
#   scale_fill_viridis_c(option = "H", alpha = 0.6, na.value = NA, limits = c(min_erro, max_erro)) +
#   theme_void() +
#   coord_cartesian(expand = FALSE)



# ggplot(Es, aes(x = t, y = a)) +
#   geom_tile(aes(fill = Emt)) +
#   geom_text(
#     aes(label = TeX(Em_str, output = "character")),
#     parse = TRUE,
#     size = 7/.pt
#   ) +
#   geom_hline(yintercept = as) +
#   geom_vline(xintercept = ts) +
#   scale_fill_viridis_c(option = "H", alpha = 0.6, na.value = NA, limits = c(min_erro, max_erro)) +
#   theme_void() +
#   coord_cartesian(expand = FALSE)


# Es %>%
#   mutate(Eft = ifelse(t %in% c(4, 6), Eft, NA)) %>%
# ggplot(aes(x = t, y = a)) +
#   geom_tile(aes(fill = Eft)) +
#   geom_text(
#     aes(label = TeX(Ef_str, output = "character")),
#     parse = TRUE,
#     size = 7/.pt
#   ) +
#   geom_hline(yintercept = as) +
#   geom_vline(xintercept = ts) +
#   scale_fill_viridis_c(option = "H", alpha = 0.6, na.value = NA, limits = c(min_erro, max_erro)) +
#   theme_void() +
#   coord_cartesian(expand = FALSE)



# Es %>%
#   mutate(Eft = ifelse((t == 4 & a == 5) | (t == 6 & a == 7), Eft, NA)) %>%
# ggplot(aes(x = t, y = a)) +
#   geom_tile(aes(fill = Eft)) +
#   geom_text(
#     aes(label = TeX(Ef_str, output = "character")),
#     parse = TRUE,
#     size = 7/.pt
#   ) +
#   geom_hline(yintercept = as) +
#   geom_vline(xintercept = ts) +
#   scale_fill_viridis_c(option = "H", alpha = 0.6, na.value = NA, limits = c(min_erro, max_erro)) +
#   theme_void() +
#   coord_cartesian(expand = FALSE)


# first, we are going to expand our tibble with the things that
# we need for drawing arrows from all the cells of different
# ages at time tt, to all the cells they can reach in db years
#' @param Es the current tibble of life table stuff
#' @param tt that starting time of the arrows
#' @param db the difference in birth years
#' @param fr age of first reproduction (i.e. 3 for females, 1 for males)
#' @param astar if not NA, then filter out the arrow data except for
#' the ones going up from astar.
cum_surv_probs <- function(Es, tt = 4, db = 2, fr = 3, astar = NA) {
  sp <- distinct(Es, a, S_a) %>% pull(S_a)  # the survival probs
  ma <- max(Es$a)  # max age
  cp <- function(t, a) { # this function computes the cumulative prob of surviving
    # to age a by time t, given that you were alive at time tt
    ret <- NA
    if(t >= tt && t <= tt + db && a - (t - tt) >= fr && a + db - (t - tt) <= ma) {
      astart = a - (t - tt);
      if(astart == a) {
        ret <- 1.0
      } else {
        ret <- prod(sp[(astart + 1):a])
      }
    }
    ret
  }
  Es %>%
    mutate(
      cp_S = map2_dbl(.x = t, .y = a, .f = cp)  # cp_S is the cumulative product of survival
    ) %>%
    mutate(
      cp_S = case_when(
        is.na(astar) | a - (t - tt) == astar ~ cp_S,
        TRUE ~ NA_real_
      ),
      arrow_start = t == tt,
      arrow_end = t == tt + db,
      ae_na = ifelse(arrow_end, 1, NA),
      ae_na = ifelse(a < fr + db, NA, ae_na)  # make sure they don't appear in the cells that don't produce any
    )
}

Esastar <- cum_surv_probs(Es, tt = 4, db = 2, fr = 3, astar = 5)


first_surv_seg_fig <- Esastar %>%
  mutate(
    Eft = ifelse((t == 4 & a == 5) | (t == 6 & a == 7), Eft, NA),
    ae_na = ifelse((t == 4 & a == 5) | (t == 6 & a == 7), ae_na, NA)
  ) %>%
  ggplot(aes(x = t, y = a)) +
  geom_tile(aes(fill = Eft)) +
  geom_text(
    aes(label = TeX(Ef_str, output = "character")),
    parse = TRUE,
    size = 7/.pt
  ) +
  geom_hline(yintercept = as) +
  geom_vline(xintercept = ts) +
  scale_fill_viridis_c(option = "H", alpha = 0.6, na.value = NA, limits = c(min_erro, max_erro)) +
  theme_void() +
  coord_cartesian(expand = FALSE) +
  scale_linewidth_continuous(limits = c(0, 1), range = c(0, 3)) +
  guides(linewidth = "none") +
  geom_segment(
    mapping = aes(
      x = t - 0.5 * !arrow_start,
      xend = t + 0.5 * !arrow_end,
      y = a - 0.5 * !arrow_start,
      yend = a + 0.5 * !arrow_end,
      linewidth = cp_S
    ),
    colour = "gray40",
    alpha = 0.5
  )

first_surv_seg_fig














Esa_m <- cum_surv_probs(Es, tt = 5, db = 6, fr = 1)

Esa_mstar <- Esa_m %>%
  mutate(
    Emt = ifelse((t == 5 & a <= 2) | (t == 11 & a %in% c(7, 8)), Emt, NA),
    ae_na = ifelse((t == 5 & a <= 2) | (t == 11 & a %in% c(7, 8)), ae_na, NA)
  )

ggplot(Esa_mstar, aes(x = t, y = a)) +
  geom_tile(aes(fill = Emt)) +
  geom_text(
    aes(label = TeX(Em_str, output = "character")),
    parse = TRUE,
    size = 7/.pt
  ) +
  geom_hline(yintercept = as) +
  geom_vline(xintercept = ts) +
  scale_fill_viridis_c(option = "H", alpha = 0.6, na.value = NA, limits = c(min_erro, max_erro)) +
  theme_void() +
  coord_cartesian(expand = FALSE) +
  scale_linewidth_continuous(limits = c(0, 1), range = c(0, 3)) +
  guides(linewidth = "none") +
  geom_segment(
    mapping = aes(
      x = t - 0.5 * !arrow_start,
      xend = t + 0.5 * !arrow_end,
      y = a - 0.5 * !arrow_start,
      yend = a + 0.5 * !arrow_end,
      linewidth = cp_S
    ),
    colour = "gray40",
    alpha = 0.5
  )








first_surv_seg_fig +
  geom_point(aes(x = ae_na * t + 0.3, y = ae_na * a + 0.3), shape = 21, fill = NA, stroke = 0.2)












Es %>%
  mutate(Eft = ifelse(t %in% c(4, 6), Eft, NA)) %>%
  mutate(Eft = ifelse((t == 4 & a >= 7) | (t==6 & a < 5), NA, Eft)) %>%
  cum_surv_probs(tt = 4, db = 2, fr = 3) %>%
  ggplot(aes(x = t, y = a)) +
  geom_tile(aes(fill = Eft)) +
  geom_text(
    aes(label = TeX(Ef_str, output = "character")),
    parse = TRUE,
    size = 7/.pt
  ) +
  geom_hline(yintercept = as) +
  geom_vline(xintercept = ts) +
  scale_fill_viridis_c(option = "H", alpha = 0.6, na.value = NA, limits = c(min_erro, max_erro)) +
  theme_void() +
  coord_cartesian(expand = FALSE) +
  scale_linewidth_continuous(limits = c(0, 1), range = c(0, 3)) +
  guides(linewidth = "none") +
  geom_segment(
    mapping = aes(
      x = t - 0.5 * !arrow_start,
      xend = t + 0.5 * !arrow_end,
      y = a - 0.5 * !arrow_start,
      yend = a + 0.5 * !arrow_end,
      linewidth = cp_S
    ),
    colour = "gray40",
    alpha = 0.5
  ) +
  geom_point(aes(x = ae_na * t + 0.3, y = ae_na * a + 0.3), shape = 21, fill = NA, stroke = 0.2)






