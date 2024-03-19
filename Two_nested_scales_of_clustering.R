
Two_nested_scales_of_clustering=function(entity,minx,miny,maxx,maxy,window,spat_gra)
{
  library(minpack.lm)
  library(spatstat)
  library(zoo)
  max_dist = max(dist(entity[, 1:2])) * 0.5
  window_ = window
  if (window > max_dist) {
    window_ = max_dist
  }
  point = ppp(entity[, 1], entity[, 2], window = owin(c(minx,maxx), c(miny, maxy)))
  G_r = envelope(point, pcf, nsim = 99, r = seq(0, max_dist,spat_gra), correction = "Ripley")
  
  pcf_ = pcf(point, r = seq(0, max_dist,spat_gra), correction = "Ripley")
  
  r = pcf_$r
  pcf_ = pcf_$iso
  pcf_ = cbind(pcf_, r)
  pcf_ = pcf_[-1, ]
  r = r[-1]
  x = pcf_[, 2]
  y = pcf_[, 1]
  fun = function(x, sigma_s, sigma_l, ps, pl) {
    1 + exp(-0.5 * x^2/(2 * sigma_s^2))/(4 * pi * sigma_s^2 * 
                                           ps) + exp(-0.5 * x^2/(2 * sigma_s^2 + 2 * sigma_l^2))/(2 * 
                                                                                                    pi * (2 * sigma_s^2 + 2 * sigma_l^2) * pl)
  }
  nls.out <- nlsLM(y ~ fun(x, sigma_s, sigma_l, ps, pl), 
                   start = list(sigma_s = 0.1, sigma_l = 1, ps = 0.1, 
                                pl = 0.01))
  result = summary(nls.out)
  sigma_s = result$parameters[1, 1]
  sigma_l = result$parameters[2, 1]
  pl = result$parameters[4, 1]
  ps = result$parameters[3, 1]
  print("OK")
  fit_y = 1 + exp(-0.5 * pcf_[, 2]^2/(2 * sigma_s^2))/(4 * 
                                                         pi * sigma_s^2 * ps) + exp(-0.5 * pcf_[, 2]^2/(2 * 
                                                                                                          sigma_s^2 + 2 * sigma_l^2))/(2 * pi * (2 * sigma_s^2 + 
                                                                                                                                                   2 * sigma_l^2) * pl)
  
  square_R = 1 - sum((fit_y - pcf_[, 1])^2)/sum((pcf_[, 1] - mean(pcf_[, 1]))^2)
  sigma_s = round(sigma_s, 6)
  sigma_l = round(sigma_l, 6)
  pl = round(pl, 6)
  ps = round(ps, 6)
  scale_s = 1.96 * 2 * sigma_s
  scale_l = 1.96 * 2 * sigma_l
  square_R = round(square_R, 3)
  
  plot(G_r, legend = FALSE, xlim = c(0, window_),main = "Two nested scales of clustering")
  abline(h = 1, col = "red", lty = 2,lwd=2)
  lines(r, fit_y, lty = 2, col = "blue",lwd=2)
  text = paste0("Sigma_s:  ", sigma_s, "\n", "Sigma_l:  ", 
                sigma_l, "\n", "ps:  ", ps, "\n", "pl:  ", pl, "\n", 
                "scale_s:  ", scale_s, "\n", "scale_l:  ", scale_l, 
                "\n", "square_R:  ", square_R)
  text(1.1 * mean(r), 0.8 * max(pcf_[, 1]), text, adj = c(0))
}
