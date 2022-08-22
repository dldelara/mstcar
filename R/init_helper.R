getislands = function(mod) {
  check = 1:mod$params$dNd[3]
  network_true = list()
  p = 0
  while(length(check)) {
    network = NULL
    net     = check[1]
    network = c(net, mod$nb$neigh[[net]])
    if (!any(mod$nb$neigh[[net]] == 0)) {
      i = 2
      while(i <= length(network)) {
        network = append(network, mod$nb$neigh[[network[i]]][!(mod$nb$neigh[[network[i]]] %in% network)])
        net = append(net, network[i])
        i = i + 1
      }
    }
    p = p + 1
    network_true[[p]] = sort(net)
    check = check[-which(check %in% net)]
  }
  return(network_true)
}
logit = function(x) return(log(x / (1 - x)))
