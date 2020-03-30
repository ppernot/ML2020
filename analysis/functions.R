dbeCalc = function(atoms) {
  # Estimate DBE from composition
  atoms['C'] - (atoms['H'] + atoms['Cl']) / 2 + atoms['N'] / 2 + 1
}

composTab = function(
  D,
  latoms = NULL,
  dataRepo = '../data/'
) {
  # Build chemical composition tab by reading geometry files
  # for a list of systems. A list of atoms (columns of the table)
  # might be provided to ensure consistency of different tables

  formulae = list()
  lD = length(D)
  for (id in 1:lD) {
    num = sprintf("%4.4i", D[id])
    M = readLines(paste0(dataRepo, 'geometry/frag_', num, '.xyz'))[-c(1, 2)]
    atoms = unlist(lapply(
      M,
      FUN = function(x)
        unlist(strsplit(x, ' '))[1]
    ))
    formulae[[id]] = table(atoms)
  }
  if (is.null(latoms)) {
    latoms = c()
    for (id in 1:lD) {
      latoms = c(latoms, names(formulae[[id]]))
    }
    latoms = unique(latoms)
  }
  compos = matrix(0, ncol = length(latoms), nrow = length(formulae))
  colnames(compos) = latoms
  for (id in 1:lD) {
    for (at in names(formulae[[id]]))
      compos[id, at] = formulae[[id]][at]
  }

  return(list(compos = compos, latoms = latoms))
}


plotBar = function(
  X,
  xlab = 'X',
  ylab = 'Y',
  mean = NULL,
  gPars
) {

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    yaxs = 'i',
    lend = 2
  )

  barplot(
    X,
    col = cols_tr[5],
    border = cols[5],
    ylim = c(0, 1.1 * max(X, na.rm = TRUE)),
    xlab = xlab,
    ylab = ylab
  )
  grid()
  # Redo plot to overlap grid
  barplot(X,
          col = cols_tr[5],
          border = cols[5],
          add = TRUE)

  if (!is.null(mean))
    abline(h = mean, col = cols[2], lty = 2)

  box()

}

plotnHDBE <- function(
  x,
  y,
  z,
  scale = 1,
  xlim = range(x),
  ylim = range(y),
  main = '',
  gPars
) {

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2
  )

  plot(
    x, y,
    type = 'n',
    xlab = 'nH',
    xlim = xlim,
    ylab = 'DBE',
    ylim = ylim,
    main = main
  )
  grid()
  points(x,
         y,
         pch = 1,
         cex = 3 * z / scale,
         col = cols[5])
  box()
}
