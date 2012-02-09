  !ModelConf
  Ns         = 5
  Nimp       = 1
  Nspin      = 1
  ! model      = "hm"
  ! order      = "no"
  ! searchmode = "nofunx"
  ! lat_label  = "n"

  !Hvars    
  d    = 1.d0
  ts   = 0.5d0
  tsp  = 0.d0
  beta = 50.d0
  xmu  = 0.d0
  u    = 2.d0
  tpd  = 0.d0
  tpp  = 0.d0
  ed0  = 0.d0
  ep0  = 0.d0
  v    = 0.d0

  !Loops
  nloop = 100
  ! iantif=0
  chiflag=.false.

  !parameters
  NL     = 1024
  Nw     = 1024
  Ltau   = 256
  Nfit   = 1024
  eps   = 0.01d0
  weigth= 0.8d0
  nread  = 0.d0
  nerr = 1.d-4
  ndelta=0.1d0
  wini = -4.d0
  wfin = 4.d0
  heff = 0.d0
  Nx   = 100
  Ny   = 100
  ! doping2deg=0.0
  cutoff=1.e-8
  eps_error=1.d-5
  nsuccess=2

  !ReadUnits
  Hfile="hamiltonian.restart"
  GMimp1file="impG_1_iw.ed"
  GRimp1file="impG_1_realw.ed"
  GMimp2file="impG_2_iw.ed"
  GRimp2file="impG_2_realw.ed"
  CTimp1file="impChi_1_tau.ed"
  CTimp2file="impChi_2_tau.ed"
  CTimpAfile="impChi_12_tau.ed"
  Ofile="observables.ed"

