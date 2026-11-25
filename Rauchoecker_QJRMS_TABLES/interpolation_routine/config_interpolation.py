
#%% coordinates
x = "xh"
xf = "xf"
zold = "zh"
zfold = "zf"        # only if staggered variables are contained

y = "yh"
yf = "yf"

z_true = "zhval"    # true height of each grid cell (for terrain-following coordinates)

dcoords = [x, y, zold, xf, yf, zfold]

#%% instantaneous output
u = "uinterp"
v = "vinterp"
w = "winterp"
t = "th"
p = "prs"

#%% time + spanwise averaged output
u_mean = "uinterptavg"
v_mean = "vinterptavg"
w_mean = "winterptavg"
t_mean = "thtavg1"
p_mean = "prstavg"
ust = "usttavg"

upup_res = "upup_res"
vpvp_res = "vpvp_res"
wpwp_res = "wpwp_res"
upvp_res = "upvp_res"
upwp_res = "upwp_res"
vpwp_res = "vpwp_res"

upup_sfs = "upup_sgs"
vpvp_sfs = "vpvp_sgs"
wpwp_sfs = "wpwp_sgs"
upvp_sfs = "upvp_sgs"
upwp_sfs = "upwp_sgs"
vpwp_sfs = "vpwp_sgs"

uptp_res = "uptp_res"
vptp_res = "vptp_res"
wptp_res = "wptp_res"
uptp_sfs = "uptp_sgs"
vptp_sfs = "vptp_sgs"
wptp_sfs = "wptp_sgs"
tptp_res = "tptp_res"
tptp_sfs = "tptp_sgs"

eps = "disstentavg"
eps_t = "epsttavg"
rho = "rhotavg"

upupup = "upupup"
upupvp = "upupvp"
upupwp = "upupwp"
vpvpup = "vpvpup"
vpvpvp = "vpvpvp"
vpvpwp = "vpvpwp"
wpwpup = "wpwpup"
wpwpvp = "wpwpvp"
wpwpwp = "wpwpwp"
upuptp = "upuptp"
vpvptp = "vpvptp"
wpwptp = "wpwptp"
uppp = "uppp"
vppp = "vppp"
wppp = "wppp"

#%% varlist to select output variables with original names
varlist_stats = [u_mean, v_mean, w_mean, t_mean, ust, p_mean,
                 upup_res, vpvp_res, wpwp_res, upvp_res, upwp_res, vpwp_res,
                 upup_sfs, vpvp_sfs, wpwp_sfs, upvp_sfs, upwp_sfs, vpwp_sfs,
                 uptp_res, vptp_res, wptp_res, uptp_sfs, vptp_sfs, wptp_sfs,
                 tptp_res, tptp_sfs,
                 rho, eps, eps_t, uppp, vppp, wppp,
                 upupup, upupvp, upupwp, upuptp,
                 vpvpup, vpvpvp, vpvpwp, vpvptp,
                 wpwpup, wpwpvp, wpwpwp, wpwptp,
                 zold,
                ]

#%% Dictionary for time+spanwise averaged statistics
vars_stats = {
    u_mean: "u_mean",
    v_mean: "v_mean",
    w_mean: "w_mean",
    t_mean: "t_mean",
    p_mean: "p_mean",
    upup_res: "upup_res",
    vpvp_res: "vpvp_res",
    wpwp_res: "wpwp_res",
    upvp_res: "upvp_res",
    upwp_res: "upwp_res",
    vpwp_res: "vpwp_res",
    upup_sfs: "upup_sfs",
    vpvp_sfs: "vpvp_sfs",
    wpwp_sfs: "wpwp_sfs",
    upvp_sfs: "upvp_sfs",
    upwp_sfs: "upwp_sfs",
    vpwp_sfs: "vpwp_sfs",
    uptp_res: "uptp_res",
    vptp_res: "vptp_res",
    wptp_res: "wptp_res",
    uptp_sfs: "uptp_sfs",
    vptp_sfs: "vptp_sfs",
    wptp_sfs: "wptp_sfs",
    tptp_res: "tptp_res",
    tptp_sfs: "tptp_sfs",
    ust: "ust",
    eps: "eps",
    eps_t: "eps_t",
    rho: "rho",
    upupup: "upupup",
    upupvp: "upupvp",
    upupwp: "upupwp",
    vpvpup: "vpvpup",
    vpvpvp: "vpvpvp",
    vpvpwp: "vpvpwp",
    wpwpup: "wpwpup",
    wpwpvp: "wpwpvp",
    wpwpwp: "wpwpwp",
    upuptp: "upuptp",
    vpvptp: "vpvptp",
    wpwptp: "wpwptp",
    uppp: "uppp",
    vppp: "vppp",
    wppp: "wppp",
}

#%% Dictionary and variable list for snapshots
vars_snapshots = {
    u: "u",
    v: "v",
    w: "w",
    t: "t",
    p: "p",
    x: "x",
    y: "y",
    }

varlist_snapshots = [u, v, w, t, p]
