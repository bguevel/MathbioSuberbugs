"""Paper-style PDE figures plus a proposed two-population extension.

Based on Simsek et al. (2024), Molecular Systems Biology,
'A rich-get-richer mechanism drives patchy dynamics and resistance evolution
in antibiotic-treated bacteria', Appendix Table S1, page 28.

Run: pip install numpy matplotlib numba
     python bacteria_full_plots.py

Creates a `bacteria_figures/` directory of PNG plots. The original four-field
model is used for paper-inspired baseline figures. The five-field variant is
NOT the published model and has separate, illustrative resistance parameters.
All results are simulated, not data digitized from the article.
"""
from dataclasses import dataclass, replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numba import njit

OUT = Path(__file__).resolve().parent / 'bacteria_figures'
OUT.mkdir(exist_ok=True)

@dataclass(frozen=True)
class Parameters:
    # Published, dimensionless Appendix Table S1
    gamma: float = 1.90
    cm: float = 1.0e6
    kappa_b: float = 0.2
    db: float = 0.012
    xi: float = 0.8
    Ath: float = 5.0
    Dc: float = 2.0e-16
    Da: float = 0.08
    Db: float = 0.0001
    Dn: float = 1.0
    n0: float = 5.0
    # Proposed extension assumptions (not in Appendix Table S1)
    Ath2: float = 20.0
    resistant_growth_factor: float = 0.95
    rho: float = 2.0e-5
    initial_resistant_fraction: float = 0.01

P = Parameters()
N = 41             # Fast run: 41x41, dx=0.18; optional: N=73, DX=0.1
DX = 0.18          # Seed separation remains 0.9 in domain coordinates
DT = 0.006         # Explicit diffusion stability: dt <= dx²/(4*max(D))
T_END = 84.0
SEED = 7
A0_VALUES = np.array([0., 20., 40., 50., 60., 80.])
DOSE_SCAN = np.arange(0., 81., 10.)
PATCH_SURVIVAL_THRESHOLD = 100.0  # Operational definition, not specified by paper


def seed_sites():
    spacing = int(round(0.9 / DX))
    if abs(spacing * DX - 0.9) > 1e-9:
        raise ValueError('DX must evenly divide 0.9')
    starts = (N - 1 - spacing * 5) // 2
    if starts < 0:
        raise ValueError('Grid too small for 6x6 seeding layout')
    return np.array([(starts + i*spacing, starts + j*spacing)
                     for i in range(6) for j in range(6)], dtype=np.int64)

SITES = seed_sites()


def initial_fields(a0, extended=False, seed=SEED, par=P):
    rng = np.random.default_rng(seed)
    cells = np.zeros((N,N), dtype=np.float64)
    samples = np.maximum(rng.normal(30., 12., 36), 0.)
    for (i,j), value in zip(SITES, samples):
        cells[i,j] = value
    f = par.initial_resistant_fraction if extended else 0.
    return (cells*(1-f), cells*f, np.full((N,N), par.n0),
            np.full((N,N), float(a0)), np.zeros((N,N)))


@njit(cache=True)
def evolve(c1, c2, n, a, b, final_time, dt, dx,
           gamma, cm, kb, db, xi, Ath, Dc, Da, Db, Dn,
           Ath2, cost, rho, extended, record_every):
    """JIT explicit 2D Euler with boundary-copy zero-normal-gradient faces.

    Nutrient uptake is (g1*c1 + g2*c2); nutrient recycling is proportional to
    antibiotic-induced lysis. Bacterial diffusivity can be effectively zero.
    """
    q1, q2, qn, qa, qb = (np.empty_like(c1), np.empty_like(c2),
                          np.empty_like(n), np.empty_like(a), np.empty_like(b))
    nx, ny = c1.shape
    steps = round(final_time/dt)
    nt = steps // record_every + 2
    hist = np.zeros((nt, 4)) # time, total c1, total c2, total antibiotic
    h = 0
    dx2 = dx*dx
    for step in range(steps+1):
        if step % record_every == 0 or step == steps:
            s1 = 0.; s2 = 0.; sa = 0.
            for i in range(nx):
                for j in range(ny):
                    s1 += c1[i,j]; s2 += c2[i,j]; sa += a[i,j]
            hist[h,0] = step*dt
            hist[h,1] = s1
            hist[h,2] = s2
            hist[h,3] = sa/(nx*ny)
            h += 1
        if step == steps:
            break
        for i in range(nx):
            im = i-1 if i>0 else 0
            ip = i+1 if i<nx-1 else nx-1
            for j in range(ny):
                jm = j-1 if j>0 else 0
                jp = j+1 if j<ny-1 else ny-1
                x1, x2 = c1[i,j], c2[i,j]
                x = x1+x2
                g1 = n[i,j]/(1.+n[i,j]) * max(0., 1.-x/cm)
                g2 = cost*g1 if extended else 0.
                ell1 = gamma*g1 if a[i,j]>=Ath else 0.
                ell2 = gamma*g2 if (extended and a[i,j]>=Ath2) else 0.
                conv = rho*x1*x2 if extended else 0.
                l1 = (c1[im,j]+c1[ip,j]+c1[i,jm]+c1[i,jp]-4*x1)/dx2
                l2 = (c2[im,j]+c2[ip,j]+c2[i,jm]+c2[i,jp]-4*x2)/dx2
                ln = (n[im,j]+n[ip,j]+n[i,jm]+n[i,jp]-4*n[i,j])/dx2
                la = (a[im,j]+a[ip,j]+a[i,jm]+a[i,jp]-4*a[i,j])/dx2
                lb = (b[im,j]+b[ip,j]+b[i,jm]+b[i,jp]-4*b[i,j])/dx2
                q1[i,j] = max(0., x1+dt*((g1-ell1)*x1-conv+Dc*l1))
                q2[i,j] = max(0., x2+dt*((g2-ell2)*x2+conv+Dc*l2))
                qn[i,j] = max(0., n[i,j]+dt*((xi*ell1-g1)*x1+(xi*ell2-g2)*x2+Dn*ln))
                qa[i,j] = max(0., a[i,j]+dt*(-kb*b[i,j]*a[i,j]+Da*la))
                qb[i,j] = max(0., b[i,j]+dt*(ell1*x1+ell2*x2-db*b[i,j]+Db*lb))
        c1, q1 = q1, c1
        c2, q2 = q2, c2
        n, qn = qn, n
        a, qa = qa, a
        b, qb = qb, b
    return c1, c2, n, a, b, hist[:h]


def solve(a0, extended=False, par=P, seed=SEED, record_every=250):
    c1, c2, n, a, b = initial_fields(a0, extended, seed, par)
    return evolve(c1,c2,n,a,b,T_END,DT,DX,par.gamma,par.cm,
                  par.kappa_b,par.db,par.xi,par.Ath,par.Dc,par.Da,
                  par.Db,par.Dn,par.Ath2,par.resistant_growth_factor,
                  par.rho,extended,record_every)


def patch_values(field):
    return np.array([field[i,j] for i,j in SITES])


def gini(values):
    s = np.sort(np.maximum(np.asarray(values,float),0.))
    if s.sum() <= 1e-12:
        return np.nan # Undefined for total extinction
    m = len(s)
    return (2*np.dot(np.arange(1,m+1),s)/(m*s.sum())-(m+1)/m)


def metrics(result):
    c1,c2 = result[0],result[1]
    vals = patch_values(c1+c2)
    return (np.mean(vals > PATCH_SURVIVAL_THRESHOLD)*100.,
            gini(vals), np.max(vals), np.sum(vals),
            np.sum(c2)/(np.sum(c1)+np.sum(c2)+1e-30))


def savefig(fig, filename):
    fig.savefig(OUT/filename, dpi=155, bbox_inches='tight')
    plt.close(fig)
    print('Wrote',filename, flush=True)


def baseline_figures(cache):
    # Figure 1C-style density maps, comparable conditions and one realization
    fig,axes=plt.subplots(1,len(A0_VALUES),figsize=(16,3.4),layout='constrained')
    vmax=max(np.max(cache[float(d)][0]) for d in A0_VALUES)
    for ax,d in zip(axes,A0_VALUES):
        im=ax.imshow(cache[float(d)][0],vmin=0,vmax=vmax,cmap='magma',origin='lower')
        ax.set_title(f'$a_0={d:g}$')
        ax.set_xticks([]);ax.set_yticks([])
    fig.colorbar(im,ax=axes,shrink=.75,label='Bacterial concentration')
    fig.suptitle('Original four-field PDE: final bacterial distributions, t = 84')
    savefig(fig,'01_paper_spatial_dose.png')

    # Figure EV2 + Figure 1D-style quantification
    data=np.array([metrics(cache[float(d)]) for d in DOSE_SCAN])
    fig,axs=plt.subplots(1,3,figsize=(13,3.5),layout='constrained')
    for ax,y,ylabel in zip(axs,[data[:,0],data[:,1],data[:,2]],
                           ['Surviving seed sites (%)','Gini of 36 seed-site densities',
                            'Largest final seed-site density']):
        ax.plot(DOSE_SCAN,y,'o-',lw=1.8)
        ax.set(xlabel='Initial antibiotic $a_0$',ylabel=ylabel)
        ax.grid(alpha=.22)
    axs[0].set_ylim(-5,105);axs[1].set_ylim(-.02,1.02)
    fig.suptitle('Original model: patch survival, inequality, and maximum growth')
    savefig(fig,'02_paper_survival_gini_max.png')

    # Figure EV3-style initial versus final density scatter
    fig,axs=plt.subplots(1,3,figsize=(12,3.6),layout='constrained')
    initial=patch_values(initial_fields(0.)[0])
    for ax,d in zip(axs,[0.,40.,60.]):
        vals=patch_values(cache[d][0])
        ax.scatter(initial, vals,s=25,alpha=.8)
        ax.set(xlabel='Initial seed-site density',ylabel='Final seed-site density',
               title=f'$a_0={d:g}$')
        ax.grid(alpha=.2)
    fig.suptitle('Original model: initial vs final bacterial density (36 sites)')
    savefig(fig,'03_paper_initial_vs_final.png')

    # Time series, useful check on transients
    fig,ax=plt.subplots(figsize=(7,4))
    for d in [0.,40.,60.,80.]:
        h=cache[d][5]
        ax.plot(h[:,0],h[:,1],label=f'$a_0={d:g}$')
    ax.set(xlabel='Time',ylabel='Total bacteria (sum over grid)',
           title='Original model: total bacterial population trajectories')
    ax.legend();ax.grid(alpha=.2)
    savefig(fig,'04_paper_time_series.png')


def sensitivity_figure():
    # Figure EV4-inspired parameter analysis; Gini at each antibiotic dose.
    variations = [('gamma',np.array([1.3,1.9,2.5])),
                  ('kappa_b',np.array([0.10,0.20,0.30])),
                  ('Ath',np.array([2.,5.,8.]))]
    doses=np.array([0.,20.,40.,60.,80.])
    fig,axs=plt.subplots(1,3,figsize=(13.7,4),layout='constrained')
    for ax,(key,vals) in zip(axs,variations):
        z=np.zeros((len(vals),len(doses)))
        for k,value in enumerate(vals):
            pp=replace(P,**{key:float(value)})
            for j,dose in enumerate(doses):
                g=metrics(solve(dose,par=pp))[1]
                z[k,j]=0.0 if np.isnan(g) else g
        im=ax.imshow(z,origin='lower',aspect='auto',vmin=0,vmax=1,cmap='viridis')
        ax.set_xticks(range(len(doses)),[f'{v:g}' for v in doses])
        ax.set_yticks(range(len(vals)),[f'{v:g}' for v in vals])
        ax.set(xlabel='Initial antibiotic $a_0$',ylabel=key,
               title=f'Effect of ${key}$ on final Gini')
    fig.colorbar(im,ax=axs,label='Gini (0 assigned for extinction)')
    fig.suptitle('Original model: parameter sensitivity (Figure EV4-inspired)')
    savefig(fig,'05_paper_parameter_sensitivity.png')


def diffusion_figure():
    # Appendix Figure S7: bacterial dispersal and patchiness.
    dcs=[2e-16,1e-3,1e-2,1.0]
    doses=[0.,40.,60.,80.]
    fig,axs=plt.subplots(len(dcs),len(doses),figsize=(11,11),layout='constrained')
    vmax=250.0
    for i,dc in enumerate(dcs):
        pp=replace(P,Dc=dc)
        for j,d in enumerate(doses):
            # Explicit Euler stable only if dt<=dx²/(4*Dc).
            # Dc=1 is stable here with dt=.006 and dx=.18.
            z=solve(d,par=pp)[0]
            im=axs[i,j].imshow(z,origin='lower',cmap='magma',vmin=0,vmax=vmax)
            axs[i,j].set_xticks([]);axs[i,j].set_yticks([])
            if i==0: axs[i,j].set_title(f'$a_0={d:g}$')
            if j==0: axs[i,j].set_ylabel(f'$D_c={dc:g}$')
    fig.colorbar(im,ax=axs,shrink=.6,label=f'Density (color scale capped at {vmax:g})')
    fig.suptitle('Original model: effect of bacterial diffusion (Appendix Figure S7-inspired)')
    savefig(fig,'06_paper_diffusion.png')


def extension_figures():
    results={float(d):solve(float(d),extended=True) for d in DOSE_SCAN}
    fig,axs=plt.subplots(2,len(A0_VALUES),figsize=(16,6),layout='constrained')
    top=max(np.max(results[float(d)][0]) for d in A0_VALUES)
    bot=max(np.max(results[float(d)][1]) for d in A0_VALUES)
    for j,d in enumerate(A0_VALUES):
        q=results[float(d)]
        im1=axs[0,j].imshow(q[0],origin='lower',cmap='viridis',vmin=0,vmax=top)
        im2=axs[1,j].imshow(q[1],origin='lower',cmap='magma',vmin=0,vmax=bot)
        axs[0,j].set_title(f'$a_0={d:g}$')
        for ax in axs[:,j]:ax.set_xticks([]);ax.set_yticks([])
    axs[0,0].set_ylabel('Sensitive $c_1$');axs[1,0].set_ylabel('Resistant $c_2$')
    fig.colorbar(im1,ax=axs[0,:],shrink=.6,label='Sensitive density')
    fig.colorbar(im2,ax=axs[1,:],shrink=.6,label='Resistant density')
    fig.suptitle('Five-field extension: sensitive and resistant bacteria, t = 84')
    savefig(fig,'07_extension_spatial_sensitive_resistant.png')

    frac=[];sensitive=[];resistant=[];survival=[]
    for dose in DOSE_SCAN:
        q=results[float(dose)]
        frac.append(metrics(q)[4])
        sensitive.append(q[0].sum());resistant.append(q[1].sum())
        survival.append(metrics(q)[0])
    fig,axs=plt.subplots(1,3,figsize=(13.8,3.7),layout='constrained')
    axs[0].plot(DOSE_SCAN,np.array(frac)*100,'o-')
    axs[0].set(ylabel='Resistant share of total bacteria (%)',ylim=(-3,103))
    axs[1].plot(DOSE_SCAN,sensitive,'o-',label='Sensitive')
    axs[1].plot(DOSE_SCAN,resistant,'s-',label='Resistant')
    axs[1].set(ylabel='Total final bacterial density (grid sum)')
    axs[1].legend()
    axs[2].plot(DOSE_SCAN,survival,'o-')
    axs[2].set(ylabel='Surviving sites (%)',ylim=(-3,103))
    for ax in axs:
        ax.set_xlabel('Initial antibiotic $a_0$');ax.grid(alpha=.2)
    fig.suptitle('Proposed resistance extension: dose dependence')
    savefig(fig,'08_extension_resistant_fraction_vs_dose.png')

    fig,axs=plt.subplots(1,2,figsize=(11,3.8),layout='constrained')
    for d in [0.,40.,60.,80.]:
        h=results[d][5]
        axs[0].plot(h[:,0],h[:,1],label=f'$a_0={d:g}$')
        axs[1].plot(h[:,0],h[:,2]/(h[:,1]+h[:,2]+1e-30)*100,
                    label=f'$a_0={d:g}$')
    axs[0].set(xlabel='Time',ylabel='Sensitive count (grid sum)')
    axs[1].set(xlabel='Time',ylabel='Resistant fraction (%)',ylim=(-3,103))
    for ax in axs:ax.legend();ax.grid(alpha=.2)
    fig.suptitle('Proposed resistance extension: transient dynamics')
    savefig(fig,'09_extension_resistance_time_series.png')

    # Sensitivity to rho: rho*c1*c2 is an interaction, NOT spontaneous mutation.
    rhos=np.array([0.,5e-6,2e-5,1e-4])
    doses=np.array([0.,20.,40.,60.,80.])
    z=np.zeros((len(rhos),len(doses)))
    for i,rho in enumerate(rhos):
        for j,d in enumerate(doses):
            z[i,j]=100*metrics(solve(d,extended=True,par=replace(P,rho=rho)))[4]
    fig,ax=plt.subplots(figsize=(7.2,4.5),layout='constrained')
    im=ax.imshow(z,origin='lower',aspect='auto',cmap='viridis',vmin=0,vmax=100)
    ax.set_xticks(range(len(doses)),[f'{d:g}' for d in doses])
    ax.set_yticks(range(len(rhos)),[f'{r:.0e}' for r in rhos])
    ax.set(xlabel='Initial antibiotic $a_0$',ylabel='Interaction coefficient $\\rho$',
           title='Extension sensitivity: resistant fraction (%)')
    fig.colorbar(im,ax=ax,label='Final resistant fraction (%)')
    savefig(fig,'10_extension_rho_sensitivity.png')


def main():
    import argparse
    cli = argparse.ArgumentParser()
    cli.add_argument('--part', choices=['all','baseline','sensitivity','diffusion','extension'], default='all')
    part = cli.parse_args().part
    if DT > DX**2/(4*max(P.Dn,P.Da,P.Db,1.0)):
        raise ValueError('Explicit diffusion time step violates 2D stability limit')
    if part in ('all','baseline'):
        print('Computing baseline antibiotic sweep...',flush=True)
        cached={float(d):solve(float(d)) for d in sorted(set(DOSE_SCAN)|set(A0_VALUES))}
        baseline_figures(cached)
    if part in ('all','sensitivity'):
        print('Computing parameter sensitivity...',flush=True)
        sensitivity_figure()
    if part in ('all','diffusion'):
        print('Computing bacterial diffusion plots...',flush=True)
        diffusion_figure()
    if part in ('all','extension'):
        print('Computing resistant-bacteria extension...',flush=True)
        extension_figures()
    print(f'Finished. All plots are in: {OUT}')


if __name__ == '__main__':
    main()
