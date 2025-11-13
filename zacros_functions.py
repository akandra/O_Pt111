#!/usr/bin/env python
# coding: utf-8

# # Zacros related functions

# In[ ]:


import numpy as np
import subprocess, os
import matplotlib.pyplot as plt
import ase, ase.io
from ase.visualize import view
from pathlib import Path
import pandas as pd
from IPython.display import clear_output, display
import time
from scipy.spatial import cKDTree,ConvexHull

# Constants
en_file_suffix  = 'energy.dat'
rdf_file_suffix = 'rdf.dat'
csd_file_suffix = 'csd.dat'
ccd_file_suffix = 'ccd.dat'
acc_file_suffix = 'acc.dat'

# RDF

def rdf(confs, coverage, v1, v2, r_max=None, dr=0.1, g_ref=None, lattice_constant=None, ax=None, plot=False):
   """
   Compute 2D radial distribution function g(r) with periodic boundary conditions.

   Parameters
   ----------
   confs : list of (N,2) array_like
         Cartesian coordinates of points.
   v1, v2 : array_like, shape (2,)
         Box vectors spanning the periodic cell.
   r_max : float or None
         Maximum radius to compute g(r). If None use min(box lengths, box diagonal)/2.
   dr : float
         Bin width.

   Returns
   -------
   dist : (M,) array
         Bin midpoints.
   g : (M,) array
         Radial distribution function values.
   """
      
   # choose r_max if not provided
   if r_max is None:
         l1 = np.linalg.norm(v1)
         l2 = np.linalg.norm(v2)
         l3 = np.linalg.norm(v1 + v2)
         r_max = min(l1, l2, l3) / 2.0

   cell_matrix = np.vstack([v1, v2]).T
   cell_matrix_inv = np.linalg.inv(cell_matrix)

   nbins = int(np.ceil(r_max / dr))
   edges = np.linspace(0.0, r_max, nbins + 1)
   g = np.zeros(nbins)
   for conf in confs:  
      counts = np.zeros(nbins, dtype=int)
      for i1, r1 in enumerate(conf[:-1]):
         for r2 in conf[i1+1:]:
            r = r2 - r1
            # apply PBC
            frac_r = cell_matrix_inv @ r
            r = cell_matrix @ (frac_r - np.rint(frac_r))
            d = np.linalg.norm(r)
            if (d>0) and (d<=r_max): counts[int(d/dr)] += 1
      # Get 2D lattice ideal gas normalization factor 
      if g_ref is not None: 
         counts_n = np.zeros_like(counts, dtype=float)
         np.divide(counts, g_ref, out=counts_n, where=g_ref!=0)
         g = g + counts_n / len(conf) / coverage  # normalize by coverage
      else:
         g = g + counts / len(conf) / coverage  # normalize by coverage

   dist = 0.5 * (edges[:-1] + edges[1:])

   # Normalization by number of configurations 
   # factor of 2 for unordered pairs
   g = 2*g / len(confs)

   # plotting
   if plot:
      if ax is None:
            fig, ax = plt.subplots()
      if lattice_constant is None:
         ax.set_xlabel(r'$r (\mathrm{\AA})$')
         lattice_constant = 1
      else:
         ax.set_xlabel(r'$r / a_0$')
      ax.set_ylabel(r'g(r)')
      ax.plot(dist/lattice_constant, g, marker='o',color='k')
      #ax.set_ylim(bottom=0)
      ax.grid(True)
      if plot:
            plt.show()
      
   return dist/lattice_constant, g 

# Cluster Size and Circularity Functions

class UnionFind_a:
    def __init__(self, size):
      
        # Initialize the parent array with each element as its own representative
        self.parent = list(range(size))
    
    def find(self, i):
      
        # If i itself is root or representative
        if self.parent[i] == i:
            return i
          
        # Else recursively find the representative of the parent
        return self.find(self.parent[i])
    
    def union(self, i, j):
      
        # Representative of set containing i
        irep = self.find(i)
        
        # Representative of set containing j
        jrep = self.find(j)
        
        # Make the representative of i's set be the representative of j's set
        self.parent[irep] = jrep


class UnionFind:
    def __init__(self, n):
        self.parent = np.arange(n)
        self.rank = np.zeros(n, dtype=int)

    def find(self, a):
        p = self.parent
        while p[a] != a:
            p[a] = p[p[a]]
            a = p[a]
        return a

    def union(self, a, b):
        ra = self.find(a); rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        else:
            self.parent[rb] = ra
            if self.rank[ra] == self.rank[rb]:
                self.rank[ra] += 1

def cluster_points_periodic(points, v1, v2, cutoff):
    """
    Cluster 2D points with periodic boundary conditions.

    Parameters
    ----------
    points : (N,2) array_like
        Cartesian coordinates of points (float).
    v1, v2 : The box is spanned by v1 and v2.
    cutoff : float
        Distance cutoff for connectivity (e.g. 3rd NN distance).

    Returns
    -------
    labels : (N,) int
        Cluster label for each original point (0..nclusters-1).
    clusters : list of ndarray
        Indices of points in each cluster.
    sizes : list of int 
        Sizes of clusters (only if return_sizes True).
    """
    pts = np.asarray(points, dtype=float)
    if pts.size == 0:
      return np.array([], dtype=int), [], []

    N = len(pts)
    # Build augmented points = original points shifted by translations i*v1 + j*v2 with i,j in {-1,0,1}
    shifts = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
    aug_pts = np.zeros((N * len(shifts), 2), dtype=float)
    orig_idx = np.zeros(N * len(shifts), dtype=int)

    k = 0
    for si, (i, j) in enumerate(shifts):
        shift_vec = i * v1 + j * v2
        aug_pts[k:k+N] = pts + shift_vec
        orig_idx[k:k+N] = np.arange(N)
        k += N

    # KD-tree on augmented points
    tree = cKDTree(aug_pts)
    pairs = tree.query_pairs(cutoff, output_type='ndarray')  # array of shape (M,2)

    uf = UnionFind(N)
    for a, b in pairs:
        ia = orig_idx[a]
        ib = orig_idx[b]
        if ia != ib:
            uf.union(ia, ib)

    # Extract roots and relabel to contiguous labels
    roots = np.array([uf.find(i) for i in range(N)])
    unique_roots, inv = np.unique(roots, return_inverse=True)
    labels = inv
    clusters = [np.nonzero(labels == k)[0] for k in range(len(unique_roots))]
    sizes = [len(c) for c in clusters]

    return labels, clusters, sizes

def cluster_circularity_periodic(points, v1, v2, cutoff):
    """
    Cluster points with PBC and compute cluster area, perimeter and circularity.
    Returns a dict with:
      labels, clusters, sizes, areas, perimeters, circularities
    Notes:
      - Area/perimeter computed from convex hull of each cluster after unwrapping PBC.
      - Circularity = 4*pi*area / perimeter^2 (0 if area==0 or perimeter==0).
    """
    pts = np.asarray(points, dtype=float)
    labels, clusters, sizes = cluster_points_periodic(pts, v1, v2, cutoff)
    if len(clusters) == 0:
        return dict(labels=labels, clusters=clusters, sizes=sizes,
                    areas=[], perimeters=[], circularities=[])
    # Build matrix V so that p = V @ frac  (frac in R^2)
    V = np.vstack([np.asarray(v1), np.asarray(v2)]).T  # shape (2,2)
    invV = np.linalg.inv(V)
    frac = (invV @ pts.T).T   # fractional coordinates (may be outside [0,1))
    frac_mod = frac - np.floor(frac)  # map to [0,1)
    areas = []
    perimeters = []
    circularities = []
    for idxs in clusters:
        n = len(idxs)
        if n == 0:
            areas.append(0.0); perimeters.append(0.0); circularities.append(0.0); continue
        fi = frac_mod[idxs]  # shape (n,2)
        # compute cluster mean and unwrap each point by nearest integer shift towards mean
        mean_f = fi.mean(axis=0)
        # compute integer shifts that bring each point close to mean
        shifts = np.rint(mean_f - fi).astype(int)  # shape (n,2)
        fi_unwrapped = fi + shifts
        coords = (V @ fi_unwrapped.T).T  # unwrapped Cartesian coords
        if n == 1:
            areas.append(0.0); perimeters.append(0.0); circularities.append(0.0); continue
        if n == 2:
            d = np.linalg.norm(coords[0] - coords[1])
            area = 0.0
            perim = 2.0 * d
            circ = 0.0
            areas.append(area); perimeters.append(perim); circularities.append(circ); continue
        # n >= 3: convex hull
        try:
            hull = ConvexHull(coords)
            hull_pts = coords[hull.vertices]
            # perimeter
            diffs = np.diff(np.vstack([hull_pts, hull_pts[0]]), axis=0)
            perim = np.sqrt((diffs**2).sum(axis=1)).sum()
            area = hull.volume  # area for 2D
        except Exception:
            # fallback: compute convex hull by monotone chain if ConvexHull fails
            # simple fallback: approximate area/perimeter from bounding box
            xmin, ymin = coords.min(axis=0)
            xmax, ymax = coords.max(axis=0)
            area = (xmax - xmin) * (ymax - ymin)
            perim = 2.0 * ((xmax - xmin) + (ymax - ymin))
        circ = 4.0 * np.pi * area / (perim**2) if perim > 0 and area > 0 else 0.0
        areas.append(float(area)); perimeters.append(float(perim)); circularities.append(float(circ))
    return dict(labels=labels, clusters=clusters, sizes=sizes,
                areas=areas, perimeters=perimeters, circularities=circularities)

def cluster_circularity_periodic_avg(confs, v1, v2, cutoff, nbins=20):
    """
    Compute cluster area, perimeter and circularity averaged over confs.

    Returns a dict with:
      labels, clusters, sizes, areas, perimeters, circularities
      averaged over all configurations in confs.
    Notes:
      - Area/perimeter computed from convex hull of each cluster after unwrapping PBC.
      - Circularity = 4*pi*area / perimeter^2 (0 if area==0 or perimeter==0).
    """

    sizes = []
    circularities = []
    for points in confs:
      res = cluster_circularity_periodic(points, v1, v2, cutoff=cutoff)
      sizes.extend(res['sizes'])
      circularities.extend(res['circularities'])

    sizes = np.asarray(sizes, dtype=int)
    sizes_hist = np.zeros(sizes.max())
    for s in sizes:
        sizes_hist[s-1] += 1
    sizes_hist /= len(confs)

    edges = np.linspace(0, 1, nbins+1)
    circ_hist = np.zeros(len(edges)-1)
    bin_index = np.digitize(circularities, edges, right=False) - 1
    for bin_idx in bin_index: circ_hist[bin_idx] += 1
    circ_hist /= len(confs)
    circ_bin_centers = 0.5 * (edges[:-1] + edges[1:])

    return np.arange(1,sizes.max()), sizes_hist, circ_bin_centers, circ_hist

def plot_cluster_size_and_circularity(sizes, circularities, bins=None):
    """
    Quick plot helper: histogram of sizes and circularity distribution.
    """
    import matplotlib.pyplot as plt
    sizes = np.asarray(sizes)
    circularities = np.asarray(circularities)
    if bins is None:
        bins = np.arange(1, sizes.max()+2) if sizes.size>0 else 10
    fig, ax = plt.subplots(1,2,figsize=(10,4))
    if sizes.size>0:
        ax[0].hist(sizes, bins=bins, color='C0', log=True)
        ax[0].set_xlabel('Cluster size (sites)'); ax[0].set_ylabel('Number of clusters (log)')
    else:
        ax[0].text(0.5,0.5,'no clusters',ha='center',va='center')
    if circularities.size>0:
        ax[1].hist(circularities[circularities>0], bins=30, color='C1')
        ax[1].set_xlabel('Circularity'); ax[1].set_ylabel('Counts')
    else:
        ax[1].text(0.5,0.5,'no circularities',ha='center',va='center')
    plt.tight_layout()
    plt.show()

def plot_avgs(df, temperature=300, cov=0.01, ints='',fraction_eq=0.5,figsize=(8,6)):
   """
   Plots rdf, energy vs. time, csd and ccd
   """
   def get_xy(fname):
      try:
         return np.loadtxt(fname, unpack=True)
      except :
         return None

   # Get df rows satisfying the conditions
   conditions =   (df['int'] == ints) & \
                  (df['temperature'] == temperature) & \
                  (df['cov'] == cov)
   
   if len(df[conditions])==0: 
      return print('No data exist for the conditions specified')

   # Get run data dir name
   rundir = df[conditions].dir.iloc[0] 

   fig, axes = plt.subplots(3,2,figsize=figsize)
   fig.suptitle(fr'$T={temperature}$K, $\theta={cov}$ML, {df[conditions].int.iloc[0]}')

   # RDF
   fname = rundir.parent.parent / 'results' / f'{rundir.name}_{rdf_file_suffix}'
   x,y = get_xy(fname)

   ax = axes[1,0]
   if x is None: 
      ax.text(0.5, 0.5, f"{fname.name} not found", ha='center', transform=ax.transAxes)
   else:
      ax.plot(x, y, marker='o', color='k', linestyle='-', markersize=4)
      ax.grid()
   ax.set_xlabel(r'$r/a_0$')
   ax.set_ylabel('Reduced $g(r)$')

   # CSD
   fname = rundir.parent.parent / 'results' / f'{rundir.name}_{csd_file_suffix}'
   x,y,z = get_xy(fname)

   ax = axes[0,1]
   if x is None: 
      ax.text(0.5, 0.5, f"{fname.name} not found", ha='center', transform=ax.transAxes)
   else:
      ax.plot(x, y, marker='o', color='k', linestyle='-', markersize=4)
      ax.grid()
   ax.set_xlabel(r'Cluster size')
   ax.set_ylabel('CSD average counts')
   
   # ACC
   ax = axes[1,1]
   if x is None: 
      ax.text(0.5, 0.5, f"{fname.name} not found", ha='center', transform=ax.transAxes)
   else:
      ax.plot(x, z, marker='o', color='k', linestyle='-', markersize=4)
      ax.grid()
   ax.set_xlabel(r'Cluster size')
   ax.set_ylabel('ACC average counts')

   # Energy vs. time
   fname = rundir.parent.parent / 'results' / f'{rundir.name}_{en_file_suffix}'
   x,y = get_xy(fname)

   ax = axes[0,0]
   if x is None: 
      ax.text(0.5, 0.5, f"{fname.name} not found", ha='center', transform=ax.transAxes)
   else:
      ax.plot(x, y, color='k', linestyle='-')
      ax.axvline(x=fraction_eq*x[-1], color='r', linestyle='--', linewidth=2)
      ax.grid()
   ax.set_xlabel('Time (s)')
   ax.set_ylabel('Energy (eV)')

   # ACC histogram
   fname = rundir.parent.parent / 'results' / f'{rundir.name}_{acc_file_suffix}'
   x,y = get_xy(fname)

   ax = axes[2,0]
   if x is None: 
      ax.text(0.5, 0.5, f"{fname.name} not found", ha='center', transform=ax.transAxes)
   else:
      ax.plot(x, y, marker='o', color='k', linestyle='-', markersize=4)
      ax.grid()
   ax.set_xlabel(r'Accessibility')
   ax.set_ylabel('Average counts')

   fig.tight_layout()
   plt.show()
   
   return rundir

# ## Lattice plotting function

# In[ ]:

def plot_cluster_size_distribution(sizes):
  """
  Plots the cluster size distribution given a list of cluster sizes.
  
  Parameters:
  ----------
  sizes : list or array-like
      List of cluster sizes (number of sites in each cluster).
  """   

  sizes = np.array(sizes)
  if sizes.size == 0:
    print("No clusters found.")
  else:
    counts = np.bincount(sizes)          # index = cluster size
    counts = counts[1:]               # ignore size=0 clusters
    s = np.arange(1, len(counts)+1)

    pmf = counts / counts.sum()
    cdf = np.cumsum(pmf)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # histogram (log-scale y helps if broad distribution)
    axes[0].bar(s, counts, color='C0')
    #axes[0].set_yscale('log')
    axes[0].set_xlabel('Cluster size (sites)')
    axes[0].set_ylabel('Number of clusters')
    axes[0].set_title('Cluster size distribution')
    axes[0].grid(True)

    # PMF + CDF
    axes[1].plot(s, pmf, marker='o', label='PMF')
    axes[1].set_xlabel('Cluster size (sites)')
    axes[1].set_ylabel('P(size)')
    ax2 = axes[1].twinx()
    ax2.plot(s, cdf, color='0.4', linestyle='--', label='CDF')
    ax2.set_ylabel('Cumulative P')
    axes[1].grid(True)
    axes[1].legend(loc='upper left')
    ax2.legend(loc='upper right')

    plt.tight_layout()
    print(f"n_clusters: {len(sizes)}, mean size: {sizes.mean():.2f}, median: {np.median(sizes):.2f}")
    plt.show()

def lattice_plot(lattice_input_file, reps=None, idx = None, show_axes = False, pause=-1, show=True, close=False, show_sites_ids=False, file_name=None, figsize=(8,6)):
    """
    Visualizes a lattice defined in a Zacros lattice input file.
    Parameters
    ----------
    lattice_input_file : str or Path
        Path to the Zacros lattice input file.
    reps : list or tuple, optional
        Number of repetitions of the unit cell in x and y directions. If None, it uses the values from the lattice input file.
    idx : int, optional
        If provided, visualizes the configuration at the specified index from history_output.txt.
    show_axes : bool, optional
        If True, shows the axes of the plot. Default is False.
    pause : float, optional
        Time in seconds to pause after showing the figure. Default is -1, which means it will wait indefinitely.
    show : bool, optional
        If True, shows the figure on the screen. Default is True.
    close : bool, optional
        If True, closes the figure window after pause time. Default is False.
    show_sites_ids : bool, optional
        If True, shows the binding sites id on the figure. Default is False.
    file_name : str, optional
        If provided, saves the figure to the specified file name. Default is None.
    Returns
    -------
    None
    """

    # Read lattice input file

    try:
        with open(lattice_input_file, 'r') as f:
            content = [line for line in f.readlines() if (line.strip() and not line.startswith('#'))]
    except FileNotFoundError:
        raise FileNotFoundError(f"Lattice input file '{lattice_input_file}' not found.")
    finally:
        content = [line.split('#')[0] for line in content]

    wind_rose = {   "self": (0,0),
                    "north": (0,1),
                    "northeast": (1,1),
                    "east": (1,0),
                    "southeast": (1,-1) }

    for i,line in enumerate(content):
            if 'cell_vectors' in line:
                unit_cell = np.array([ [float(x) for x in content[i+1].split()],
                              [float(x) for x in content[i+2].split()] ])
            if 'repeat_cell' in line:
                repeat_cell = np.array([ int(x) for x in line.split()[1:3] ])
            if 'n_site_types' in line:
                n_site_types = int(line.split()[1])
            if 'n_cell_sites' in line:
                n_cell_sites = int(line.split()[1])
            if 'site_types' in line:
                site_types_names = line.split()[1:]
            if 'site_coordinates' in line:
                site_coordinates = []
                for j in range(n_cell_sites):
                    site_coordinates.append([float(x) for x in content[i+1+j].split()[:2]])
                site_uc_coordinates = np.array(site_coordinates) @ unit_cell


    if reps is not None:
        repeat_cell = np.array(reps)

    site_coordinates = []
    site_types = []
    for i in range(repeat_cell[0]):
        for j in range(repeat_cell[1]):
            shift = np.array([i,j]) @ unit_cell
            for coords,st in zip(site_uc_coordinates , site_types_names):
                site_types.append(st)
                site_coordinates.append(coords + shift)
                site_types.append(st)

    # Read configurations from history_output.txt file

    try:
        with open(lattice_input_file.parent / 'history_output.txt', 'r') as f:
            content = f.readlines()
    except FileNotFoundError:
        print(f" File '{lattice_input_file.parent / 'history_output.txt'}' not found.")
    finally:

        # Get the species list
        species = content[1].split()[1:]

        # Get number of configurations
        for line in content:
            if 'configuration' in line:
                n_confs = int(line.split()[1])

        n_sites = repeat_cell[0] * repeat_cell[1] * n_cell_sites
        confs = np.zeros((n_confs, n_sites), dtype=int)
        for i_conf in range(n_confs):
            for i in range(n_sites):
                confs[i_conf, i] = int(content[7 + i_conf*(n_sites+1) +i].split()[2])

    v1 = repeat_cell[0] * unit_cell[0]
    v2 = repeat_cell[1] * unit_cell[1]
    xvalues = [0.0, v1[0], v1[0] + v2[0], v2[0], 0.0]
    yvalues = [0.0, v1[1], v1[1] + v2[1], v2[1], 0.0]

    markers =   ["v", "s", "o", "D", "p", "^", "+", "x", "*", "P", "H", "X", "d", "h", ",", ".", "<", ">", "1", "2"]
    colors = ["lightgray", "r", "g", "b", "m", "c", "k",
        "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
        "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
        "gold", "turquoise", "lime", "indigo"]


    fig, ax = plt.subplots(figsize=figsize)

    x_len = abs(v2[0] - v1[0])
    ax.set_xlim([0.0 - 0.1 * x_len, v1[0] + v2[0] + 0.1 * x_len])
    y_len = abs(v2[1] - v1[1])
    ax.set_ylim([0.0 - 0.1 * y_len, v1[1] + v2[1] + 0.1 * y_len])
    ax.set_aspect(1.0)

    ax.set_xlabel(r'x ($\AA$)')
    ax.set_ylabel(r'y ($\AA$)')

    ax.plot(xvalues, yvalues, color='k', linestyle="dashed", linewidth=1, zorder=1)

    if idx is not None:
        confs = [confs[idx]]
    for ic,conf in enumerate(confs):

        ax.cla()

        for i, st_i in enumerate(sorted(list(set(site_types)))):
            for occ in range(len(species)+1):
                xvalues = [x for (x, y), st, spec in zip(site_coordinates, site_types, conf) if (st == st_i and spec == occ)]
                yvalues = [y for (x, y), st, spec in zip(site_coordinates, site_types, conf) if (st == st_i and spec == occ)]

                ax.scatter(
                    xvalues,
                    yvalues,
                    color=colors[occ],
                    marker=markers[i],
                    s=1.5*440 / np.sqrt(len(site_coordinates)),
                    zorder=2,
                    label=st_i if occ == 0 else species[occ - 1]
                )
                if show_sites_ids:
                    for i, (x, y) in enumerate(site_coordinates):
                        plt.annotate(str(i), (x, y), ha="center", va="center",
                            fontsize=14 / len(site_coordinates)**0.25, 
                            zorder=100)

        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.tight_layout()

        display(fig)                  # show current frame in notebook
        clear_output(wait=True)       # remove previous frame output
        if pause > 0:
          time.sleep(pause)               # delay between frames
    #display(fig)                      # keep final frame visible    #plt.show()

        # if not show_axes: ax.set_axis_off()

        # if file_name is not None:
        #     plt.savefig(file_name)

        # if show:
        #     if pause == -1:
        #         plt.show()

        #     else:
        #         plt.pause(pause)

        #         if close:
        #             plt.close("all")


# ## Data and input files related functions

# In[ ]:


def get_dft_energy(path):
    """
    Extracts the DFT energy from a results file.

    Parameters:
    dir : Path
        Directory containing the 'short_results.txt' file.

    Returns:
    float: The DFT energy value from the last line of the results file.
        If there's only one line of data, returns that single energy value.
    """

    try:
        result = np.loadtxt(path / "short_results.txt", skiprows=1, unpack=True)
        return result[3] if len(result.shape)==1 else result[3][-1]
    except:
        return None


# In[ ]:


def replicate_structure(sites, sizex, sizey, Nrepx, Nrepy):
    """
    Replicates a structure by mapping site indices from the original structure to a larger repeated structure.

    Parameters:
    sites (list): List of site indices in the original structure.
    sizex (int): Number of columns in the original structure.
    sizey (int): Number of rows in the original structure.
    Nrepx (int): Number of times to replicate the structure in the x-direction.
    Nrepy (int): Number of times to replicate the structure in the y-direction.

    Returns:
    list: List of site indices in the replicated structure.
    """

    rep_sites = []
    for l in sites:
        j = (l-1) % sizey + 1
        i = (l-1) // sizey + 1
        for x in range(Nrepx):
            for y in range(Nrepy):
                new_l = j + y*sizey + Nrepy*sizey*(i-1+x*sizex)
                rep_sites.append(new_l)
    return rep_sites


# In[ ]:


def plot_numbered_cells(sites, num_cols, num_rows, rep_col, rep_row):
    """
    Visualizes a grid of cells with numbered positions, highlighting specified sites.

    Parameters:
    sites (list): List of site indices to highlight in the original structure.
    num_cols (int): Number of columns in the original structure.
    num_rows (int): Number of rows in the original structure.
    rep_col (int): Number of times to replicate the structure in the column direction.
    rep_row (int): Number of times to replicate the structure in the row direction.

    The function creates a plot where:
    - Red circles (o) mark the specified sites
    - Light blue squares mark empty sites
    - Each cell is numbered sequentially
    - The plot title shows the dimensions of the replicated structure
    """
    fig, ax = plt.subplots()
    new_sites = replicate_structure(sites, num_cols, num_rows, rep_col, rep_row)
    for i in range(rep_row*num_rows):  
        for j in range(rep_col*num_cols):
            s = j*rep_row*num_rows + i + 1
            if s in new_sites:
                ax.plot(j, i, 'o', color='red', markersize=20)
            else:
                ax.plot(j, i, 's', color='lightblue', markersize=40)  # Plot cells as squares
            ax.text(j, i, f'{s}', color='black', ha='center', va='center')  # Add numbers column-wise
    ax.set_xlim(-1, rep_col*num_cols)
    ax.set_ylim(-1, rep_row*num_rows)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(f"{rep_col*num_cols}x{rep_row*num_rows}")
    ax.axis('off')  # Remove coordinate axes
    plt.show()


def get_xy(lattice_input_file, idx=None):
    """
    Produces cartesian coordinates of adsorbates from history_output.txt
    Parameters
    ----------
    lattice_input_file : str or Path
        Path to the Zacros lattice input file.
    idx : int, optional : which configuration to extract (default is None, meaning all snapshots)
    Returns
    -------
    coords : np.array
        (List of) Array(s) of cartesian coordinates of adsorbates in the specified configuration.
    """

    # Read lattice input file

    try:
        with open(lattice_input_file, 'r') as f:
            content = [line for line in f.readlines() if (line.strip() and not line.startswith('#'))]
    except FileNotFoundError:
        raise FileNotFoundError(f"Lattice input file '{lattice_input_file}' not found.")
    finally:
        content = [line.split('#')[0] for line in content]

    for i,line in enumerate(content):
            if 'cell_vectors' in line:
                unit_cell = np.array([ [float(x) for x in content[i+1].split()],
                              [float(x) for x in content[i+2].split()] ])
            if 'repeat_cell' in line:
                repeat_cell = np.array([ int(x) for x in line.split()[1:3] ])
            if 'n_site_types' in line:
                n_site_types = int(line.split()[1])
            if 'n_cell_sites' in line:
                n_cell_sites = int(line.split()[1])
            if 'site_types' in line:
                site_types_names = line.split()[1:]
            if 'site_coordinates' in line:
                site_coordinates = []
                for j in range(n_cell_sites):
                    site_coordinates.append([float(x) for x in content[i+1+j].split()[:2]])
                site_uc_coordinates = np.array(site_coordinates) @ unit_cell


    site_coordinates = []
    site_types = []
    for i in range(repeat_cell[0]):
        for j in range(repeat_cell[1]):
            shift = np.array([i,j]) @ unit_cell
            for coords,st in zip(site_uc_coordinates , site_types_names):
                site_types.append(st)
                site_coordinates.append(coords + shift)
                site_types.append(st)

    # Read configurations from history_output.txt file

    try:
        with open(lattice_input_file.parent / 'history_output.txt', 'r') as f:
            content = f.readlines()
    except FileNotFoundError:
        print(f" File '{lattice_input_file.parent / 'history_output.txt'}' not found.")
    finally:

        # Get the species list
        species = content[1].split()[1:]
        # Set number of sites
        n_sites = repeat_cell[0] * repeat_cell[1] * n_cell_sites

        if idx is None:
          # Get number of configurations
          for line in content:
              if 'configuration' in line:
                  n_confs = int(line.split()[1])
          confs = np.zeros((n_confs, n_sites), dtype=int)
          for i_conf in range(n_confs):
            for i in range(n_sites):
              confs[i_conf, i] = int(content[7 + i_conf*(n_sites+1) +i].split()[2])
        else:
          confs = np.zeros((1, n_sites), dtype=int)
          for i in range(n_sites):
            confs[0, i] = int(content[7 + idx*(n_sites+1) +i].split()[2])

    # adsorbate coordinates
    ads_coords = []
    for conf in confs:
      ads_coords.append(np.array([(x,y) for (x, y), st, spec in zip(site_coordinates, site_types, conf) if (spec > 0)]))

    return ads_coords



# In[ ]:


def make_lattice_input(dir, Nrepx, Nrepy, simcell, Nrepx_lattice, Nrepy_lattice, header=""):
    """
    Creates a lattice input file for Zarcos.

    Parameters:
    ----------
    dir : Path
        Directory path where the lattice_input.dat file will be created
    Nrepx : int
        Number of repetitions in x direction for the unit cell
    Nrepy : int 
        Number of repetitions in y direction for the unit cell
    simcell : array_like
        Simulation cell vectors matrix
    Nrepx_lattice : int
        Number of times to repeat the cell in x direction in the lattice
    Nrepy_lattice : int
        Number of times to repeat the cell in y direction in the lattice
    header : str, optional
        Header text to add at the top of the file (default: "")

    Returns
    -------
    None
        Writes lattice_input.dat file to the specified directory
    """
    unitcell = simcell[:2][:,0:2] / np.array([Nrepx, Nrepy])

    cellvec1, cellvec2 = unitcell

    lattice_input_content = [
    f"# {header}\n",
    f"lattice periodic_cell\n",
    f"\n",
    f"cell_vectors       # in row format (Angstroms)\n",
    f"\n",
    f"   {cellvec1[0]:14.10f}   {cellvec1[1]:14.10f} \n",
    f"   {cellvec2[0]:14.10f}   {cellvec2[1]:14.10f}\n",
    f"\n",
    f"repeat_cell       {Nrepx_lattice:.0f} {Nrepy_lattice:.0f}\n",
    f"\n",
    f"n_site_types      1\n",
    f"site_type_names   fcc\n",
    f"\n",
    f"n_cell_sites      1\n",
    f"site_types        fcc\n",
    f"\n",
    f"site_coordinates   # fractional coordinates (x,y) in row format\n",
    f"\n",
    f"   0.333333333333333   0.333333333333333\n",
    f"\n",
    f"neighboring_structure\n",
    f"   \n",
    f"   1-1  north\n",
    f"   1-1  east\n",
    f"   1-1  southeast\n",
    f"\n",
    f"end_neighboring_structure\n",
    f"\n",
    f"end_lattice\n",
    ]

    with open(dir / "lattice_input.dat", "w") as file:
        for line in lattice_input_content:
            file.write(line)

    return



# In[ ]:


def make_state_input(dir, adsorbates, random=True, header=""):
    """
    Creates a state input file for Zarcos.

    Parameters:
    ----------
    dir : Path
        Directory path where the state_input.dat file will be created
    adsorbates : list of tuples
        List of (species, amount) tuples where species is a string and amount is an integer
    header : str, optional
        Header text to add at the top of the file (default: "")

    Returns
    -------
    None
        Writes state_input.dat file to the specified directory
    """

    state_input_content = [
    f"# {header}\n",
    f"\n",
    f"initial_state\n",
    f"\n",
    ]

    if random:
        for s in adsorbates:
            line =  [
                        f"seed_multiple {s[0]} {s[1] }\n",
                        f"  site_types fcc\n",
                        f"end_seed_multiple\n",
                    ]
            state_input_content.extend(line)

    else:
        ads_list = []
        for s in adsorbates:
            ads_list.extend([(s[0],i+1) for i in range(len(ads_list),len(ads_list)+s[1])])
        line = [f"seed_on_sites {a[0]} {a[1]:3.0f}\n" for a in ads_list]
        state_input_content.extend(line)

    state_input_content.append(f"\n")
    state_input_content.append(f"end_initial_state\n")

    with open(dir / "state_input.dat", "w") as file:
        for line in state_input_content:
            file.write(line)

    return



# In[ ]:


def make_mechanism_input(dir, header=""):

    lines = [
    f"# {header}\n",
    f"\n",
    f"mechanism\n",
    f"\n",
    f"reversible_step O_hopping\n",
    f"  sites 2\n",
    f"  neighboring 1-2\n",
    f"  initial # (entitynumber, species, dentate)\n",
    f"    1 O*    1\n",
    f"	  2 *     1\n",
    f"  final\n",
    f"    2 *     1\n",
    f"	  1 O*    1\n",
    f"  site_types fcc fcc\n",
    f"  pre_expon   5.0e9\n",
    f"  pe_ratio    1.0\n",
    f"  activ_eng   0.43\n",
    f"  prox_factor 0.5\n",
    f"end_reversible_step\n",
    f"\n",
    f"end_mechanism\n",
    f"\n"
      ]

    with open(dir / "mechanism_input.dat", "w") as file:
        file.writelines(lines)

    return



# In[ ]:


def make_energy(conf_dir, energy):
    """
    Creates an energy file for Zarcos containing the formation energy

    Parameters:
    ----------
    conf_dir : Path
        Directory where to write the energy file
    energy : float
        Formation energy to be written to the file
    """
    # Write energy to file
    with open(conf_dir / "energy", "w") as f:
        f.write(f"{energy:.6f}\n")
    return



# In[ ]:


def make_calculation_input(dir, number_of_confs):

    with open(dir / "calculation_input.dat", "w") as f:
        f.write(f"n_config {number_of_confs}\n")
        f.write(f"\n")
        f.write(f"n_surf_species            1\n")
        f.write(f"surf_specs_names        O*\n")
        f.write(f"surf_specs_dent         1\n")
        f.write(f"\n")
        f.write(f"n_site_types              1\n")
        f.write(f"site_type_names           fcc\n")
        f.write(f"\n")
        f.write(f"debug_report_global_energetics\n")
        f.write(f"debug_check_lattice\n")
        f.write(f"\n")
        f.write(f"finish\n")
    return



# In[ ]:


def make_simulation_input(dir, temperature=300, snapshots='on event 1', max_steps=100, wall_time=60, header=""):

    """ 
    Creates a simulation input file for Zarcos.
    with parameters for temperature and header.
    """

    with open(dir / "simulation_input.dat", "w") as f:
        f.write(f"# {header}\n")
        f.write(f"\n")
        f.write(f"random_seed               314159265")
        f.write(f"\n")
        f.write(f"temperature               {float(temperature)}\n")
        f.write(f"pressure                  1.00\n")
        f.write(f"\n")
        f.write(f"n_gas_species             0\n")
        f.write(f"\n")
        f.write(f"n_surf_species            1\n")
        f.write(f"surf_specs_names          O*\n")
        f.write(f"surf_specs_dent           1\n")
        f.write(f"kmc_propagation_method first_reaction binary_heap\n")
        f.write(f"\n")
        f.write(f"override_array_bounds & & 100 145\n")
        f.write(f"\n")
        f.write(f"snapshots                 {snapshots}\n")
        f.write(f"process_statistics        off #on realtime 2e-1\n")
        f.write(f"species_numbers           off #on realtime 2e-1\n")
        f.write(f"energetics_lists          off #on realtime 2e-1\n")
        f.write(f"process_lists             off #on realtime 2e-1\n")
        f.write(f"\n")
        f.write(f"event_report              off\n")
        f.write(f"on_sites_seeding_report   off\n")
        f.write(f"\n")
        f.write(f"max_steps                 {max_steps}\n")
        f.write(f"# max_time                  1e-3\n")
        f.write(f"\n")
        f.write(f"wall_time                 {wall_time} # in seconds\n")
        f.write(f"\n")
        f.write(f"no_restart\n")
        f.write(f"\n")
        f.write(f"# debug_report_processes\n")
        f.write(f"# debug_report_global_energetics\n")
        f.write(f"# debug_check_processes\n")
        f.write(f"# debug_check_lattice\n")
        f.write(f"\n")
        f.write(f"finish\n")

    return


# In[ ]:


def get_dft_data(dft_dir):
    """
    Retrieves or processes DFT calculation data from a directory.

    This function first attempts to load existing data from 'dft_data.csv'.
    If the file doesn't exist, it processes DFT data from scratch by scanning
    the directory for 'short_results.txt' files and extracting relevant information.

    Parameters
    ----------
    dft_dir : Path
        Directory containing DFT calculation results in a hierarchical structure.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing DFT calculation data with columns:
        - path: Path to the results file
        - system: System name (e.g., 'Pt(111)')
        - cell_size: Size of the cell (e.g., '2x2')
        - kpoints: k-points setting used
        - functional: DFT functional used
        - type: Calculation type (slab, bulk, O_ads, etc.)
        - adsorbates: List of adsorbate site indices (None for slab/bulk)
        - energy: DFT energy value

    Notes
    -----
    The function caches results in 'dft_data.csv' for faster subsequent access.
    """
    # If file doesn't exist, process the data from scratch
    paths_to_dft_data = Path(dft_dir).rglob("short_results.txt")
    paths_to_dft_data = [p for p in paths_to_dft_data if "ads" in str(p)]
    df = pd.DataFrame([p for p in paths_to_dft_data], columns=['path'])

    # Extract information from paths
    df['system'] =      df['path'].apply(lambda x: x.parents[4].name.split("_")[0])
    df['cell_size'] =   df['path'].apply(lambda x: x.parents[4].name.split("_")[1])
    df['kpoints'] =     df['path'].apply(lambda x: x.parents[3].name)
    df['functional'] =  df['path'].apply(lambda x: x.parents[2].name)
    df['adsorbates'] =  df['path'].apply(lambda x: [int(i) for i in x.parents[0].name.split("_")])
    df['energy'] =      df['path'].apply(lambda x: get_dft_energy(x.parent))
    df['slab_energy'] = df['path'].apply(lambda x: get_dft_energy(x.parents[2] / "slab"))

    return df


# In[ ]:


def get_cefit_output(wdir):
    """
    Loads and parses the cefit_output.txt file.

    Parameters
    ----------
    wdir : Path or str
        Directory containing the cefit_output.txt file

    Returns
    -------
    cl_names, cl_energies, cl_multiplicity
    """
    try:
        with open(Path(wdir) / 'cefit_output.txt', 'r') as f:
            lines = f.readlines()

        cl_names = lines[0].strip().split()
        cl_energies = [float(line.split()[0]) for line in lines[1:]]
        cl_multiplicity = [int(line.split()[1]) for line in lines[1:]]

        return cl_names, cl_energies, cl_multiplicity

    except FileNotFoundError:
        print(f"Could not find cefit_output.txt in {wdir}")
        return None


# We calculate the adsorbate interaction energy for configurations prepared for Zacros as follows:
# 
# $\qquad V = n_xn_y \left( E_\text{DFT} - E_\text{slab} - N_\text{ads}E_\text{ref}\right),$
# 
# where $E_\text{DFT}$ is the energy of the configuration produced by a DFT code, 
# 
# $E_\text{slab}$ is the DFT energy of the configuration without adsorbates,  
# 
# $ E_\text{ref} = E^{(8\times8)}_\text{DFT} - E^{(8\times8)}_\text{slab}$
# is the reference energy defined as adsorption energy of a single O at $(8\times8)$ slab (Florian Nitz's DFT calculations show that this slab is large enough to make image and lateral interaction negligible), 
# 
# $N_\text{ads}$ is the number of adsorbates in the configuration, 
# 
# $n_x$ and $n_y$ are the number of repetitions of the DFT configuration in $x$ and $y$ directions, resp.
# 
# For example, DFT calculations for an singe O-atom adsorbed at fcc site of the $(3\times3)$ slab give
# $ E_\text{DFT} = -192.4571\,\text{eV}$ and $ E_\text{slab} = -186.7272\,\text{eV}$, so that
# 
# $\qquad E_\text{DFT} - E_\text{slab} = -5.7299 \,\text{eV},$
# 
# and the reference energy defined above is $ E_\text{ref} = -5.8092\,\text{eV}$. Then ($N_\text{ads}=1$)
# 
# $\qquad E_\text{DFT} - E_\text{slab} - E_\text{ref} = 0.0793 \,\text{eV}$
# 
# is  the lateral interaction energy due to the small size of the $(3\times3)$ simulation cell. Since the closest O-atom that the adsorbate sees via periodic boundary conditions belongs to the 5th shell and the total number of 5th-shell neighbors is 6, the lateral interaction energy per pair of O atoms at that distance is approximately
# 
# $\qquad 0.0793/6 = 0.0132  \,\text{eV}$
# 
# **Question**: is factor $n_xn_y$ used in the 1st formula correct from point of view of counting the O-O pairs

# In[ ]:


0.0793/6


# In[ ]:


def produce_fit(wdir, dfz, Nrepx_target, Nrepy_target, shell_list, Eref, show=False):
    #
    # Produce energetics and calculation input files
    #
    make_energetics_input(Path(wdir), shell_list)
    make_calculation_input(Path(wdir), len(dfz))

    # Loop over the configurations
    conf_counter = 1
    E_form = []
    n_ads = []
    for idx, row in dfz.iterrows():

        # Get cell size as a list of integers
        Nrepx, Nrepy = [ int(n) for n in row['cell_size'].split('x')]

        # Check if the cell size can be replicated to the target lattice size
        if (Nrepx_target % Nrepx != 0) or (Nrepy_target % Nrepy) != 0:
            print(f"WARNING! {row['cell_size']} cell cannot be replicated to {Nrepx_target}x{Nrepy_target} lattice. Skipping this configuration.")
            continue
        # Get the tagret-origin repetitions ratio in x and y directions
        Nrepx_ratio = Nrepx_target // Nrepx
        Nrepy_ratio = Nrepy_target // Nrepy

        # Read VASP configuration files
        sini    = ase.io.read(row['path'].parent / "POSCAR")
        sopt    = ase.io.read(row['path'].parent / "CONTCAR")

        # Get the adsorbate indices
        O_ind   = np.where(np.array(sopt.get_chemical_symbols()) == "O")[0]

        # skip configuration if adsorbate atoms moved by more than 0.5 Angstrom during geometry optimization
        delta_r = np.sum((sopt.positions[O_ind] - sini.positions[O_ind])**2, axis=1)**0.5
        if np.any(delta_r > 0.5):
            print(f"WARNING! Large displacement in folder {row['path'].parent.name} for adsorbate atoms. Skipping this configuration.")
            continue


        # Get configuration directory name and create it
        conf_dir = Path(wdir) / f"Conf{conf_counter}"
        conf_dir.mkdir(exist_ok=True)

        # Get the adsorbate Zacros positions
        ads_pos = [ ('O*',pos) for pos in replicate_structure(np.array(row['adsorbates']) + 1, Nrepx, Nrepy, Nrepx_ratio, Nrepy_ratio) ]

        # Calculate formation energy per adsorbate from DFT data accounting for the replication
        E_form.append(Nrepx_ratio*Nrepy_ratio*(row['energy'] - row['slab_energy'] - Eref*len(row['adsorbates'])))
        n_ads.append(Nrepx_ratio*Nrepy_ratio*len(row['adsorbates']))

        # Write the Zacros input files for ce_fit
        make_lattice_input(conf_dir, Nrepx, Nrepy, sini.cell, Nrepx_target, Nrepy_target, header=f"{Nrepx_target}x{Nrepy_target} Pt(111)")
        make_state_input(conf_dir, ads_pos)
        make_energy(conf_dir, E_form[-1])

        # Print the configuration directory name
        print(f"Configuration {conf_counter} created in {conf_dir}")
        # Show the configuration
        if show: plot_numbered_cells(np.array(row['adsorbates']) + 1, Nrepx, Nrepy, 1, 1)
        # Increment the configuration counter
        conf_counter += 1

    #
    # Do CE fit
    #

    # Get the full path to ce_fit.x
    ce_fit_path = zacros_path / "ce_fit.x"

    # Change to the zacros run directory where input files are located
    original_dir = os.getcwd()
    os.chdir(wdir)

    try:
        # Run ce_fit.x
        result = subprocess.run([str(ce_fit_path)], 
                            check=True,
                            capture_output=True,
                            text=True)
        print(result.stdout)

    except subprocess.CalledProcessError as e:
        print("Error running ce_fit.x:")
        print(e.stderr)

    finally:
        # Change back to original directory
        os.chdir(original_dir)

    #
    # Get the fit results
    #

    # Read the output file
    with open(wdir / 'general_output.txt', 'r') as f:
        content = f.readlines()

    # Check the termination status
    if len(content) > 0:
        if  ('> Normal termination <' not in content[-1]):
            for i,line in enumerate(content):
                if  'error' in line.lower():
                    print(f"Error message in general_output.txt:\n")
                    [print("   " + l) for l in content[i:]]
                    break
    else:
        print("Output file is empty.")

    E_fit = [float(line.strip().split()[-1]) for line in content if 'Total adlayer energy' in line][len(E_form):]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Fit quality check')

    #
    # Dispay quality of the fit
    #
    ax1.scatter(E_form, E_fit, color='black', label='Conf Energy')
    ax1.plot([min(E_form), max(E_form)], [min(E_form), max(E_form)], 'r--', label=r'$E_{{fit}}$ = $E_{{DFT}}$')
    ax1.legend()
    ax1.set_ylabel(r'$E_{{fit}}$ (eV)')
    ax1.set_xlabel(r'$E_{{DFT}}$ (eV)')
    ax1.set_ylabel(r'$E_{{fit}}$ (eV)')

    x = range(1, conf_counter)
    ax2.scatter(x, [e/n for e,n in zip(E_form,n_ads)], color='black', marker='s', facecolors='none', label='DFT')
    ax2.scatter(x, [e/n for e,n in zip(E_fit,n_ads)], color='red', facecolors='none',label='CE fit')
    ax2.legend()
    ax2.set_xlabel('Configuration number')
    ax2.set_ylabel('Formation Energy per ads (eV)')

    plt.tight_layout()
    plt.savefig(wdir / 'fit_quality.pdf', bbox_inches='tight')
    plt.show()

    return E_form, E_fit   


# ## Cluster definitions

# In[ ]:


#   Generates Zacros cluster definition for a single site with one adsorbate.
def cluster_1_site():
  content = [
    f"cluster O_fcc\n",
    f"sites 1\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  site_types fcc\n",
    f"  graph_multiplicity 1\n",
    f"end_cluster\n"]

  return content

#   Generates Zacros 2-site cluster with 2 adsorbates.
def cluster_2_site():
  content = [
    f"cluster O_fcc-nn1\n",
    f"sites 2\n",
    f"neighboring 1-2\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"end_cluster\n"]

  return content

#   Generates Zacros 3-site cluster with 2 adsorbates.
def cluster_3_site_2nn():
  content = [
    f"cluster O_fcc-nn2\n",
    f"sites 3\n",
    f"neighboring 1-2 2-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  & &  &\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:-120.0\n",
    f"end_cluster\n"]

  return content

def cluster_3_site_3nn():
  content = [
    f"cluster O_fcc-nn3\n",
    f"sites 3\n",
    f"neighboring 1-2 2-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  & &  &\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:180.0\n",
    f"end_cluster\n"]

  return content

#   Generates Zacros 4-site clusters with 2 adsorbates.
def cluster_4_site_4nn():
  content = [
    f"cluster O_fcc-nn4\n",
    f"sites 4\n",
    f"neighboring 1-2 2-3 3-4\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:180.0 2-3-4:-120.0\n",
    f"end_cluster\n"]

  return content

def cluster_4_site_5nn():
  content = [
    f"cluster O_fcc-nn5\n",
    f"sites 4\n",
    f"neighboring 1-2 2-3 3-4\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:180.0 2-3-4:180.0\n",
    f"end_cluster\n"]

  return content

#   Generates Zacros 5-site clusters with 2 adsorbates.
def cluster_5_site_6nn():
  content = [
    f"cluster O_fcc-nn6\n",
    f"sites 5\n",
    f"neighboring 1-2 2-3 3-4 4-5\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:180.0 2-3-4:-120.0  3-4-5:180.0\n",
    f"end_cluster\n"]

  return content

def cluster_5_site_7nn():
  content = [
    f"cluster O_fcc-nn7\n",
    f"sites 5\n",
    f"neighboring 1-2 2-3 3-4 4-5\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:180.0 2-3-4:180.0  3-4-5:-120.0\n",
    f"end_cluster\n"]

  return content

def cluster_5_site_8nn():
  content = [
    f"cluster O_fcc-nn8\n",
    f"sites 5\n",
    f"neighboring 1-2 2-3 3-4 4-5\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:180.0 2-3-4:180.0  3-4-5:180.0\n",
    f"end_cluster\n"]

  return content

def cluster_6_site_9nn():
  content = [
    f"cluster O_fcc-nn9\n",
    f"sites 6\n",
    f"neighboring 1-2 2-3 3-4 4-5 5-6\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  & &  &\n",
    f"  2 O* 1\n",
    f"  site_types fcc fcc fcc fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:180.0 2-3-4:180.0  3-4-5:-120.0 4-5-6:180.0\n",
    f"end_cluster\n"]

  return content

#   Generates Zacros 3-site cluster with 3 adsorbates.
def cluster_3_site_3():
  content = [
    f"cluster O_fcc-3-3\n",
    f"sites 3\n",
    f"neighboring 1-2 2-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 3\n",
    f"  angles 1-2-3:180.0\n",
    f"end_cluster\n"]

  return content


# In[ ]:


def make_energetics_input(dir, cluster_list, eng_list=None):
    """
    Creates an energetics input file for Zacros containing cluster definitions.

    Parameters
    ----------
    dir : Path
        Directory where to write the energetics_input.dat file
    cluster_list : list
        List of integers indicating which clusters to include in the input file.
        The integers correspond to:
        0: single site cluster
        1: 2-site cluster with 1st nearest neighbor
        2: 3-site cluster with 2nd nearest neighbor
        3: 3-site cluster with 3rd nearest neighbor
        4: 4-site cluster with 4th nearest neighbor
        5: 4-site cluster with 5th nearest neighbor
        6: 5-site cluster with 6th nearest neighbor
        7: 5-site cluster with 7th nearest neighbor
        8: 5-site cluster with 8th nearest neighbor
        9: 6-site cluster with 9th nearest neighbor
       10: 3-site cluster with 3 adsorbates
    eng_list : list, optional
        List of floats indicating the cluster energies corresponding to the clusters in shell_list.
        If None, no energies are included in the input file (default: None)

    Returns
    -------
    None
        Writes energetics_input.dat file to the specified directory

    """

    dispatcher = {0: cluster_1_site,
                  1: cluster_2_site,
                  2: cluster_3_site_2nn,
                  3: cluster_3_site_3nn,
                  4: cluster_4_site_4nn,
                  5: cluster_4_site_5nn,
                  6: cluster_5_site_6nn,
                  7: cluster_5_site_7nn,
                  8: cluster_5_site_8nn,
                  9: cluster_6_site_9nn,
                 10: cluster_3_site_3
                 }

    with open(dir / "energetics_input.dat", "w") as f:
        f.write('# O at Pt(111)\n')
        f.write('# For structures and values see Dropbox:\n')
        f.write('# "Kinetics of Surface Reactions/zacros/O_Pt111/O_Pt111 structures.pptx"\n')
        f.write('\n')
        f.write('energetics\n')
        f.write('\n')

        for i, s in enumerate(cluster_list):
            content = dispatcher[s]()
            if eng_list is not None:
                content.insert(-1, f"  cluster_eng   {eng_list[i]:.6f}\n")
            [f.write(line) for line in content]
            f.write('\n')

        f.write('end_energetics\n')


    return

