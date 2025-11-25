#!/usr/bin/env python
# coding: utf-8

# Functions to generate Zacros cluster definitions.

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
def cluster_1_1_2():
  content = [
    f"cluster O_fcc-1-1-2\n",
    f"sites 3\n",
    f"neighboring 1-2 2-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:-120.0\n",
    f"end_cluster\n"]

  return content

def cluster_1_1_3():
  content = [
    f"cluster O_fcc-1-1-3\n",
    f"sites 3\n",
    f"neighboring 1-2 2-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 2\n",
    f"  angles 1-2-3:180.0\n",
    f"end_cluster\n"]

  return content

def cluster_1_1_1a1():
  content = [
    f"cluster O_fcc-1-1-1a1\n",
    f"sites 3\n",
    f"neighboring 1-2 1-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 1\n",
    f"  angles 1-2-3:-60.0\n",
    f"  no_mirror_images\n",
    f"  absl_orientation 1-2:0.0\n",
    f"end_cluster\n"]

  return content

def cluster_1_1_1a2():
  content = [
    f"cluster O_fcc-1-1-1a2\n",
    f"sites 3\n",
    f"neighboring 1-2 1-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 1\n",
    f"  angles 1-2-3:-60.0\n",
    f"  no_mirror_images\n",
    f"  absl_orientation 1-2:180.0\n",
    f"end_cluster\n"]

  return content

def cluster_1_1_1a3():
  content = [
    f"cluster O_fcc-1-1-1a3\n",
    f"sites 3\n",
    f"neighboring 1-2 1-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 1\n",
    f"  angles 1-2-3:-60.0\n",
    f"  no_mirror_images\n",
    f"  absl_orientation 1-2:120.0\n",
    f"end_cluster\n"]

  return content

def cluster_1_1_1a4():
  content = [
    f"cluster O_fcc-1-1-1a4\n",
    f"sites 3\n",
    f"neighboring 1-2 1-3\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc\n",
    f"  graph_multiplicity 1\n",
    f"  angles 1-2-3:-60.0\n",
    f"  no_mirror_images\n",
    f"  absl_orientation 1-2:-120.0\n",
    f"end_cluster\n"]

  return content

def cluster_1_2_3a1():
  content = [
    f"cluster O_fcc-1-2-3a1\n",
    f"sites 4\n",
    f"neighboring 1-2 1-3 2-3 3-4\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  & &  &\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc fcc\n",
    f"  graph_multiplicity 1\n",
    f"  angles 1-2-3:-60.0  2-3-4:180.0\n",
    f"  no_mirror_images\n",
    f"  absl_orientation 1-2:60.0\n",
    f"end_cluster\n"]

  return content

def cluster_1_2_3a2():
  content = [
    f"cluster O_fcc-1-2-3a2\n",
    f"sites 4\n",
    f"neighboring 1-2 1-3 2-3 3-4\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  & &  &\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc fcc\n",
    f"  graph_multiplicity 1\n",
    f"  angles 1-2-3:-60.0  2-3-4:120.0\n",
    f"  no_mirror_images\n",
    f"  absl_orientation 1-2:60.0\n",
    f"end_cluster\n"]

  return content

def cluster_1_2_3a3():
  content = [
    f"cluster O_fcc-1-2-3a3\n",
    f"sites 4\n",
    f"neighboring 1-2 1-3 2-3 3-4\n",
    f"lattice_state\n",
    f"  1 O* 1\n",
    f"  2 O* 1\n",
    f"  & &  &\n",
    f"  3 O* 1\n",
    f"  site_types fcc fcc fcc fcc\n",
    f"  graph_multiplicity 1\n",
    f"  angles 1-2-3:-60.0  2-3-4:120.0\n",
    f"  no_mirror_images\n",
    f"  absl_orientation 1-2:60.0\n",
    f"end_cluster\n"]

  return content


cluster O_triplet_1-2-3b
  sites 4
  neighboring 1-2 1-3 2-3 3-4
  lattice_state
    1 O*   1
    2 O*   1
    & &    &
    3 O*   1
	
  variant b1
    site_types          fcc_t fcc_t fcc_t fcc_t
    graph_multiplicity  1
    angles              1-2-3:60.0  2-3-4:180.0 
    no_mirror_images
    absl_orientation    1-2:180.0
    cluster_eng         0.016
  end_variant

  variant b2
    site_types          fcc_t fcc_t fcc_t fcc_t
	graph_multiplicity  1
    angles              1-2-3:-60.0  2-3-4:180.0 
    no_mirror_images
    absl_orientation    1-2:0.0
    cluster_eng         0.016
  end_variant

  variant b3
    site_types          fcc_t fcc_t fcc_t fcc_t
    graph_multiplicity  1
    angles              1-2-3:60.0  2-3-4:180.0 
    no_mirror_images
    absl_orientation    1-2:-60.0
    cluster_eng         0.016
  end_variant

  variant b4
    site_types          fcc_t fcc_t fcc_t fcc_t
    graph_multiplicity  1
    angles              1-2-3:-60.0  2-3-4:180.0 
    no_mirror_images
    absl_orientation    1-2:-120.0
    cluster_eng         0.016
  end_variant

  variant b5
    site_types          fcc_t fcc_t fcc_t fcc_t
    graph_multiplicity  1
    angles              1-2-3:60.0  2-3-4:180.0 
    no_mirror_images
    absl_orientation    1-2:60.0
    cluster_eng         0.016
  end_variant
  
  variant b6
    site_types          fcc_t fcc_t fcc_t fcc_t
    graph_multiplicity  1
    angles              1-2-3:-60.0  2-3-4:180.0 
    no_mirror_images
    absl_orientation    1-2:120.0
    cluster_eng         0.016
  end_variant

end_cluster


############################################################################

end_energetics

