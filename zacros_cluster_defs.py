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

